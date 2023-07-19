// gcc -fopenmp -Wall -m64 -I"${MKLROOT}/include" omp_zpotrf.c -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group ${CUDA_LINK_FLAGS} -liomp -lpthread -lm -ldl
//#include <math.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>

//#include <time.h>

#include "zpotrf_omp.h"

#define MIN(a,b) (a < b ? a : b)
#define ZPOTRF_BLOCKSIZE (128)

/*
  own version of zpotrf using omp parallelization
  PRACE_Advanced_OpenMP_Tutorial.pdf && https://fdesprez.github.io/teaching/par-comput/lectures/slides/L4-OpenMP-pt2-2p.pdf
*/
int zpotrf_omp(const char uplo, const int n,
               double complex* restrict const a, const int lda, int* shared_info) {
	#pragma omp single
	{
	double complex cone = 1.;
	double complex cmone = -1.;
	
	// test input parameters
	const int upper = (uplo == 'U' || uplo == 'u');
	*shared_info = 0;
	if(!upper && !(uplo == 'L' || uplo == 'l')) {
		*shared_info = -1;
	} else if( n < 0) {
		*shared_info = -2;
	} else if(lda < (n>1?n:1)) { // lda < max(1,n)
		*shared_info = -4;
	}
	
	// blocksize 128 seems to work best, followed by 64
	const int nb = ZPOTRF_BLOCKSIZE; // ILAENV(&ione, "zpotrf", &uplo, &n, &imone, &imone, &imone);
	int info;
	if(*shared_info != 0 || n == 0)
	{
		// nothing, early return
	} else if(nb <= 1 || nb >= n) {
		// Use unblocked code
		//info = LAPACKE_zpotrf(LAPACK_COL_MAJOR, uplo, n, (lapack_complex_double*) a, lda);
		zpotf2_(&uplo, &n, (lapack_complex_double*) a, &lda, shared_info);
	} else {
		// Use blocked code
		if(upper) {
			for(int j = 0; j < n; j += nb) {
				int jb = MIN(nb, n-j);
				// Diagonal Block factorization
				#pragma omp task depend(inout: a[j + j*lda]) //firstprivate(jb, info)
				{
					if(!*shared_info) {
						zpotf2_(&uplo, &jb, (lapack_complex_double*) &a[j + j*lda], &lda, &info);
						//info = LAPACKE_zpotrf(LAPACK_COL_MAJOR, 'U', jb,
						//                      (lapack_complex_double*) &a[j + j*lda], lda);
						if(info) {
							*shared_info = info;
						}
					}
				}
				
				// Triangular systems
				for (int i = j+nb; i < n; i += nb) {
					const int ib = MIN(nb, n-i);
					#pragma omp task depend(in: a[j + j*lda])\
					                 depend(inout: a[j + i*lda]) //firstprivate(ib, jb)
					{
						if(!*shared_info)
							cblas_ztrsm(CblasColMajor, CblasLeft, CblasUpper, CblasConjTrans,
							            CblasNonUnit, jb, ib, &cone,
							            (lapack_complex_double*)  &a[j + j*lda],
							            lda, &a[j + i*lda], lda);
					}
				}
				
				// Update trailing matrix
				for (int i = j+nb; i < n; i += nb) {
					const int ib = MIN(nb, n-i);
					for (int k = j+nb; k < i; k += nb) {
						const int kb = MIN(nb, n-k);
						#pragma omp task depend(in: a[j + k*lda], a[j + i*lda])\
						                 depend(inout: a[k + i*lda]) //firstprivate(ib, jb, kb)
						{
							if(!*shared_info)
								cblas_zgemm (CblasColMajor, CblasConjTrans, CblasNoTrans,
								             kb, ib, jb,
								             &cmone, (lapack_complex_double*) &a[j + k*lda], lda,
								                     (lapack_complex_double*) &a[j + i*lda], lda,
								             &cone,  (lapack_complex_double*) &a[k + i*lda], lda);
						}
					}
					#pragma omp task depend(in: a[j + i*lda])\
					                 depend(inout: a[i + i*lda]) //firstprivate(ib, jb)
					{
						if(!*shared_info)
							cblas_zherk(CblasColMajor, CblasUpper, CblasConjTrans,
							            ib, jb, //const MKL_INT k,
							            -1., (lapack_complex_double*) &a[j + i*lda], lda,
							             1., (lapack_complex_double*) &a[i + i*lda], lda);
					}
				}
			}
		} else {
			// Compute the Cholesky factorization A = L*L'
			for(int j = 0; j < n; j += nb) {
				int jb = MIN(nb, n-j);
				// Diagonal Block factorization
				#pragma omp task depend(inout: a[j + j*lda]) //firstprivate(jb, info)
				{
					if(!*shared_info) {
						zpotf2_(&uplo, &jb, (lapack_complex_double*) &a[j + j*lda], &lda, &info);
						//info = LAPACKE_zpotrf(LAPACK_COL_MAJOR, 'L', jb,
						//                      (lapack_complex_double*) &a[j + j*lda], lda);
						if(!info) {
							*shared_info = info;
						}
					}
				}
				
				// Triangular systems
				for (int i = j+nb; i < n; i += nb) {
					const int ib = MIN(nb, n-i);
					#pragma omp task depend(in: a[j + j*lda])\
					                 depend(inout: a[i + j*lda]) //firstprivate(ib, jb)
					{
						if(!*shared_info)
							cblas_ztrsm(CblasColMajor, CblasRight, CblasLower, CblasConjTrans,
							            CblasNonUnit, ib, jb, &cone, &a[j + j*lda],
							            lda, &a[i + j*lda], lda);
					}
				}
				
				// Update trailing matrix
				for (int i = j+nb; i < n; i += nb) {
					const int ib = MIN(nb, n-i);
					for (int k = j+nb; k < i; k += nb) {
						const int kb = MIN(nb, n-k);
						#pragma omp task depend(in: a[i + j*lda], a[k + j*lda])\
						                 depend(inout: a[i + k*lda]) //firstprivate(ib, jb, kb)
						{
							if(!*shared_info)
								cblas_zgemm (CblasColMajor, CblasNoTrans, CblasConjTrans,
								             ib, kb, jb,
								             &cmone, &a[i + j*lda], lda,
								                     &a[k + j*lda], lda,
								             &cone,  &a[i + k*lda], lda);
						}
					}
					#pragma omp task depend(in: a[i + j*lda])\
					                 depend(inout: a[i + i*lda]) //firstprivate(ib, jb)
					{
						if(!*shared_info)
							cblas_zherk(CblasColMajor, CblasLower, CblasNoTrans,
							            ib, jb, //const MKL_INT k,
							            -1., &a[i + j*lda], lda,
							             1., &a[i + i*lda], lda);
					}
				}
			}
		}
	}
	} // pragma omp single
	// implicit taskwait
	return *shared_info;
}

/*
void printMatrix(int m, int n, const double complex*A, int lda, const char* praetext)
{
	printf("%s\n", praetext);
	for(int row = 0 ; row < m ; row++){
		for(int col = 0 ; col < n ; col++){
			double complex Areg = A[row + col*lda];
			printf("%g%+gi\t", creal(Areg), cimag(Areg));
		}
		printf("\n");
	}
}

int main(int argc, char* argv[]) {
	int n = atoi(argv[1]);
	
	struct timespec time_start, time_end;
	int print_matrix_in  = 0;
	int print_matrix_out = 0;
	complex cone = 1.;
	complex czero = 0.;
	complex double* a = NULL; // = (complex double*) malloc(n*n*sizeof(double complex));
	complex double* b = NULL; // = (complex double*) malloc(n*n*sizeof(double complex));
	
	int err =
	posix_memalign((void **) &a, 64, n*n*sizeof(double complex));
	if(err) {
		fprintf(stderr, "malloc failed\n");
	}
	err =
	posix_memalign((void **) &b, 64, n*n*sizeof(double complex));
	if(err) {
		fprintf(stderr, "malloc failed\n");
	}
	
	const char uplo = 'l';
	
	srand(time(NULL));
	
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			a[i + j*n] = rand()/(double) RAND_MAX + I*rand()/(double) RAND_MAX;
		}
	}
	cblas_zgemm (CblasColMajor, CblasConjTrans, CblasNoTrans, n, n, n,
	             &cone, a, n, a, n, &czero, b, n);
	
	memcpy(a, b, sizeof(complex double)*n*n);
	
	if(print_matrix_in)
		printMatrix(n, n, a, n, "a=");
	
	clock_gettime(CLOCK_MONOTONIC_RAW, &time_start);
	int info =
	LAPACKE_zpotrf(LAPACK_COL_MAJOR, uplo, n, (lapack_complex_double*) a, n);
	
	clock_gettime(CLOCK_MONOTONIC_RAW, &time_end);
	double runtime = (double) (time_end.tv_sec - time_start.tv_sec) + (time_end.tv_nsec - time_start.tv_nsec) / 1e9;
	printf("runtime lapack_zpotrf: %g\n", runtime);
	printf("info=%i\n", info);
	if(print_matrix_out)
		printMatrix(n, n, a, n, "lapack_zpotrf(a)=");
	
	if(print_matrix_in)
		printMatrix(n, n, b, n, "\nb=");
	
	
	#pragma omp parallel num_threads(4)
	{
		#pragma omp single
			clock_gettime(CLOCK_MONOTONIC_RAW, &time_start);
		
		
		own_zpotrf(uplo, n, b, n, &info);
		
		#pragma omp single
			clock_gettime(CLOCK_MONOTONIC_RAW, &time_end);
	}
	
	runtime = (double) (time_end.tv_sec - time_start.tv_sec) + (time_end.tv_nsec - time_start.tv_nsec) / 1e9;
	printf("runtime own_zpotrf: %g\n", runtime);
	
	printf("info=%i\n", info);
	
	if(print_matrix_out)
		printMatrix(n, n, b, n, "own_zpotrf(b=a)=");
	
	
	int identical = 1;
	double max_val = 0.;
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			if(fabs(a[i + j*n] - b[i+j*n]) > 1e-10) {
				printf("difference in %i, %i: %g\n", i, j, fabs(a[i + j*n] - b[i+j*n]));
				identical = 0;
			}
			if(fabs(b[i+j*n]) > max_val) {
				max_val = fabs(b[i+j*n]);
			}
		}
	}
	if(identical) {
		printf("\nidentical results\n");
	}
	
	printf("max_val = %g\n", max_val);
	
	
	free(a);
	free(b);
	return 0;
}
*/
