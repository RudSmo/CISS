#include "sim_MCE_util.h"

#include <stdio.h>
#include <string.h>

#include <mkl.h>

/*
prints the overlap matrix. the hermitean property is assumed and only the
upper triangle (CblasColMajor) storage referenced
*/
void MCE1_print_overlap_matrix(const double complex* const restrict overlap_matrix,
                               const size_t num_configs, const size_t max_num_configs) {
	#pragma omp single nowait
	{
	for(size_t k = 0; k < num_configs; k++) {
		for(size_t j = 0; j < num_configs; j++) {
			double complex z;
			if(k <= j)
				z = overlap_matrix[k + max_num_configs*j];
			else
				z = conj(overlap_matrix[j + max_num_configs*k]);
			printf("%g+i%g\t", creal(z), cimag(z));
		}
		printf("\n");
	}
	}
}


/*
prints a 2d vector of complex numbers. could be e.g. hoz or el_coeff
*/
void MCE1_print_vec2(const double complex* const restrict vec,
                     const size_t num_configs, const size_t max_num_configs,
                     const size_t num_sites) {
	#pragma omp single nowait
	{
	for(size_t n = 0; n < num_sites; n++) {
		for(size_t j = 0; j < num_configs; j++) {
			double complex z = vec[n*max_num_configs + j];
			printf("%g+i%g\t", creal(z), cimag(z));
		}
		printf("\n");
	}
	}
}


/*
prints all relevant information in case of a failure of LAPACKE_zpotrf or zpotrs
(e.g. singular overlap matrix or invalid values in either of the vectors)
*/
void print_lapack_solve_error(const size_t num_configs, const size_t max_num_configs,
                              const size_t num_sites,
                              const double complex* const overlap_matrix,
                              const double complex* const vec2_eom_rhs,
                              const double complex* const vec2_hoz,
                              const double complex* const vec2_el_coeff,
                              const int info) {
	#pragma omp single nowait
	{
	printf("solving linear equation failed with info %i\n", info);
	#if PRINT_FULL_DEBUG
		printf("\ninput 5 = overlap_matrix:\n");
		MCE1_print_overlap_matrix(overlap_matrix, num_configs, max_num_configs);
		
		printf("\ninput 7 = vec2_eom_rhs:\n");
		MCE1_print_vec2(vec2_eom_rhs, num_configs, max_num_configs, num_sites);
		
		printf("\nhoz:\n");
		MCE1_print_vec2(vec2_hoz, num_configs, max_num_configs, num_sites);
		printf("\nel_coeff:\n");
		MCE1_print_vec2(vec2_el_coeff, num_configs, max_num_configs, num_sites);
	#endif
	}
}

/*
calculates eigenvectors and eigenvalues of given lattice distortions (vec_x)
in Holstein model (1d or 2d)
used values: vec_x[config_offset + ld_x*n], n= 0..pb->num_sites-1

output:
eigenvectors in pb->BO_eigvecs
eigenvalues starting at (double*) &pb->vec2_buffer[omp_get_thread_num()*pb->num_sites],
written continuously as double array.

returns 0 on success, non-zero value indicates some lapack issue
*/
int MCE_fill_BO_basis(struct param_and_buffer* const pb,
                      const double coupling,
                      const double* const vec_x,
                      const size_t ld_x,
                      size_t config_offset,
                      int force_signs_eigvecs) {
	int info;
	// find separted buffers. vec2_buffer is complex, therefore space for diag AND offdiag
	// and buffersize is max(num_threads, max_num_configs)*num_sites.
	// buffer2 is used for coefficients after projection
	#ifdef _OPENMP
		const int thread_ID = omp_get_thread_num();
	#else
		const int thread_ID = 0;
	#endif
	double* const restrict H_diag = (double*) &pb->vec2_buffer[thread_ID*pb->num_sites];
	double* const restrict H_offdiag = &H_diag[pb->num_sites];
	
	if(pb->num_sites_y == 1) { // 1D and no PBC
		cblas_dcopy(pb->num_sites, // copy real part
		            &vec_x[config_offset], ld_x,
		            H_diag, 1);
		cblas_dscal(pb->num_sites, coupling, H_diag, 1);
		for(size_t i = 0; i < pb->num_sites-1; i++) {
			H_offdiag[i] = -pb->hopping;
		}
		#if PRINT_STORE_OBS_DEBUG
			printf("H_diag   =");
			for(size_t i = 0; i < pb->num_sites; i++) {
				printf("%g ", H_diag[i]);
			}
			printf("\nH_offdiag=");
			for(size_t i = 0; i < pb->num_sites-1; i++) {
				printf("%g ", H_offdiag[i]);
			}
			printf("\n");
		#endif
		// diagonalize
		info
		= LAPACKE_dsteqr(LAPACK_COL_MAJOR, 'I', pb->num_sites,
		                 H_diag, H_offdiag, pb->BO_eigvecs, pb->num_sites);
	} else { // TODO PBC?
		memset(pb->BO_eigvecs, 0, pb->num_sites*pb->num_sites*sizeof(double));
		
		// diagonal part from coupling
		cblas_dcopy(pb->num_sites, &vec_x[config_offset], ld_x,
		            pb->BO_eigvecs, pb->num_sites+1);
		cblas_dscal(pb->num_sites, coupling, pb->BO_eigvecs, pb->num_sites+1); // +1 to follow diagonal instead of row
		
		// offdiagonal part from hopping
		for(size_t y = 0; y < pb->num_sites_y; y++) {
			for(size_t x = 0; x < pb->num_sites_x; x++) {
				const size_t site = y*pb->num_sites_x + x;
				if(x != pb->num_sites_x-1) {
					// from right (left by symmetry)
					pb->BO_eigvecs[site + (site+1)*pb->num_sites] = -pb->hopping;
					//pb->BO_eigvecs[site+1 + site*pb->num_sites] = -pb->hopping;
				}
				// from above (below by symmetry)
				if(y != pb->num_sites_y-1) {
					pb->BO_eigvecs[site + (site+pb->num_sites_x)*pb->num_sites] = -pb->hopping;
					//pb->BO_eigvecs[site+pb->num_sites_x + site*pb->num_sites] = -pb->hopping;
				}
			}
		}
		// diagonalize
		info =
		LAPACKE_dsytrd(LAPACK_COL_MAJOR, CblasUpper, pb->num_sites, pb->BO_eigvecs, pb->num_sites, H_diag, H_offdiag, pb->tau);
		if(!info) {
			info =
			LAPACKE_dorgtr(LAPACK_COL_MAJOR, CblasUpper, pb->num_sites, pb->BO_eigvecs, pb->num_sites, pb->tau);
		}
		if(!info) {
			info =
			LAPACKE_dsteqr(LAPACK_COL_MAJOR, 'V', pb->num_sites, H_diag, H_offdiag, pb->BO_eigvecs, pb->num_sites);
		}
	}
	if(force_signs_eigvecs && !info) {
		// FIXME fix sign of eigenvectors to stop them from flipping
		// TODO get good reference vector set and project to it
		for(size_t n = 0; n < pb->num_sites; n++) {
		if(cblas_ddot(pb->num_sites,
		              &pb->BO_eigvecs[n*pb->num_sites], 1,
		              &((double*)pb->vec2_conf_overlap)[n*pb->num_sites], 1) < 0) {
				cblas_dscal(pb->num_sites, -1, &pb->BO_eigvecs[n*pb->num_sites], 1);
			}
		}
	}
	return info;
}
