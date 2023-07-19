#ifndef ZPOTRF_OMP_H
	#define ZPOTRF_OMP_H

#include <mkl.h>
#include <complex.h>


/*
  own version of zpotrf using omp parallelization
  PRACE_Advanced_OpenMP_Tutorial.pdf && https://fdesprez.github.io/teaching/par-comput/lectures/slides/L4-OpenMP-pt2-2p.pdf
*/
int zpotrf_omp(const char uplo, const int n,
               double complex* restrict const a, const int lda, int* shared_info);

#endif // ZPOTRF_OMP_H
