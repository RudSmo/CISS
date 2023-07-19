#ifndef MCE_UTIL_H
#define MCE_UTIL_H

#include <stddef.h>
#include <complex.h>
#include "sim_MCE_defs.h"

/*
prints the overlap matrix. the hermitean property is assumed and only the
upper triangle (CblasColMajor) storage referenced
*/
void MCE1_print_overlap_matrix(const double complex* const restrict overlap_matrix,
                               const size_t num_configs, const size_t max_num_configs);


/*
prints a 2d vector of complex numbers. could be e.g. hoz or el_coeff
*/
void MCE1_print_vec2(const double complex* const restrict vec,
                     const size_t num_configs, const size_t max_num_configs,
                     const size_t num_sites);

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
                              const int info);

/*
if force_signs_eigvecs, pb->vec2_conf_overlap, the signs are fixed by ensuring
a positive scalar product with a basis prestored in pb->vec2_conf_overlap
pb->vec2_buffer is overwritten
coupling given seperat to allow for calling it with scaled real part of hoz, e.g. <x_l>
the desired eigenvectors are stored column wise LAPACK_COL_MAJOR in pb->BO_eigvecs
returns info form lapack
*/
int MCE_fill_BO_basis(struct param_and_buffer* const pb,
                      const double coupling,
                      const double* const vec_x,
                      const size_t ld_x,
                      size_t config_offset,
                      int force_signs_eigvecs);

#endif // MCE_UTIL_H
