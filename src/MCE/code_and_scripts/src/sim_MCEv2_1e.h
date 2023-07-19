#ifndef SIM_MCEV2_1E_H
#define SIM_MCEV2_1E_H

#include "sim_MCE_defs.h"

#include <stddef.h>
#include <complex.h>

void MCE2_to_MCE1_wavefunction(struct param_and_buffer* const pb,
                               const double complex* const restrict MCE2_vec2_hoz,
                               const double complex* const restrict MCE2_vec2_el_coeff,
                               const double complex* const restrict MCE2_vec_coupling_coeff,
                               double complex* const restrict MCE1_vec2_hoz,
                               double complex* const restrict MCE1_vec2_el_coeff);

void MCE2_joint_derivative_for_rkf(struct param_and_buffer* const pb,
                                   const double complex* restrict vec2_hoz_old,
                                   const double complex* restrict vec2_el_coeff,
                                   const double complex* restrict vec_coupling_coeff,
                                   double complex* restrict vec2_d_hoz,
                                   double complex* restrict vec2_d_el_coeff,
                                   double complex* restrict vec_d_coupling_coeff);

#if ALLOW_CLONING_MCEv2
void MCE2_reduce_configs(struct param_and_buffer* pb,
                         double complex* restrict vec2_hoz,
                         double complex* restrict vec2_el_coeff,
                         double complex* restrict vec_coupling_coeff);

void MCE2_clone_configs(struct param_and_buffer* pb,
                        double complex* restrict vec2_hoz,
                        double complex* restrict vec2_el_coeff,
                        double complex* restrict vec_coupling_coeff);
#endif

#endif // SIM_MCEV2_1E_H
