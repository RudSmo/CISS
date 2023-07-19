#ifndef SIM_MCEV1_1E_H
#define SIM_MCEV1_1E_H

#include "sim_MCE_defs.h"

#include <stddef.h>
#include <complex.h>

double MCE1_store_observables(const size_t time, const size_t obs_BO_nel_save_dist,
                              struct param_and_buffer* const pb,
                              const double complex* const restrict vec2_hoz,
                              const double complex* const restrict vec2_el_coeff);

void MCE1_joint_derivative_for_rkf(struct param_and_buffer* const pb,
                                   const double complex* restrict vec2_hoz_old,
                                   const double complex* restrict vec2_el_coeff,
                                   double complex* restrict vec2_d_hoz,
                                   double complex* restrict vec2_d_el_coeff);

#if ALLOW_CLONING_MCEv1
void MCE1_clone_configs(struct param_and_buffer* pb,
                        double complex* restrict vec2_hoz,
                        double complex* restrict vec2_el_coeff);
#endif

#endif // SIM_MCEV1_1E_H
