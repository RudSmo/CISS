#ifndef SIM_MCE_H
#define SIM_MCE_H

#include <stddef.h>
#include <complex.h>

// to customize the behavior of the MCE program, take a look at sim_MCE_defs.h

int MCE_simulate(size_t num_configs,
                 size_t max_num_configs,
                 int allow_config_delete,
                 int allow_config_clone,
                 size_t num_sites_x,
                 size_t num_sites_y,
                 size_t num_electrons,
                 size_t num_times,
                 double dt,
                 double hopping,
                 double omega,
                 double gamma,
                 double complex* restrict vec2_hoz,  // hoz = harmonic oscillator z
                 double complex* restrict vec2_el_coeff,
                 double complex* restrict vec_coupling_coeff,
                 double complex* restrict vec2_obs_y,
                 double* restrict vec2_obs_nph,
                 double* restrict vec2_obs_nel,
                 double* restrict vec_obs_el_ekin,
                 double* restrict vec2_obs_el_ph_inter,
                 double* restrict vec2_obs_BO_nel,
                 int obs_BO_nel_save_dist,
                 size_t* restrict fullStorage_num_configs,
                 double complex* restrict fullStorage_hoz,
                 double complex* restrict fullStorage_el_coeff,
                 int printProgress);

#endif // SIM_MCE_H
