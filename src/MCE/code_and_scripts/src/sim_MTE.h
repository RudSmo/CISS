#ifndef EHRENFEST_SIM_H
#define EHRENFEST_SIM_H

#include <stddef.h>
#include <complex.h>

void MTE_simulate(size_t num_samples,
                  size_t num_sites_x,
                  size_t num_sites_y,
                  size_t num_electrons,
                  size_t num_times,
                  double dt_save,
                  double tol_energy_error_per_site,
                  const double* restrict y_initial,
                  const complex double* restrict psi_initial,
                  double hopping,
                  double omega,
                  double gamma,
                  double* restrict y_mean_a_t,	// a_t = all_times
                  double* restrict density_a_t,
                  double* restrict y_squared_a_t,
                  int printProgress,
                  double* restrict y_full_a_t,
                  double* restrict el_kinetic_energy_a_t,
                  double* restrict energy_interaction_a_t);

#endif // EHRENFEST_SIM_H
