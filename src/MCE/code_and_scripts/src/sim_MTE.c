//gcc -I"${MKLROOT}/include" -m64 -Wall -O3 -fopenmp ehrenfest_sim.c -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
//gcc -I"${MKLROOT}/include" -m64 -Wall -g -fopenmp ehrenfest_sim.c -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
#include "sim_MTE.h"
#include <string.h>
#include <stdio.h>

//#define MKL_INT int
//#define MKL_COMPLEX8 complex
//#define MKL_COMPLEX16 complex double
#include <mkl.h>


#include <time.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>

#include <omp.h>

#define MAX_NUM_THREADS 4
#define GSL_RANGE_CHECK_OFF
#define HAVE_INLINE
#define SQRT2 1.4142135623730950488

static const double complex cplx_neg_one = -1.;

// TODO 2d Hamiltonian
// stores the hamiltonian as symmetric tridiagonal matrix, diagonal first, followed by offdiagonal
static inline void MTE_fill_hamiltonian(double* const restrict Hel_diag,
                                        double* const restrict Hel_offdiag,
                                        const size_t num_sites_x,
                                        const double* const restrict x,
                                        const double coupling,
                                        const double hopping) {
	cblas_dcopy(num_sites_x, x, 2, Hel_diag, 1);
	cblas_dscal(num_sites_x, -coupling, Hel_diag, 1);
	for(size_t i = 0; i < num_sites_x-1; i++) {
		Hel_offdiag[i] = -hopping;
	}
}

// LAPACK_COL_MAJOR and just upper part
static inline void MTE_fill_hamiltonian2d(double* const restrict Hel,
                                          const size_t num_sites_x,
                                          const size_t num_sites_y,
                                          const double* const restrict x,
                                          const double coupling,
                                          const double hopping) {
	const size_t num_sites = num_sites_x*num_sites_y;
	memset(Hel, 0, num_sites*num_sites*sizeof(double));
	
	// diagonal part from coupling
	cblas_dcopy(num_sites, x, 2, Hel, num_sites+1);
	cblas_dscal(num_sites, -coupling, Hel, num_sites+1);
	
	// offdiagonal part from hopping
	for(size_t y = 0; y < num_sites_y; y++) {
		for(size_t x = 0; x < num_sites_x; x++) {
			const size_t site = y*num_sites_x + x;
			if(x != num_sites_x-1) {
				// from right (left by symmetry)
				Hel[site + (site+1)*num_sites] = -hopping;
				//Hel[site+1 + site*num_sites] = -hopping;
			}
			// from above (below by symmetry)
			if(y != num_sites_y-1) {
				Hel[site + (site+num_sites_x)*num_sites] = -hopping;
				//Hel[site+num_sites_x + site*num_sites] = -hopping;
			}
		}
	}
	/*
	#pragma omp single
	{
		for(size_t i = 0; i < num_sites; i++) {
			for(size_t j = 0; j < num_sites; j++) {
				printf("%+6f\t", Hel[i+j*num_sites]);
			}
			printf("\n");
		}
		printf("\n\n");
		exit(1);
	}
	*/
}


// adds the electron kinetic energy divided by hopping coefficent to e_kin_sum (thread local)
static inline void MTE_add_el_kinetic_energy(const size_t num_sites_x, const size_t num_sites_y, const size_t num_electrons,
                                             const complex double* const restrict psi,
                                             double* const restrict e_kin_sum,
                                             complex double* const restrict psi_buffer) {
	const size_t num_sites = num_sites_x*num_sites_y;
	if(num_sites_y == 1) {
		// 1d
		for(size_t e = 0; e < num_electrons; e++) {
			memset(&psi_buffer[e*num_sites], 0, num_sites*sizeof(complex double));
			
			// super-diagonal 
			cblas_zaxpy(num_sites_x-1, &cplx_neg_one, &psi[e*num_sites+1], 1, &psi_buffer[e*num_sites+0], 1);
			
			// sub-diagonal
			cblas_zaxpy(num_sites_x-1, &cplx_neg_one, &psi[e*num_sites+0], 1, &psi_buffer[e*num_sites+1], 1);
			
			// conj(psi).T @ psi_buffer
			complex double e_kin;
			cblas_zdotc_sub(num_sites, &psi[e*num_sites], 1, &psi_buffer[e*num_sites], 1, &e_kin);
			e_kin_sum[e] += creal(e_kin);
		}
	} else {
		// 2d
		for(size_t e = 0; e < num_electrons; e++) {
			memset(&psi_buffer[e*num_sites], 0, num_sites*sizeof(complex double));
			for(size_t y = 0; y < num_sites_y; y++) {
				const size_t offset = e*(num_sites_y + y)*num_sites_x;
				
				// hop from right
				cblas_zaxpy(num_sites_x-1, &cplx_neg_one, &psi[offset+1], 1, &psi_buffer[offset+0], 1);
				// hop from left
				cblas_zaxpy(num_sites_x-1, &cplx_neg_one, &psi[offset+0], 1, &psi_buffer[offset+1], 1);
				
				if(y != num_sites_y-1) {
					// hop from above
					cblas_zaxpy(num_sites_x, &cplx_neg_one, &psi[offset+num_sites_x], 1, &psi_buffer[offset+0], 1);
					// hop from below
					cblas_zaxpy(num_sites_x, &cplx_neg_one, &psi[offset+0], 1, &psi_buffer[offset+num_sites_x], 1);
				}
				
			}
			// conj(psi).T @ psi_buffer
			complex double e_kin;
			cblas_zdotc_sub(num_sites, &psi[e*num_sites], 1, &psi_buffer[e*num_sites], 1, &e_kin);
			e_kin_sum[e] += creal(e_kin);
		}
	}
}


static inline void MTE_psi_basis_change(complex double* const restrict psi_new,
                                        const complex double* const restrict psi_old,
                                        const double* const restrict eigvecs, const CBLAS_TRANSPOSE transpose,
                                        const size_t num_sites, const size_t num_electrons) {
	for(size_t e = 0; e < num_electrons; e++) {
		// real part
		cblas_dgemv(CblasColMajor, transpose, num_sites, num_sites,
		            1., eigvecs, num_sites, ((const double*) &psi_old[e*num_sites])+0, 2,
		            0., ((double*) &psi_new[e*num_sites])+0, 2);
		// imag part
		cblas_dgemv(CblasColMajor, transpose, num_sites, num_sites,
		            1., eigvecs, num_sites, ((const double*) &psi_old[e*num_sites])+1, 2,
		            0., ((double*) &psi_new[e*num_sites])+1, 2);
	}
}


static inline void MTE_update_electronic_density(const complex double* const restrict psi,
                                                 double* const restrict density,
                                                 const size_t num_sites, const size_t num_electrons) {
	for(size_t s=0; s < num_sites; s++) {
		double abs2 = 0.;
		for(size_t e=0; e < num_electrons; e++) {
			complex double p = psi[e*num_sites+s];
			abs2 += creal(p)*creal(p) + cimag(p)*cimag(p);
		}
		density[s] = abs2;
	}
}


static inline void MTE_scale_psi_eig(const double* const restrict Hel_diag, 
                                     complex double* restrict psi_eig,
                                     const double dt,
                                     const size_t num_sites, const size_t num_electrons) {
                                 
	for(size_t i=0; i < num_sites; i++) {
		const complex double phase_fact = cexp(-1.i*dt*Hel_diag[i]);
		for(size_t e=0; e < num_electrons; e++) {
			psi_eig[e*num_sites+i] = phase_fact*psi_eig[e*num_sites+i];
		}
	}
}


void MTE_store_observables(const size_t time,
                           const size_t num_samples,
                           const size_t num_sites_x,
                           const size_t num_sites_y,
                           const size_t num_electrons,
                           
                           double* restrict y_mean_a_t,
                           const double* const restrict curr_y_full,
                           double* const restrict curr_y_mean_sum,
                           
                           double* restrict y_squared_a_t,
                           double* const restrict curr_y_sq_sum,
                           
                           double* restrict density_a_t,
                           const double* const restrict curr_density,
                           double* const restrict curr_density_sum,
                           
                           double* const restrict el_kinetic_energy_a_t,
                           double* const restrict e_kin_sum,
                           const complex double* const restrict vec_psi_full,
                           complex double* const restrict psi_buffer,
                           
                           double* restrict energy_interaction_a_t,
                           double* restrict const interaction_sum) {
	const size_t num_sites = num_sites_x*num_sites_y;
	// shift pointers and reset buffers
	y_mean_a_t  += time*2*num_sites;
	density_a_t += time*num_sites;
	
	memset(curr_y_mean_sum, 0, sizeof(double)*2*num_sites);
	memset(curr_density_sum, 0, sizeof(double)*num_sites);
	if(y_squared_a_t) {
		y_squared_a_t += time*2*num_sites;
		memset(curr_y_sq_sum, 0, sizeof(double)*2*num_sites);
	}
	if(el_kinetic_energy_a_t) {
		memset(e_kin_sum, 0, num_electrons*sizeof(double));
	}
	if(energy_interaction_a_t) {
		memset(interaction_sum, 0, num_sites*sizeof(double));
	}
	#pragma omp for schedule(static) nowait
	for(size_t n = 0; n < num_samples; n++) {
		const double* y_full = curr_y_full + n*2*num_sites;
		for(size_t i = 0 ; i < 2*num_sites; i++) {
			curr_y_mean_sum[i] += y_full[i];
		}
		
		if(y_squared_a_t) {
			for(size_t i = 0 ; i < 2*num_sites; i++) {
				curr_y_sq_sum[i] += y_full[i]*y_full[i];
			}
		}
		
		const double* d = curr_density + n*num_sites;
		for(size_t i = 0 ; i < num_sites; i++) {
			curr_density_sum[i] += d[i];
			if(energy_interaction_a_t) {
				interaction_sum[i] += d[i] * y_full[2*i];
			}
		}
		//cblas_daxpy(num_sites, 1., curr_density + n*num_sites, 1, curr_density_sum, 1);
		
		if(el_kinetic_energy_a_t) {
			MTE_add_el_kinetic_energy(num_sites_x, num_sites_y, num_electrons,
			                          vec_psi_full + n*num_electrons*num_sites,
			                          e_kin_sum, psi_buffer);
		}
	}
	
	
	#pragma omp critical(observable_mean_addition)
	{
		cblas_daxpy(2*num_sites, 1./num_samples, curr_y_mean_sum, 1, y_mean_a_t, 1);
		if(y_squared_a_t) {
			cblas_daxpy(2*num_sites, 1./num_samples, curr_y_sq_sum, 1, y_squared_a_t, 1);
		}
		cblas_daxpy(num_sites, 1./num_samples, curr_density_sum, 1, density_a_t, 1);
		
		if(el_kinetic_energy_a_t) {
			cblas_daxpy(num_electrons, 1./num_samples, e_kin_sum, 1, el_kinetic_energy_a_t + time*num_electrons, 1);
		}
		if(energy_interaction_a_t) {
			cblas_daxpy(num_sites, -SQRT2/num_samples, interaction_sum, 1, energy_interaction_a_t + time*num_sites, 1);
		}
	}
}



/** calculates the multitrajectrory ehrenfest propagation of the holstein model
 * the arrays y_mean_a_t and density_a_t are used to return the results, but must be initialized to zero before calling this function
 */
void MTE_simulate(const size_t num_samples,
                  const size_t num_sites_x,
                  const size_t num_sites_y,
                  const size_t num_electrons,
                  const size_t num_times,
                  const double dt_save,
                  const double tol_energy_error_per_site,
                  const double* const restrict y_initial,
                  const complex double* const restrict psi_initial,
                  const double hopping,
                  const double omega,
                  const double gamma,
                  double* const restrict y_mean_a_t,	// a_t = all_times
                  double* const restrict density_a_t,
                  double* const restrict y_squared_a_t,
                  const int printProgress,
                  double* const restrict y_full_a_t,
                  double* const restrict el_kinetic_energy_a_t,
                  double* const restrict energy_interaction_a_t) {
	// regular stepsize, if larger then dt_save, never used
	double dt_step_reg = sqrt(10*tol_energy_error_per_site*hopping/omega/gamma/gamma);
	
	if(printProgress) {
		printf("num_sites_x=%lu, num_sites_y = %lu\n", num_sites_x, num_sites_y);
		printf("dt_save = %g, dt_step_reg = %g\n", dt_save, dt_step_reg);
	}
	
	const size_t num_sites = num_sites_x*num_sites_y;
	
	memset(y_mean_a_t,  0, num_times*num_sites*2*sizeof(double));
	memset(density_a_t, 0, num_times*num_sites*sizeof(double));
	if(el_kinetic_energy_a_t) {
		memset(el_kinetic_energy_a_t, 0, num_times*num_electrons*sizeof(double));
	}
	if(energy_interaction_a_t) {
		memset(energy_interaction_a_t, 0, num_times*num_sites*sizeof(double));
	}
	
	const double coupling = gamma*sqrt(2.);
	
	// main buffer for current state
	complex double* vec_psi_full     = (complex double*) malloc(num_samples*num_electrons*num_sites*sizeof(complex double));
	complex double* vec_psi_eig_full = (complex double*) malloc(num_samples*num_electrons*num_sites*sizeof(complex double));
	double* vec_y_full               = (double*) malloc(num_samples*num_sites*2*sizeof(double));
	double* vec_density              = (double*) malloc(num_samples*num_sites*sizeof(double));
	
	
	// initialize for all samples with initial state
	memcpy(vec_y_full, y_initial, num_samples*num_sites*2*sizeof(double));
	memcpy(vec_psi_full, psi_initial, num_samples*num_electrons*num_sites*sizeof(complex double));
	
	// (maybe) store initial full configuration
	if(y_full_a_t) {
		memcpy(y_full_a_t, vec_y_full, num_samples*num_sites*2*sizeof(double));
	}
	
	// limit the number of threads to number of samples.
	// otherwise storage of observables maybe not safe due to race condition
	// (and waste of ressources anyway)
	#ifdef _OPENMP
	omp_set_num_threads(num_samples < MAX_NUM_THREADS? num_samples : MAX_NUM_THREADS);
	#endif
	
	#pragma omp parallel
	{
	// buffers for caluclations of observables
	double* vec_y_sum          = (double*) malloc(num_sites*2*sizeof(double));
	double* vec_density_sum    = (double*) malloc(num_sites*sizeof(double));
	double* vec_y_square       = y_squared_a_t ? (double*) malloc(num_sites*2*sizeof(double)) : NULL;
	complex double* psi_buffer = el_kinetic_energy_a_t?(complex double*) malloc(num_electrons*num_sites*sizeof(complex double)):NULL;
	double* e_kin_sum          = el_kinetic_energy_a_t?(double*) malloc(num_electrons*sizeof(double)):NULL;
	double* interaction_sum    = energy_interaction_a_t?(double*) malloc(num_sites*sizeof(double)):NULL;
	
	// buffer for tridiagonal hamiltonian
	double* Hel_diag = (double*) malloc(num_sites*sizeof(double));
	double* Hel_offdiag = (double*) malloc((num_sites-1)*sizeof(double));
	double* eigvecs = (double*) malloc(num_sites*num_sites*sizeof(double));
	
	// additional worksapce for diagonalization of Hamiltoina in 2D or with PBC
	double* tau      = num_sites_y > 1 ? (double*) malloc(num_sites*sizeof(double)) : NULL;
	
	// parameters for lapacks dstemr
	char char_I = 'I';
	char char_V = 'V';
	int Hel_size = num_sites;
	double* work = (double*) malloc(sizeof(double)*(2*num_sites-2));
	int info;
	
	
	// ########################### start of simulation ###########################
	
	// update the full density, required for storage and coupling strength
	#pragma omp for schedule(static) nowait
	for(size_t n=0; n < num_samples; n++) {
		complex double* view_psi = vec_psi_full + n*num_electrons*num_sites;
		double* view_density = vec_density + n*num_sites;
		
		MTE_update_electronic_density(view_psi, view_density, num_sites, num_electrons);
	}
	
	MTE_store_observables(0, num_samples, num_sites_x, num_sites_y, num_electrons,
	                      y_mean_a_t, vec_y_full, vec_y_sum,
	                      y_squared_a_t, vec_y_square,
	                      density_a_t, vec_density, vec_density_sum,
	                      el_kinetic_energy_a_t, e_kin_sum, vec_psi_full, psi_buffer,
	                      energy_interaction_a_t, interaction_sum);
	double time = 0.;
	const double time_deviation_tol = 1e-10;
	for(size_t t = 1; t < num_times; t++) {
		double time_target = t*dt_save;
		while(time < time_target - time_deviation_tol) {
			double dt_step;
			if(time + dt_step_reg <= time_target - time_deviation_tol) {
				dt_step = dt_step_reg;
			} else {
				dt_step = time_target - time;
			}
			double cos_omega_half_dt_step = cos(omega*.5*dt_step);
			double sin_omega_half_dt_step = sin(omega*.5*dt_step);
			
			// Hamiltonian is not time dependent, otherwise belongs to end of loop
			time += dt_step;
			
			#pragma omp for schedule(static) nowait
			for(size_t n = 0; n < num_samples; n++) {
				// those two vectors are interweaved, always use stride 2 to access
				double* const x = vec_y_full + n*2*num_sites+0;
				double* const p = vec_y_full + n*2*num_sites+1;
				
				complex double* view_psi = vec_psi_full + n*num_electrons*num_sites;
				complex double* view_psi_eig = vec_psi_eig_full + n*num_electrons*num_sites;
				double* view_density = vec_density + n*num_sites;
				
				
				// half step harmonic oscillator
				cblas_drot(num_sites, x, 2, p, 2, cos_omega_half_dt_step, sin_omega_half_dt_step);
				
				// half step coupling
				cblas_daxpy(num_sites, .5*dt_step*coupling, view_density, 1, p, 2);
				
				
				// full step quantum electrons via ED
				if(num_sites_y == 1) {
					// 1d, no PBC, just open boundaries
					MTE_fill_hamiltonian(Hel_diag, Hel_offdiag, num_sites_x, x, coupling, hopping);
					// eigvals in Hel_diag, eigenvectors in eigvecs
					dsteqr_(&char_I, &Hel_size, Hel_diag,  Hel_offdiag,
						    eigvecs, &Hel_size, work, &info);
					/*
					info =
					lapacke_dsteqr_work(LAPACK_COL_MAJOR, 'I', num_sites, Hel_diag,  Hel_offdiag,
						                eigvecs, num_sites, work);
					*/
				} else {
					// 2d (could also do PBC, also in 1D)
					MTE_fill_hamiltonian2d(eigvecs, num_sites_x, num_sites_y, x, coupling, hopping);
					LAPACKE_dsytrd(LAPACK_COL_MAJOR, 'U', num_sites, eigvecs, num_sites, Hel_diag, Hel_offdiag, tau);
					LAPACKE_dorgtr(LAPACK_COL_MAJOR, 'U', num_sites, eigvecs, num_sites, tau);
					dsteqr_(&char_V, &Hel_size, Hel_diag,  Hel_offdiag,
						    eigvecs, &Hel_size, work, &info);
				}
				
				// psi_eig = eigvec @ (exp(eigval...) * (eigvecs.T @ psi_initial))
				// function does it for all electrons stored as (num_electrons, num_sites), 
				MTE_psi_basis_change(view_psi_eig, view_psi, eigvecs, CblasTrans, num_sites, num_electrons);
				MTE_scale_psi_eig(Hel_diag, view_psi_eig, dt_step, num_sites, num_electrons);
				MTE_psi_basis_change(view_psi, view_psi_eig, eigvecs, CblasNoTrans, num_sites, num_electrons);
				
				MTE_update_electronic_density(view_psi, view_density, num_sites, num_electrons);
				
				
				// second half step coupling
				cblas_daxpy(num_sites, .5*dt_step*coupling, view_density, 1, p, 2);
				
				// second half step harmonic oscillator
				cblas_drot(num_sites, x, 2, p, 2, cos_omega_half_dt_step, sin_omega_half_dt_step);
			}
		}
		
		MTE_store_observables(t, num_samples, num_sites_x, num_sites_y, num_electrons,
		                      y_mean_a_t, vec_y_full, vec_y_sum,
		                      y_squared_a_t, vec_y_square,
		                      density_a_t, vec_density, vec_density_sum,
		                      el_kinetic_energy_a_t, e_kin_sum, vec_psi_full, psi_buffer,
		                      energy_interaction_a_t, interaction_sum);
		
		#pragma omp barrier
		#pragma omp single
		{
			if(y_full_a_t) {
				memcpy(&y_full_a_t[t*num_samples*num_sites_x*2], vec_y_full, sizeof(double)*num_samples*num_sites_x*2); 
			}
			
		}
		
		#pragma omp single nowait
		if(printProgress && t%50 == 0) {
			printf("%g\n", (100.*t)/num_times);
		}
		
	}
	
	free(work);
	free(tau);
	free(eigvecs);
	free(Hel_offdiag);
	free(Hel_diag);
	free(vec_y_sum);
	free(vec_density_sum);
	free(vec_y_square);
	free(psi_buffer);
	free(e_kin_sum);
	free(interaction_sum);
	} // end of omp parallel
	
	free(vec_psi_eig_full);
	free(vec_psi_full);
	free(vec_density);
	free(vec_y_full);
}

/*
#ifndef PY_VERSION_HEX // compile the main only, if not using cython
int main() {
	
	size_t num_samples = 4000;
	size_t num_sites = 2;
	size_t num_electrons = 1;
	size_t num_times = 200;
	double dt = 0.1;
	double* y_initial = (double*) calloc(num_samples*num_sites*2, sizeof(double));
	y_initial[0] = 1;
	y_initial[3] = 1;
	
	double y[] = {-8.93836900e-01, -1.32308935e+00, -6.82151505e-03, -1.77738098e+00,
				  -3.04824688e-01, -4.84881687e-01, -9.92888159e-02, -4.33259869e-01,
				   2.72740343e-01,  5.91564665e-01,  6.99431922e-02, -8.45020536e-01,
				   9.31910459e-01,  1.20417783e+00, -3.80717809e-01, 2.61532176e-01,
				   2.14787479e-01, -1.31256759e+00,  5.51228333e-01, 1.47151769e-03,
				  -4.13331443e-01,  1.74214270e-01, -5.03656106e-01,  7.75652632e-01,
				   1.45334111e+00, -1.73846001e-01};
	for(size_t i=0; i<2*num_sites; i++) {
		y_initial[i] = y[i];
	}
	
	complex double* psi_initial = (complex double*) calloc(num_samples*num_electrons*num_sites, sizeof(complex double));
	printf("placing electrons at:\n");
	for(size_t i=0; i<num_electrons; i++) {
		size_t electron_site = 2*i;
		printf("%lu ", electron_site);
		for(size_t n=0; n<num_samples; n++) {
			psi_initial[(n*num_electrons+i)*num_sites + electron_site] = 1.+0i;
		}
	}
	printf("\n\n");
	double hopping = 1.;
	double omega = 1.;
	double gamma = .0;
	double* y_all_times = (double*) malloc(num_times*num_sites*2* sizeof(double));
	double* density_all_times  = (double*) malloc(num_times*num_sites* sizeof(double));
	double* e_kin = (double*) malloc(sizeof(double)*num_times*num_electrons);
	double* e_inter = (double*) malloc(sizeof(double)*num_times*num_sites);
	
	simulate_ehrenfest_holstein(num_samples, num_sites, num_electrons, num_times, dt,
		                        y_initial, psi_initial,
		                        hopping, omega, gamma,
		                        y_all_times, density_all_times, NULL, 1,
		                        0, NULL, 0, NULL, NULL,		// histogram parameters
		                        NULL, e_kin, e_inter); // full storage parameter
	
	printf("ekin\teinter\n");
	for(size_t t=0; t<num_times; t++) {
		for(size_t e=0; e<num_electrons; e++) {
			printf("%g\t", e_kin[t*num_electrons + e]);
		}
		for(size_t n=0; n<num_sites; n++) {
			printf("%g\t", e_inter[t*num_sites + n]);
		}
		printf("\n");
	}
	
	
	free(e_inter);
	free(e_kin);
	free(density_all_times);
	free(y_all_times);
	free(psi_initial);
	free(y_initial);
	
	return 0;
}
#endif
*/
