#include "sim_MCE.h"
#include "sim_MCE_defs.h"
#include "sim_MCEv1_1e.h"
#include "sim_MCEv2_1e.h"


#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include <signal.h>

#include <mkl.h>
#include <mkl_vsl.h>

#include <math.h>
#include <unistd.h>

#if MEASURE_TIME
	#include <time.h>
#endif

// setting up varaibles from sim_MCE_defs.h

const double complex cplx_zero = 0.;
const double complex cplx_one = 1.;
const double complex cplx_neg_one = -1.;

// Runge Kutta parameters
static const double eps_abs = 1e-10;

static const double error_tolerance_step_repeat = 1.5; // should be slightly larger than 1
static const double time_step_safety = 0.9;            // should be slightly smaller than 1
static const double h_min = 1e-5;
//static const int step_num_check_cloning = 100;

static const double time_diff_check_cloning = .1;	// is devided by MAX(omega, hopping)

int MCE_num_time_steps;

// indicates whether the program shall terminate, triggered by SIGTERM or SIGINT
volatile int sigterm_received = 0;
extern const size_t BO_occupation_quad_integral_num_pointsMAX;

void handle_sig_term_int(int sig) {
	fprintf(stderr, "terminate/interrupt handler invoked\n");
	if(sig == SIGTERM || sig == SIGINT) {
		sigterm_received = 1;
	}
	fflush(stdout);
	fflush(stderr);
}


// TODO consider switching to dormand-prince
// https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method

// all coefficients of RKF45
static const size_t num_rk_steps = 6;
//const double c[7] = {0., 1./4., 3./8., 12./13., 1., 1./2.};
static const double mce_rkf_a[6][5] = {{   0.},
                               {   1./4.},
                               {   3./32.,       9./32.},
                               {1932./2197., -7200./2197.,  7296./2197.},
                               { 439./216. ,    -8.      ,  3680./513. , -845./4104.},
                               {  -8./27.  ,     2.      , -3544./2565., 1859./4104., -11./40.}};

// 5th order is just the error estimate
//const double mce_rkf_b5[6] = {16./135., 0., 6656./12825., 28561./56430., -9./50., 2./55.};
static const double mce_rkf_b4[6] = {25./216., 0., 1408./2565. ,  2197./4104. , -1./5. , 0.};
static const double mce_rkf_b_err[6] = {1./360., 0., -128./4275., -2197./75240., 1./50., 2./55.};

/*
performs one integration step of RKF45 of size h.
the data of the new step are contained in vec2_hoz_buffer and vec2_el_coeff,
the error estimate is written to *error, if the calculation succeded.
returns info return of solving
If it is not zero, don't use any results and abort the simulation
*/
void MCE1_rk_step(struct param_and_buffer* pb, const double h,
                  
                  const double complex* restrict vec2_hoz,
                  const double complex* restrict vec2_el_coeff,
                  
                  double complex* restrict k_hoz[num_rk_steps],
                  double complex* restrict k_el_coeff[num_rk_steps],
                  
                  double complex* restrict vec2_hoz_buffer,
                  double complex* restrict vec2_el_coeff_buffer,
                  
                  double* error) {
	
	for(size_t partial_step = 0; partial_step < num_rk_steps && !*pb->info; partial_step++) {
		// copy old values of hoz and el_coeff
		#pragma omp for
		for(size_t n = 0; n < pb->num_sites; n++) {
			cblas_zcopy(*pb->num_configs,
			            &vec2_hoz[n*pb->max_num_configs], 1,
			            &vec2_hoz_buffer[n*pb->max_num_configs], 1);
			cblas_zcopy(*pb->num_configs,
			            &vec2_el_coeff[n*pb->max_num_configs], 1,
			            &vec2_el_coeff_buffer[n*pb->max_num_configs], 1);
			
			// add previous derivative estimates to state
			for(size_t i = 0; i < partial_step; i++) {
				if(mce_rkf_a[partial_step][i] != 0.) {
					cblas_daxpy(*pb->num_configs*2, h*mce_rkf_a[partial_step][i],
					            (double*) &k_hoz[i][n*pb->max_num_configs], 1,
					            (double*) &vec2_hoz_buffer[n*pb->max_num_configs], 1);
					cblas_daxpy(*pb->num_configs*2, h*mce_rkf_a[partial_step][i],
					            (double*) &k_el_coeff[i][n*pb->max_num_configs], 1,
					            (double*) &vec2_el_coeff_buffer[n*pb->max_num_configs], 1);
				}
			}
		}
		
		// calculate next derivative estimate
		MCE1_joint_derivative_for_rkf(pb,
		                              vec2_hoz_buffer, vec2_el_coeff_buffer,
		                              k_hoz[partial_step], k_el_coeff[partial_step]);
	}
	
	
	if(!*pb->info) {
		double max_thread_error = 0.;
		// collect results
		#pragma omp for schedule(static)
		for(size_t n = 0; n < pb->num_sites; n++) {
			// result (order 4)
			// copy old values of hoz and el_coeff
			cblas_zcopy(*pb->num_configs,
			            &vec2_hoz[n*pb->max_num_configs], 1,
			            &vec2_hoz_buffer[n*pb->max_num_configs], 1);
			cblas_zcopy(*pb->num_configs,
			            &vec2_el_coeff[n*pb->max_num_configs], 1,
			            &vec2_el_coeff_buffer[n*pb->max_num_configs], 1);
			
			
			// add derivative estimates to state
			for(size_t partial_step = 0; partial_step < num_rk_steps; partial_step++) {
				if(mce_rkf_b4[partial_step] != 0.) {
					cblas_daxpy(*pb->num_configs*2, h*mce_rkf_b4[partial_step],
					            (double*) &k_hoz[partial_step][n*pb->max_num_configs], 1,
					            (double*) &vec2_hoz_buffer[n*pb->max_num_configs], 1);
					cblas_daxpy(*pb->num_configs*2, h*mce_rkf_b4[partial_step],
					            (double*) &k_el_coeff[partial_step][n*pb->max_num_configs], 1,
					            (double*) &vec2_el_coeff_buffer[n*pb->max_num_configs], 1);
				}
			}
			
			// error estimate
			for(int i = 0; i < 2; i++) { // FIXME i < 2
				double* vec_diff = (double*) &pb->vec2_buffer[n*pb->max_num_configs];
				
				memset(vec_diff, 0, *pb->num_configs*sizeof(double complex));
				for(size_t partial_step = 0; partial_step < num_rk_steps; partial_step++) {
					if(mce_rkf_b_err[partial_step] != 0.) {
						cblas_daxpy(*pb->num_configs*2, h*mce_rkf_b_err[partial_step],
						            (double*) &(i?k_hoz:k_el_coeff)[partial_step][n*pb->max_num_configs], 1,
						            (double*) vec_diff, 1);
					}
				}
				
				// TODO add relative error tolerance
				
				// search max difference, real and imag part treated indepently
				for(size_t k = 0; k < *pb->num_configs*2; k++) {
					const double diff = fabs(vec_diff[k]);
					if(diff > max_thread_error) {
						max_thread_error = diff;
					}
				}
			}
		}
		
		
		#pragma omp single
			*error = 0.;
		#pragma omp critical(mce1_rkf_step_error_collection)
		{
			if(max_thread_error > *error) {
				*error = max_thread_error;
			}
		}
	}
	#pragma omp barrier
}

static double MCE_rk_error_estimate(const size_t offset, const size_t len, const double h,
                                    const size_t num_steps,
                                    const double* weights,
                                    double* const restrict vec_diff,
                                    const double* const * const vec_steps) {
	// filling buffer with error estimate per site
	memset(vec_diff, 0, len*sizeof(double));
	for(size_t partial_step = 0; partial_step < num_steps; partial_step++) {
		if(weights[partial_step] != 0.) {
			cblas_daxpy(len, h*weights[partial_step],
			            (double*) &vec_steps[partial_step][offset], 1,
			            (double*) vec_diff, 1);
		}
	}
	// TODO add relative error tolerance
	
	double max_diff = 0.;
	// search max difference, real and imag part treated indepently
	for(size_t k = 0; k < len; k++) {
		const double diff = fabs(vec_diff[k]);
		if(diff > max_diff) {
			max_diff = diff;
		}
	}
	return max_diff;
}

void MCE2_rk_step(struct param_and_buffer* pb, const double h,
                  
                  const double complex* restrict vec2_hoz,
                  const double complex* restrict vec2_el_coeff,
                  const double complex* restrict vec_coupling_coeff,
                  
                  double complex* restrict k_hoz[num_rk_steps],
                  double complex* restrict k_el_coeff[num_rk_steps],
                  double complex* restrict k_coupling_coeff[num_rk_steps],
                  
                  double complex* restrict vec2_hoz_buffer,
                  double complex* restrict vec2_el_coeff_buffer,
                  double complex* restrict vec_coupling_coeff_buffer,
                  
                  double* error) {
	
	#ifdef _OPENMP
		// calculated fairest distribution of configurations onto threads
		// *pb->num_configs can change between calls
		const int thread_num = omp_get_thread_num();
		const int num_threads = omp_get_num_threads();
		const ldiv_t ratio = ldiv(*pb->num_configs, num_threads);
		
		const size_t thread_config_start = thread_num*ratio.quot + MIN(thread_num, ratio.rem);
		const size_t thread_config_num = ratio.quot + (thread_num < ratio.rem);
		
		#if PRINT_FULL_DEBUG
			printf("rk_step: thread %i: configs %lu-%lu\n",
			       thread_num, thread_config_start,
			       thread_config_start+thread_config_num);
		#endif
	#else
		const size_t thread_config_start = 0;
		const size_t thread_config_num = *pb->num_configs;
	#endif
	
	for(size_t partial_step = 0; partial_step < num_rk_steps && !*pb->info; partial_step++) {
		// copy old values of
		// hoz
		cblas_zcopy(thread_config_num*pb->num_sites,
		            &vec2_hoz[thread_config_start*pb->num_sites], 1,
		            &vec2_hoz_buffer[thread_config_start*pb->num_sites], 1);
		// el_coeff
		cblas_zcopy(thread_config_num*pb->num_sites,
		            &vec2_el_coeff[thread_config_start*pb->num_sites], 1,
		            &vec2_el_coeff_buffer[thread_config_start*pb->num_sites], 1);
		// coupling_coeff
		cblas_zcopy(thread_config_num,
		            &vec_coupling_coeff[thread_config_start], 1,
		            &vec_coupling_coeff_buffer[thread_config_start], 1);
		
		for(size_t i = 0; i < partial_step; i++) {
			// add previous 
			if(mce_rkf_a[partial_step][i] != 0.) {
				cblas_daxpy(thread_config_num*pb->num_sites*2, h*mce_rkf_a[partial_step][i],
				            (double*) &k_hoz[i][thread_config_start*pb->num_sites], 1,
				            (double*) &vec2_hoz_buffer[thread_config_start*pb->num_sites], 1);
				
				cblas_daxpy(thread_config_num*pb->num_sites*2, h*mce_rkf_a[partial_step][i],
				            (double*) &k_el_coeff[i][thread_config_start*pb->num_sites], 1,
				            (double*) &vec2_el_coeff_buffer[thread_config_start*pb->num_sites], 1);
				
				cblas_daxpy(thread_config_num*2, h*mce_rkf_a[partial_step][i],
				            (double*) &k_coupling_coeff[i][thread_config_start], 1,
				            (double*) &vec_coupling_coeff_buffer[thread_config_start], 1);
			}
		}
		#pragma omp barrier
		#if PRINT_FULL_DEBUG
			#pragma omp single
			{
				printf("vec_coupling_coeff_buffer:\n");
				for(size_t k = 0; k < *pb->num_configs; k++) {
					printf("%g+i%g ",
					       creal(vec_coupling_coeff_buffer[k]),
					       cimag(vec_coupling_coeff_buffer[k]));
				}
				printf("\n");
			}
		#endif
		
		// calculate next derivative estimate
		MCE2_joint_derivative_for_rkf(pb,
		                              vec2_hoz_buffer, vec2_el_coeff_buffer, vec_coupling_coeff_buffer,
		                              k_hoz[partial_step], k_el_coeff[partial_step], k_coupling_coeff[partial_step]);
	}
	
	if(!*pb->info) {
		// collect results
		
		// copy old values of
		// hoz
		cblas_zcopy(thread_config_num*pb->num_sites,
		            &vec2_hoz[thread_config_start*pb->num_sites], 1,
		            &vec2_hoz_buffer[thread_config_start*pb->num_sites], 1);
		// el_coeff
		cblas_zcopy(thread_config_num*pb->num_sites,
		            &vec2_el_coeff[thread_config_start*pb->num_sites], 1,
		            &vec2_el_coeff_buffer[thread_config_start*pb->num_sites], 1);
		// coupling_coeff
		cblas_zcopy(thread_config_num,
		            &vec_coupling_coeff[thread_config_start], 1,
		            &vec_coupling_coeff_buffer[thread_config_start], 1);
		
		for(size_t partial_step = 0; partial_step < num_rk_steps; partial_step++) {
			// add previous 
			if(mce_rkf_b4[partial_step] != 0.) {
				cblas_daxpy(thread_config_num*pb->num_sites*2, h*mce_rkf_b4[partial_step],
				            (double*) &k_hoz[partial_step][thread_config_start*pb->num_sites], 1,
				            (double*) &vec2_hoz_buffer[thread_config_start*pb->num_sites], 1);
				
				cblas_daxpy(thread_config_num*pb->num_sites*2, h*mce_rkf_b4[partial_step],
				            (double*) &k_el_coeff[partial_step][thread_config_start*pb->num_sites], 1,
				            (double*) &vec2_el_coeff_buffer[thread_config_start*pb->num_sites], 1);
				
				cblas_daxpy(thread_config_num*2, h*mce_rkf_b4[partial_step],
				            (double*) &k_coupling_coeff[partial_step][thread_config_start], 1,
				            (double*) &vec_coupling_coeff_buffer[thread_config_start], 1);
			}
		}
		
		// error estimate
		double max_thread_error =
		MCE_rk_error_estimate(thread_config_start*pb->num_sites*2, thread_config_num*pb->num_sites*2, h,
		                      num_rk_steps, mce_rkf_b_err, (double*) &pb->vec2_buffer[thread_config_start*pb->num_sites],
		                      (const double* const *) k_hoz);
		//printf("  phonon error: %g\n", max_thread_error);
		
		double error_temp = 
		MCE_rk_error_estimate(thread_config_start*pb->num_sites*2, thread_config_num*pb->num_sites*2, h,
		                      num_rk_steps, mce_rkf_b_err, (double*) &pb->vec2_buffer[thread_config_start*pb->num_sites],
		                      (const double* const *) k_el_coeff);
		//printf("electron error: %g\n", error_temp);
		if(error_temp > max_thread_error)
			max_thread_error = error_temp;
		
		error_temp = 
		MCE_rk_error_estimate(thread_config_start*2, thread_config_num*2, h,
		                      num_rk_steps, mce_rkf_b_err, (double*) &pb->vec2_buffer[thread_config_start*pb->num_sites],
		                      (const double* const *) k_coupling_coeff);
		//printf("coupling error: %g\n", error_temp);
		if(error_temp > max_thread_error)
			max_thread_error = error_temp;
		
		
		#pragma omp single
			*error = 0.;
		#pragma omp critical(mce2_rkf_step_error_collection)
		{
			if(max_thread_error > *error) {
				*error = max_thread_error;
			}
		}
	}
	#pragma omp barrier
}

/*
integrates the system using Runge-Kutta-Fehlberg4(5) (RKF45)
info return of LAPACKE_zpotrf or LAPACKE_zpotrs is saved in pb_shared->info
If it is not zero, don't use any results and abort the simulation

The value of vec_coupling_coeff determines which simualtion is run:
vec_coupling_coeff == NULL -> MCEv1
vec_coupling_coeff != NULL -> MCEv2
note the different assumptions about memory layout stated in documentation of MCE_simulate
*/
void MCE_rkf_integrate(struct param_and_buffer* pb_shared,
                       const double dt, const size_t num_times, const int obs_BO_nel_save_dist,
                       double complex* restrict vec2_hoz,
                       double complex* restrict vec2_el_coeff,
                       double complex* restrict vec_coupling_coeff,
                       const int printProgress) {
	double complex* k_hoz[num_rk_steps];
	double complex* k_el_coeff[num_rk_steps];
	double complex* k_coupling_coeff[num_rk_steps];
	for(size_t i = 0; i < num_rk_steps; i++) {
		k_hoz[i] = malloc(pb_shared->num_sites*pb_shared->max_num_configs*sizeof(double complex));
		k_el_coeff[i] = malloc(pb_shared->num_sites*pb_shared->max_num_configs*sizeof(double complex));
		k_coupling_coeff[i] = NULL;
		if(vec_coupling_coeff != NULL) {
			k_coupling_coeff[i] = malloc(pb_shared->max_num_configs*sizeof(double complex));
		}
	}
	double complex* vec2_hoz_buffer = malloc(pb_shared->num_sites*pb_shared->max_num_configs*sizeof(double complex));
	double complex* vec2_el_coeff_buffer = malloc(pb_shared->num_sites*pb_shared->max_num_configs*sizeof(double complex));
	
	double complex* vec_coupling_coeff_buffer = NULL;
	if(vec_coupling_coeff != NULL) {
		vec_coupling_coeff_buffer = malloc(pb_shared->max_num_configs*sizeof(double complex));
	}
	
	#if PRINT_FULL_DEBUG
		printf("\nallocated memory addresses:\n");
		printf("vec2_hoz                 : %p\n", vec2_hoz);
		printf("vec2_hoz_buffer          : %p\n\n", vec2_hoz_buffer);
		
		printf("vec2_el_coeff            : %p\n", vec2_el_coeff);
		printf("vec2_el_coeff_buffer     : %p\n\n", vec2_el_coeff_buffer);
		
		printf("vec_coupling_coeff       : %p\n", vec_coupling_coeff);
		printf("vec_coupling_coeff_buffer: %p\n", vec_coupling_coeff_buffer);
		
		printf("\nk_hoz:");
		for(size_t i = 0; i < num_rk_steps; i++) {
			printf(" %p", k_hoz[i]);
		}
		printf("\nk_el_coeff:");
		for(size_t i = 0; i < num_rk_steps; i++) {
			printf(" %p", k_el_coeff[i]);
		}
		printf("\nk_coupling:");
		for(size_t i = 0; i < num_rk_steps; i++) {
			printf(" %p", k_coupling_coeff[i]);
		}
		printf("\n");
	#endif
	
	double error; // error estimate of Runge-Kutta integration
	#pragma omp parallel
	{
		struct param_and_buffer pb = *pb_shared;
		pb.BO_eigvecs = NULL;
		pb.tau = NULL;
		if(*pb.num_configs != pb.max_num_configs || pb.vec2_obs_BO_nel != NULL)
		{
			pb.BO_eigvecs = (double*) malloc(pb.num_sites*pb.num_sites*sizeof(double));
			pb.tau = (double*) malloc(pb.num_sites*sizeof(double));
		}
		
		if (pb.vec2_obs_BO_nel != NULL) {
			#ifdef _OPENMP
				const int thread_ID = omp_get_thread_num();
			#else
				const int thread_ID = 0;
			#endif
			vslNewStream(&pb.rng_stream, VSL_BRNG_MT19937, 12845129u + 76519u*thread_ID);
			pb.vec_buffer_local  = (double*) malloc((BO_OCCUPATION_CALC_METHOD == BO_QUAD_INTEGRAL
			                                         ? MAX(pb.num_sites, BO_occupation_quad_integral_num_pointsMAX)
			                                         : pb.num_sites)
			                                        *sizeof(double));
			pb.vec_buffer_local2 = (double*) malloc((BO_OCCUPATION_CALC_METHOD == BO_QUAD_INTEGRAL
			                                         ? MAX(pb.num_sites*pb.num_sites, BO_occupation_quad_integral_num_pointsMAX)
			                                         : pb.num_sites)
			                                        *sizeof(double));
			pb.vec_buffer_local_sq_occ = (double*) malloc(pb.num_sites*sizeof(double));
		}
		
		
		if(vec_coupling_coeff == NULL) {
			MCE1_store_observables(0, obs_BO_nel_save_dist, &pb, vec2_hoz, vec2_el_coeff);
		} else {
			MCE2_to_MCE1_wavefunction(&pb,
			                          vec2_hoz, vec2_el_coeff, vec_coupling_coeff,
			                          vec2_hoz_buffer, vec2_el_coeff_buffer);
			MCE1_store_observables(0, obs_BO_nel_save_dist, &pb, vec2_hoz_buffer, vec2_el_coeff_buffer);
		}
		
		double complex* vec2_hoz_old = vec2_hoz;
		double complex* vec2_hoz_new = vec2_hoz_buffer;
		
		double complex* vec2_el_coeff_old = vec2_el_coeff;
		double complex* vec2_el_coeff_new = vec2_el_coeff_buffer;
		
		
		double complex* vec_coupling_coeff_old = vec_coupling_coeff;
		double complex* vec_coupling_coeff_new = vec_coupling_coeff_buffer;
		
		double h = 5e-4;
		double real_time = 0.;
		double complex norm2;
		//int check_cloning = step_num_check_cloning;
		double next_clone_check_time = time_diff_check_cloning/MAX(pb.hopping, pb.omega);
		
		for(size_t time = 1; time < num_times && !*pb.info && !sigterm_received; time++) {
			#pragma omp single nowait
				MCE_num_time_steps = 0;
			
			// multiple steps might be required for one dt
			while(real_time < time*dt && !*pb.info && !sigterm_received) {
				double dt_step;
				if(real_time + h < time*dt) {
					dt_step = h;
				} else {
					dt_step = time*dt - real_time;
				}
				#pragma omp single nowait
					MCE_num_time_steps++;
				
				#if PRINT_FULL_DEBUG
					#pragma omp single nowait
						printf("starting rkf step\n");
				#endif
				#if USE_CUDA_PROFILE_MARKS
					nvtxEventAttributes_t eventAttrib = {0};
					eventAttrib.version = NVTX_VERSION;
					eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
					eventAttrib.colorType = NVTX_COLOR_ARGB;
					eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;
					
					eventAttrib.color = COLOR_RED;
					eventAttrib.message.ascii = "rk_step";
					nvtxRangePushEx(&eventAttrib);
				#endif
				if(vec_coupling_coeff == NULL) {
					MCE1_rk_step(&pb, dt_step,
					             vec2_hoz_old, vec2_el_coeff_old,
					             k_hoz, k_el_coeff,
					             vec2_hoz_new, vec2_el_coeff_new,
					             &error);
				} else {
					MCE2_rk_step(&pb, dt_step,
					             vec2_hoz_old, vec2_el_coeff_old, vec_coupling_coeff_old,
					             k_hoz, k_el_coeff, k_coupling_coeff,
					             vec2_hoz_new, vec2_el_coeff_new, vec_coupling_coeff_new,
					             &error);
				}
				
				
				if(!*pb.info) {
					const double error_ratio = error/eps_abs;
					//printf("error_ratio=%g\n", error_ratio);
					if(error_ratio < error_tolerance_step_repeat || h <  2.*h_min) {
						//#pragma omp critical
						//	printf("%i: accept error=%g\n", omp_get_thread_num(), error);
						// accept step
						double complex* swapper = vec2_hoz_new;
						vec2_hoz_new = vec2_hoz_old;
						vec2_hoz_old = swapper;
						
						swapper = vec2_el_coeff_new;
						vec2_el_coeff_new = vec2_el_coeff_old;
						vec2_el_coeff_old = swapper;
						
						swapper = vec_coupling_coeff_new;
						vec_coupling_coeff_new = vec_coupling_coeff_old;
						vec_coupling_coeff_old = swapper;
						
						if(dt_step == h) {
							real_time += h;
						} else {
							real_time = time*dt;
						}
						if(pb.gamma != 0. && *pb.num_configs != pb.max_num_configs && real_time >= next_clone_check_time) {
							
						#if USE_CUDA_PROFILE_MARKS
							nvtxRangePop();
							eventAttrib.color = COLOR_BLUE;
							eventAttrib.message.ascii = "cloning";
							nvtxRangePushEx(&eventAttrib);
						#endif
							const size_t old_num_configs = *pb.num_configs;
							if(vec_coupling_coeff == NULL) {
								#if ALLOW_CLONING_MCEv1
									if(pb.allow_config_clone) {
										MCE1_clone_configs(&pb, vec2_hoz_old, vec2_el_coeff_old);
									}
								#endif
							} else {
								#if ALLOW_CLONING_MCEv2
									if(pb.allow_config_delete) {
										MCE2_reduce_configs(&pb, vec2_hoz_old, vec2_el_coeff_old, vec_coupling_coeff_old);
									}
									if(pb.allow_config_clone) {
										MCE2_clone_configs(&pb, vec2_hoz_old, vec2_el_coeff_old, vec_coupling_coeff_old);
									}
								#endif
							}
							const size_t num_cloned = *pb.num_configs-old_num_configs+*pb.num_configs_to_delete;
							#pragma omp single nowait
							{
								if(*pb.num_configs_to_delete) {
									printf("t=%g, deleted %lu configuration(s), new num_configs=%lu\n", real_time, *pb.num_configs_to_delete, *pb.num_configs-num_cloned);
								}
								if(num_cloned)
								{
									printf("t=%g, cloned %lu configuration(s), new num_configs=%lu\n", real_time, num_cloned, *pb.num_configs);
								}
							}
							next_clone_check_time = real_time + time_diff_check_cloning;
						}
					}
					
					#if USE_CUDA_PROFILE_MARKS
						nvtxRangePop();
					#endif
					
					// adapt time step
					const double h_new = dt_step*fmin(3., fmax(1./3., time_step_safety*pow(1./error_ratio, 1./4.)));
					//printf("%g  \t%g\n", real_time, h_new);
					if(dt_step == h || h_new < dt_step) {
						h = h_new;
					}
					if(h < h_min) {
						h = h_min;
					}
				}
				#pragma omp barrier // synchronize to minimize risk of deadlock on interrupt
			}
			
			if(!*pb.info) {
				if(vec_coupling_coeff == NULL) {
					norm2 = MCE1_store_observables(time, obs_BO_nel_save_dist, &pb, vec2_hoz_old, vec2_el_coeff_old);
				} else {
					MCE2_to_MCE1_wavefunction(&pb,
					                          vec2_hoz_old, vec2_el_coeff_old, vec_coupling_coeff_old,
					                          vec2_hoz_new, vec2_el_coeff_new);
					norm2 = MCE1_store_observables(time, obs_BO_nel_save_dist, &pb, vec2_hoz_new, vec2_el_coeff_new);
				}
				
				#pragma omp single
				{
					if(printProgress) {
						printf("%5.2f/%g\t%g  \t%i  \t%g\n", time*dt, (pb.num_times-1)*dt, creal(norm2), MCE_num_time_steps, h);
						fflush(stdout);
					}
				}
			}
		}
		
		// copy latest values back to the original buffer
		if(vec2_hoz_old != vec2_hoz || vec2_el_coeff_old != vec2_el_coeff) {
			#pragma omp for nowait
			for(size_t n = 0; n < pb.num_sites; n++) {
				cblas_zcopy(*pb.num_configs,
				            &vec2_hoz_old[n*pb.max_num_configs], 1,
				            &vec2_hoz[n*pb.max_num_configs], 1);
				cblas_zcopy(*pb.num_configs,
				            &vec2_el_coeff_old[n*pb.max_num_configs], 1,
				            &vec2_el_coeff[n*pb.max_num_configs], 1);
			}
		}
		
		// free thread local buffers
		free(pb.BO_eigvecs);
		free(pb.tau);
		if (pb.vec2_obs_BO_nel != NULL) {
			vslDeleteStream(&pb.rng_stream);
			free(pb.vec_buffer_local);
			free(pb.vec_buffer_local2);
			free(pb.vec_buffer_local_sq_occ);
		}
	} // end parallel
	
	free(vec2_el_coeff_buffer);
	free(vec2_hoz_buffer);
	for(size_t i = 0; i < num_rk_steps; i++) {
		free(k_hoz[i]);
		free(k_el_coeff[i]);
		free(k_coupling_coeff[i]);
	}
	free(vec_coupling_coeff_buffer);
}

/*
coeffiecents (z=hoz, a=el_coeff) of different configurations for same site are stored continuous.
-> different sites of same configuration are far apart

if any problem occurs in solving the equation of motion (using LAPACKE_zpotrf and LAPACKE_zpotrs),
the simulation is aborted and the errorcode returned

The value of vec_coupling_coeff determines which simualtion is run:
vec_coupling_coeff == NULL -> MCEv1
vec_coupling_coeff != NULL -> MCEv2

The two versions of MCE assume different memory layouts for vec2_el_coeff and vec2_hoz!
both: row major storage
MCEv1: site_y, site_x, config
MCEv2: config, site_y, site_x
*/
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
                 int printProgress) {
	
	if(max_num_configs*sizeof(double complex)%64 != 0) {
		printf("max_num_configs should be a multiple of 4 for best perfomance.\n");
	}
	
	if(obs_BO_nel_save_dist <= 0) {
		fprintf(stderr, "invalid value %i for obs_BO_nel_save_dist. replaced by 1\n", obs_BO_nel_save_dist);
		obs_BO_nel_save_dist = 1;
	}
	
	// set up terminate/interrupt handler
	#if CLUSTER_MODE
		signal(SIGINT, handle_sig_term_int);
		signal(SIGTERM, handle_sig_term_int);
	#else
		struct sigaction new_signal_action;
		new_signal_action.sa_handler = handle_sig_term_int;
		sigemptyset(&new_signal_action.sa_mask);
		new_signal_action.sa_flags = 0;
		sigaction(SIGINT, &new_signal_action, NULL);
		sigaction(SIGTERM, &new_signal_action, NULL);
	#endif
	
	#if MEASURE_TIME
		// start timer for runtime
		struct timespec time_start, time_end;
		clock_gettime(CLOCK_MONOTONIC_RAW, &time_start);
	#endif
	
	#if PRINT_FULL_DEBUG
		printf("simulation inputs:\n");
		printf("num_times = %lu\n", num_times);
		printf("dt        = %g\n", dt);
		printf("hopping   = %g\n", hopping);
		printf("omega     = %g\n", omega);
		printf("gamma     = %g\n\n", gamma);
	#endif
	const size_t num_sites = num_sites_x*num_sites_y;
	//srandom(time(NULL));
	
	
	int num_threads = 1; // TODO MCEv2 can use more than num_sites threads (mostly parallel in num_configs)
	#ifdef _OPENMP
		#if OVERRIDE_OMP_NUM_THREADS
			if(vec2_obs_BO_nel == NULL && vec_coupling_coeff == NULL) {
				num_threads = MIN(num_sites, MAX_NUM_THREADS);
			} else {
				// TODO doen't really allow for good cloning if initial number is smaller than umber of cores
				num_threads = MIN(MAX(num_configs, num_sites), MAX_NUM_THREADS);
			}
			omp_set_num_threads(num_threads);
		#else
			#pragma omp parallel
			{
				#pragma omp single
					num_threads = MIN(num_configs, omp_get_num_threads());
			}
		#endif
	#endif
	
	if(printProgress) {
		printf("num_threads=%i\n", num_threads);
	}
	
	if(num_configs > max_num_configs) {
		fprintf(stderr, "num_configs > max_num_configs! aborting.\n");
		return 1;
	}
	
	if(vec2_obs_y == NULL || vec2_obs_nph == NULL || vec2_obs_nel == NULL) {
		fprintf(stderr, "buffers for all observables must be given!\n");
		return 1;
	}
	
	// buffer weight order in cloning
	int* order = (int*) malloc(sizeof(int)*max_num_configs);
	
	// buffers for basis reduction
	size_t num_configs_to_delete = 0;
	size_t* configs_to_delete = (size_t*) malloc(max_num_configs*sizeof(size_t));
	short* config_checked = (short*) malloc(max_num_configs*sizeof(short));
	size_t* nearest_neighbors = (size_t*) malloc(2*max_num_configs*sizeof(size_t));
	
	// buffers have size num_sites*max(max_num_configs, num_threads), because num_threads can exceed max_num_configs (allowed, but makes no sense...)
	// and cloning requires a larger buffer
	double complex* restrict vec2_buffer = NULL;
	double complex* restrict vec2_buffer2 = NULL;
	
	// buffer for weights from squared norm in coupling term
	double* restrict vec2_density = NULL;
	
	
	int info = 0;
	#if !CLUSTER_MODE
		int err_memalign = 0;
	#endif
	// buffers for electronic equation of motion
	double complex* restrict phononic_overlap_matrix = NULL;
	double complex* restrict total_overlap_matrix = NULL;
	double complex* restrict vec2_conf_overlap = NULL;
	
	#if USE_CUDA
		#if PRINT_FULL_DEBUG
			printf("cuda buffer allocation\n");
		#endif
		cudaError_t cuda_status = cudaSuccess;
		cusolverStatus_t cusolver_status = CUSOLVER_STATUS_SUCCESS;
		
		cusolverDnHandle_t cusolver_handle = NULL;
		
		cuDoubleComplex* restrict R_overlap_gpu = NULL;
		cuDoubleComplex* restrict work_gpu = NULL;
		int* restrict info_gpu = NULL;
		
		
		/*
		cuda_status = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
		
		cuda_status = cusolverDnSetStream(handle, stream);
		*/
		// TODO check status each time
		cuda_status = cusolverDnCreate(&cusolver_handle);
		if(cuda_status != cudaSuccess) {
			printf("cusolverDnCreate failed\n");
			goto abort_simulation;
		}
		
		cuda_status = cudaMallocHost((void**) &phononic_overlap_matrix,
		                             sizeof(double complex)*max_num_configs*max_num_configs);
		if(cuda_status != cudaSuccess) {
			printf("cudaMallocHost phononic_overlap_matrix failed\n");
			goto abort_simulation;
		}
		
		cuda_status = cudaMalloc((void**) &R_overlap_gpu,
		                         sizeof(double complex)*max_num_configs*max_num_configs);
		if(cuda_status != cudaSuccess) {
			printf("cudaMalloc R_overlap_gpu failed\n");
			goto abort_simulation;
		}
		
		cuda_status = cudaMalloc((void**) &info_gpu, sizeof(int)*1);
		if(cuda_status != cudaSuccess) {
			printf("cudaMalloc info_gpu failed\n");
			goto abort_simulation;
		}
		
		
		
		// setting up work buffer for potrf
		int Lwork_gpu;
		cusolver_status = cusolverDnZpotrf_bufferSize(cusolver_handle, CUBLAS_FILL_MODE_UPPER,
		                                             max_num_configs, NULL,
		                                             max_num_configs, &Lwork_gpu);
		cudaDeviceSynchronize();
		if(cusolver_status != CUSOLVER_STATUS_SUCCESS) {
			goto abort_simulation;
		}
		
		cuda_status = cudaMalloc((void**) &work_gpu, sizeof(double complex)*Lwork_gpu);
		if(cuda_status != cudaSuccess) {
			goto abort_simulation;
		}
		
	#else
		#if CLUSTER_MODE
			phononic_overlap_matrix = malloc(max_num_configs*max_num_configs*sizeof(double complex));
		#else
			err_memalign =
			posix_memalign((void **) &phononic_overlap_matrix, 64, max_num_configs*max_num_configs*sizeof(double complex));
		#endif
	#endif // USE_CUDA
	
	
	#if CLUSTER_MODE
		if(vec_coupling_coeff != NULL) {
			total_overlap_matrix = malloc(max_num_configs*max_num_configs*sizeof(double complex));
		}
	#else
		if(vec_coupling_coeff != NULL && !err_memalign) {
			err_memalign =
			posix_memalign((void **) &total_overlap_matrix, 64, max_num_configs*max_num_configs*sizeof(double complex));
		}
	#endif
	
	
	// create buffers
	const size_t max_thread_config_num = MAX(max_num_configs, num_threads);
	const size_t max_sites_config_num = MAX(max_num_configs, num_sites);
	const size_t vec2_conf_overlap_size = max_sites_config_num*max_sites_config_num;
	#if CLUSTER_MODE
		vec2_conf_overlap = malloc(vec2_conf_overlap_size*sizeof(double complex));
		vec2_buffer = malloc(num_sites*max_thread_config_num*sizeof(double complex));
		vec2_buffer2 = malloc(num_sites*max_thread_config_num*sizeof(double complex));
		vec2_density = malloc(num_sites*max_thread_config_num*sizeof(double complex));
		
		
	#else
		if(!err_memalign) {
			err_memalign =
			posix_memalign((void **) &vec2_conf_overlap, 64, vec2_conf_overlap_size*sizeof(double complex));
		}
		if(!err_memalign) {
			err_memalign =
			posix_memalign((void **) &vec2_buffer, 64, num_sites*max_thread_config_num*sizeof(double complex));
		}
		if(!err_memalign) {
			err_memalign =
			posix_memalign((void **) &vec2_buffer2, 64, num_sites*max_thread_config_num*sizeof(double complex));
		}
		if(!err_memalign) {
			err_memalign =
			posix_memalign((void **) &vec2_density, 64, num_sites*max_thread_config_num*sizeof(double complex));
		}
		
		if(err_memalign) {
			fprintf(stderr, "memaligend alloc failed\n");
			info = 1;
			goto abort_simulation;
		}
	#endif
	
	
	// initial observables
	memset(vec2_obs_y,   0, num_sites*num_times*sizeof(double complex));
	memset(vec2_obs_nph, 0, num_sites*num_times*sizeof(double));
	memset(vec2_obs_nel, 0, num_sites*num_times*sizeof(double));
	memset(vec_obs_el_ekin, 0, num_times*sizeof(double));
	memset(vec2_obs_el_ph_inter, 0, num_sites*num_times*sizeof(double));
	if(vec2_obs_BO_nel != NULL) {
		memset(vec2_obs_BO_nel, 0, (1+(num_times-1)/obs_BO_nel_save_dist)*num_sites*sizeof(double));
	}
	
	
	if(fullStorage_num_configs) {
		memset(fullStorage_num_configs, 0, num_times*sizeof(size_t));
	}
	if(fullStorage_el_coeff) {
		memset(fullStorage_el_coeff, 0, num_sites*max_num_configs*num_times*sizeof(double complex));
	}
	if(fullStorage_hoz) {
		memset(fullStorage_hoz, 0, num_sites*max_num_configs*num_times*sizeof(double complex));
	}
	
	
	#if USE_CUDA
		double runtime_cuda = 0.;
		double runtime_cpu = 0.;
	#endif
	
	struct param_and_buffer pb = {// Hamiltonian parameters
	                              hopping, omega, gamma,
	                              // system sizes
	                              &num_configs, max_num_configs,
	                              num_sites, num_sites_x, num_sites_y,
	                              num_electrons, num_times,
	                              // observables
	                              vec2_obs_y, vec2_obs_nph, vec2_obs_nel,
	                              vec_obs_el_ekin, vec2_obs_el_ph_inter,
	                              vec2_obs_BO_nel,
	                              fullStorage_num_configs, fullStorage_el_coeff, fullStorage_hoz,
	                              // buffers
	                              phononic_overlap_matrix,
	                              total_overlap_matrix,
	                              vec2_conf_overlap,
	                              vec2_buffer, vec2_buffer2, vec2_density,
	                              NULL, NULL, NULL, NULL, // rng and small thread local buffers
	                              #if USE_CUDA
	                              cusolver_handle,
	                              R_overlap_gpu,
	                              work_gpu, Lwork_gpu,
	                              info_gpu, cudaSuccess, CUSOLVER_STATUS_SUCCESS,
	                              #endif
	                              // simulation related
	                              &info, printProgress, order,
	                              allow_config_delete, allow_config_clone,
	                              NULL, NULL, // thread local buffers BO_eigvecs, tau
	                              &num_configs_to_delete, configs_to_delete,
	                              config_checked, nearest_neighbors
	                              #if USE_CUDA
	                              , &runtime_cuda, &runtime_cpu
	                              #endif
	                              };
	
	#if PRINT_FULL_DEBUG
		printf("starting simulation\n");
	#endif
	
	MCE_rkf_integrate(&pb, dt, num_times, obs_BO_nel_save_dist, vec2_hoz, vec2_el_coeff, vec_coupling_coeff, printProgress);
	
	#if PRINT_FULL_DEBUG
		printf("finished simulation\n");
	#endif
	
	#if USE_CUDA
		cuda_status = pb.cuda_status;
		cusolver_status = pb.cusolver_status;
	#endif
	
#if USE_CUDA || !CLUSTER_MODE
abort_simulation:
#endif
	
	#if USE_CUDA
		if(cuda_status != cudaSuccess) {
			printf("CUDA error %i: %s\n", cuda_status, cudaGetErrorString(cuda_status));
			info = CUDA_INFO_RETURN;
		}
		if(cusolver_status != CUSOLVER_STATUS_SUCCESS) {
			printf("cusolver error %i\n", cusolver_status);
			info = CUDA_INFO_RETURN;
		}
		cusolverDnDestroy(cusolver_handle);
		
		cudaFree(phononic_overlap_matrix);
		cudaFree(R_overlap_gpu);
		cudaFree(work_gpu);
		cudaFree(info_gpu);
	#else
		free(phononic_overlap_matrix);
	#endif
	free(total_overlap_matrix);
	free(vec2_conf_overlap);
	free(order);
	free(vec2_density);
	free(vec2_buffer2);
	free(vec2_buffer);
	free(nearest_neighbors);
	free(config_checked);
	free(configs_to_delete);
	
	#if MEASURE_TIME
		clock_gettime(CLOCK_MONOTONIC_RAW, &time_end);
		const double runtime = (double) (time_end.tv_sec - time_start.tv_sec) + (time_end.tv_nsec - time_start.tv_nsec) / 1e9;
		if(printProgress) {
			printf("runtime simulation: %g\n", runtime);
			#if USE_CUDA
				printf("runtime cuda_zpotrf: %g\nruntime  cpu_zpotrf: %g\n", runtime_cuda, runtime_cpu);
			#endif
		}
	#endif
	return info;
}
