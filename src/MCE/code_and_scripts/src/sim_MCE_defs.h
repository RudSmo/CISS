#ifndef SIM_MCE_DEFS_H
#define SIM_MCE_DEFS_H
#include <stddef.h>

//#define CLUSTER_MODE 1
#ifndef CLUSTER_MODE
	#warning CLUSTER_MODE should be given as compile flag, falling back to 1
	#define CLUSTER_MODE 1
#endif
#define ALLOW_CLONING_MCEv1 0 // TODO cleanup
#define ALLOW_CLONING_MCEv2 (1 && !CLUSTER_MODE) // sim_MCEv2_1e.c:739 vdSubI not known on cluster

#define MCEv2_REFINE_SOLUTION 0 // 0 recommended since no different in result was observed but runtime is ~50% larger

#define BO_MID_EACH 0
#define BO_MID_PAIR 1
#define BO_MC_INTEGRAL 2
#define BO_QUAD_INTEGRAL 3

#define BO_OCCUPATION_CALC_METHOD BO_QUAD_INTEGRAL

// global variable for Integral of Born-Oppenheimer occupation number,
// defined in sim_MCEv1_1e.c together with some more options
extern const size_t BO_occupation_quad_integral_num_pointsMAX;


#ifdef _OPENMP
	#include <omp.h>
	
	#define OVERRIDE_OMP_NUM_THREADS (1 && !CLUSTER_MODE)
	
	#if OVERRIDE_OMP_NUM_THREADS
		#define MAX_NUM_THREADS 4
	#endif

	#include "zpotrf_omp.h"
#endif

// USE_CUDA should be 0, probably breaks MCEv2 code (phononic_overlap_matrix on device instead of host)
#define USE_CUDA (!CLUSTER_MODE && 0)
#define USE_CUDA_PROFILE_MARKS (!CLUSTER_MODE && 0)

#if USE_CUDA
	#include <cuda_runtime.h>
	#include <cusolverDn.h>
	
	#define CUDA_INFO_RETURN -42
#endif

#if USE_CUDA_PROFILE_MARKS
	#include <nvToolsExt.h>
	
	#define COLOR_RED 0xff0000
	#define COLOR_GREEN 0x00ff00
	#define COLOR_BLUE 0x0000ff
#endif

#define SQRT2 1.4142135623730950488

#define PRINT_FULL_DEBUG 0
#define PRINT_CLONE 0
#define PRINT_STORE_OBS_DEBUG 0
#define PRINT_MCEv2_CLONING_DEBUG 0
#define PRINT_MCEv2_REDUCTION_DEBUG 0
#define MEASURE_TIME (USE_CUDA || 1)
#define PRINT_BO_PROGRESS 1		// currently just in BO_QUAD_INTEGRAL

#define EMULATE_MTE 0

#include <complex.h>

#include <mkl_vsl.h>

#define MAX(a, b) (((a) > (b))?(a):(b))
#define MIN(a, b) (((a) < (b))?(a):(b))

// global variables for Runge Kutta algorithm, initialized in sim_MCE.c
extern const double complex cplx_zero;
extern const double complex cplx_one;
extern const double complex cplx_neg_one;

// indicates whether the program shall terminate, triggered by SIGTERM or SIGINT
volatile extern int sigterm_received;

struct param_and_buffer {
	// Hamiltonian parameters
	double hopping;
	double omega;
	double gamma;
	
	// system sizes
	size_t* num_configs;
	size_t max_num_configs;
	size_t num_sites;
	size_t num_sites_x;
	size_t num_sites_y;
	
	size_t num_electrons;
	size_t num_times;
	
	// observables
	double complex* restrict vec2_obs_y;
	double* restrict vec2_obs_nph;
	double* restrict vec2_obs_nel;
	double* restrict vec_obs_el_ekin;
	double* restrict vec2_obs_el_ph_inter;
	double* restrict vec2_obs_BO_nel;
	size_t* restrict fullStorage_num_configs;
	double complex* restrict fullStorage_el_coeff;
	double complex* restrict fullStorage_hoz;
	
	// buffers
	double complex* const restrict phononic_overlap_matrix; // exists always, size max_num_configs^2
	double complex* const restrict total_overlap_matrix;    // exists for MCEv2, size max_num_configs^2
	double complex* const restrict vec2_conf_overlap;       // exists always, size max(max_num_configs^2, num_sites^2)
	double complex* const restrict vec2_buffer;             // exists always, size max(max_num_configs, num_threads)*num_sites
	double complex* const restrict vec2_buffer2;            // exists always, size max(max_num_configs, num_threads)*num_sites
	double* const restrict vec2_density;                    // exists always, size max(max_num_configs, num_threads)*num_sites
	
	double* restrict vec_buffer_local;                      // exists iff vec2_obs_BO_nel != NULL, size BO_OCCUPATION_CALC_METHOD == BO_QUAD_INTEGRAL
	                                                        //                                          ? max(num_sites, BO_occupation_quad_integral_num_points)
	                                                        //                                          : num_sites
	double* restrict vec_buffer_local2;                     // exists iff vec2_obs_BO_nel != NULL, size BO_OCCUPATION_CALC_METHOD == BO_QUAD_INTEGRAL
	                                                        //                                          ? max(num_sites^2, BO_occupation_quad_integral_num_points)
	                                                        //                                          : num_sites
	double* restrict vec_buffer_local_sq_occ;               // exists iff vec2_obs_BO_nel != NULL, size num_sites
	VSLStreamStatePtr rng_stream; // pointer to a thread local RNG
	
	#if USE_CUDA
	cusolverDnHandle_t cusolver_handle;
	cuDoubleComplex* const restrict R_overlap_gpu;
	cuDoubleComplex* const restrict work_gpu;
	int Lwork_gpu;
	int* info_gpu;
	cudaError_t cuda_status;
	cusolverStatus_t cusolver_status;
	#endif
	
	// information
	int* info;
	int printProgress;
	int* order; // buffer for order of weigths in config cloning
	
	// cloning, thread local buffers
	int allow_config_delete;
	int allow_config_clone;
	// both exist if num_configs != max_num_configs or vec2_obs_BO_nel != NULL
	double* BO_eigvecs; // also used to fill with hamiltonian, size num_sites^2
	double* tau;                                            // size num_sites
	// basis reduction, shared buffers, both exist always
	size_t* num_configs_to_delete;    // single value
	size_t* configs_to_delete;        // buffer of size max_num_configs
	short*  config_checked;
	size_t* nearest_neighbors;        // buffer of size 2*max_num_configs
	
	
	#if USE_CUDA
	// runtimes of derivative function with and without cuda
	// to automatically choose the faster one. only timed once, init with 0.
	double* runtime_cuda;
	double* runtime_cpu;
	#endif
	
};
#endif // SIM_MCE_DEFS_H
