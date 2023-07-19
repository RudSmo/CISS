#include "sim_MCEv2_1e.h"
#include "sim_MCE_util.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include <mkl.h>

#include <math.h>
#include "zpotrf_omp.h"

/*
copies and/or transforms the MCE2 wavefunction into MCE1 storage format/normalization
to be used before MCE1_store_observables
*/
void MCE2_to_MCE1_wavefunction(struct param_and_buffer* const pb,
                               const double complex* const restrict MCE2_vec2_hoz,
                               const double complex* const restrict MCE2_vec2_el_coeff,
                               const double complex* const restrict MCE2_vec_coupling_coeff,
                               double complex* const restrict MCE1_vec2_hoz,
                               double complex* const restrict MCE1_vec2_el_coeff) {
	#pragma omp for
	for(size_t n = 0 ; n < pb->num_sites; n++) {
		cblas_zcopy(*pb->num_configs,
		            &MCE2_vec2_hoz[n], pb->num_sites,
		            &MCE1_vec2_hoz[n*pb->max_num_configs], 1);
		
		#if CLUSTER_MODE
			for(size_t k = 0; k < *pb->num_configs; k++) {
				MCE1_vec2_el_coeff[k + n*pb->max_num_configs] = MCE2_vec2_el_coeff[n + k*pb->num_sites]*MCE2_vec_coupling_coeff[k];
			}
		#else
			vzMulI(*pb->num_configs,
			       (MKL_Complex16*) MCE2_vec_coupling_coeff, 1,
			       (MKL_Complex16*) &MCE2_vec2_el_coeff[n], pb->num_sites,
			       (MKL_Complex16*) &MCE1_vec2_el_coeff[n*pb->max_num_configs], 1);
		#endif
	}
}


static double inline cabs2(const double complex z) {
	return fma(creal(z), creal(z), cimag(z)*cimag(z));
}

static void MCE2_fill_overlap_matrix(struct param_and_buffer* const pb,
                                     const double complex* const restrict vec2_hoz,
                                     const double complex* const restrict vec2_el_coeff,
                                     const int also_fill_lower_half) {
	// phononic and total overlap matrix
	// buffer[k] = \sum_n |z_kn|^2  (buffer is here treated as double buffer)
	#pragma omp for
	for(size_t k = 0; k < *pb->num_configs; k++) {
		((double*) pb->vec2_buffer)[k] = cabs2(vec2_hoz[k*pb->num_sites]);
		for(size_t n = 1; n < pb->num_sites; n++) {
			((double*) pb->vec2_buffer)[k] += cabs2(vec2_hoz[k*pb->num_sites + n]);
		}
	}
	
	// fill phononic and total overlap matrix (only upper triangle including main diagonal)
	#pragma omp for schedule(dynamic)
	for(size_t j = 0; j < *pb->num_configs; j++) {
		for(size_t k = 0; k < j; k++) { // only upper triangle excluding main diagonal
			double complex exponent;
			// exponent = \sum_l conj(z_kl)*z_jl
			cblas_zdotc_sub(pb->num_sites,
			                &vec2_hoz[k*pb->num_sites], 1,
			                &vec2_hoz[j*pb->num_sites], 1,
			                &exponent);
			
			// exponent -= .5*\sum_l (|z_kl|^2 + |z_jl|^2) (sums stored in vec2_buffer)
			exponent -= .5*(((double*) pb->vec2_buffer)[k] + ((double*) pb->vec2_buffer)[j]);
			const double complex phononic_overlap = cexp(exponent);
			pb->phononic_overlap_matrix[k + j*pb->max_num_configs] = phononic_overlap;
			
			double complex electronic_overlap;
			cblas_zdotc_sub(pb->num_sites,
			                &vec2_el_coeff[k*pb->num_sites], 1,
			                &vec2_el_coeff[j*pb->num_sites], 1,
			                &electronic_overlap);
			
			pb->total_overlap_matrix[k + j*pb->max_num_configs] = electronic_overlap*phononic_overlap;
		}
		// main diag = 1
		pb->phononic_overlap_matrix[j + j*pb->max_num_configs] = 1.;
		pb->total_overlap_matrix[j + j*pb->max_num_configs] = 1.;
	}
	if(also_fill_lower_half) {
		#pragma omp for schedule(dynamic)
		for(size_t j = 0; j < *pb->num_configs; j++) {
			// fill lower triangle with complex conjugate of upper triangle
			for(size_t k = j+1; k < *pb->num_configs; k++) {
				pb->phononic_overlap_matrix[k + j*pb->max_num_configs] = conj(pb->phononic_overlap_matrix[j + k*pb->max_num_configs]);
				pb->total_overlap_matrix[k + j*pb->max_num_configs] = conj(pb->total_overlap_matrix[j + k*pb->max_num_configs]);
			}
		}
	}
	
	#if PRINT_FULL_DEBUG
		#pragma omp single
			printf("phononic overlap matrix:\n");
		MCE1_print_overlap_matrix(pb->phononic_overlap_matrix, *pb->num_configs, pb->max_num_configs);
		#pragma omp single
		{
			printf("full buffer content:\n");
			for(size_t k = 0; k < *pb->num_configs; k++) {
				for(size_t j = 0; j < *pb->num_configs; j++) {
					printf("%g+i%g ",
					creal(pb->phononic_overlap_matrix[k + j*pb->max_num_configs]),
					cimag(pb->phononic_overlap_matrix[k + j*pb->max_num_configs]));
				}
				printf("\n");
			}
		}
	#endif
	
}

// assumes that phononic and total overlap matrices are already filled
static void MCE2_fill_eom_rhs(struct param_and_buffer* const pb,
                              const double complex* const restrict vec2_hoz,
                              const double complex* const restrict vec2_el_coeff,
                              const double complex* restrict vec_coupling_coeff,
                              double complex* restrict vec_d_coupling_coeff) {
	
	// evaluate right hand side of equation of motion of dD_k/dt
	// TODO hermitean symmetry NOT exploited, unfriendly to parallelization
	#pragma omp for
	for(size_t k = 0; k < *pb->num_configs; k++) { // only upper triangle excluding main diagonal
		vec_d_coupling_coeff[k] = 0.;
		for(size_t j = 0; j < *pb->num_configs; j++) {
			if(j == k) {
				continue; // the diagonal element always evaluates as 0
			}
			/*
			const double complex phononic_overlap = pb->phononic_overlap_matrix[k+j*pb->max_num_configs];
			const double complex total_overlap = pb->total_overlap_matrix[k+j*pb->max_num_configs];
			*/
			const double complex phononic_overlap = k <= j
			                                        ? pb->phononic_overlap_matrix[k+j*pb->max_num_configs]
			                                        : conj(pb->phononic_overlap_matrix[j+k*pb->max_num_configs]);
			const double complex total_overlap = k <= j
			                                     ? pb->total_overlap_matrix[k+j*pb->max_num_configs]
			                                     : conj(pb->total_overlap_matrix[j+k*pb->max_num_configs]);
			
			double complex config_dist_factor = 0.;
			for(size_t n = 0; n < pb->num_sites; n++) {
				const double complex temp = conj(vec2_hoz[j*pb->num_sites+n]-vec2_hoz[k*pb->num_sites+n])
				                            *vec2_el_coeff[j*pb->num_sites+n];
				config_dist_factor += (phononic_overlap*conj(vec2_el_coeff[k*pb->num_sites+n])
				                       -  total_overlap*conj(vec2_el_coeff[j*pb->num_sites+n]))
				                      *temp;
			}
			
			vec_d_coupling_coeff[k] += config_dist_factor*vec_coupling_coeff[j];
		}
		vec_d_coupling_coeff[k] *= -I*pb->gamma;
	}
	
	//const int thread_ID = omp_get_thread_num();
}


void MCE2_joint_derivative_for_rkf(struct param_and_buffer* const pb,
                                   const double complex* restrict vec2_hoz,
                                   const double complex* restrict vec2_el_coeff,
                                   const double complex* restrict vec_coupling_coeff,
                                   double complex* restrict vec2_d_hoz,
                                   double complex* restrict vec2_d_el_coeff,
                                   double complex* restrict vec_d_coupling_coeff) {
	
	#if USE_CUDA_PROFILE_MARKS
		nvtxEventAttributes_t eventAttrib = {0};
		eventAttrib.version = NVTX_VERSION;
		eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
		eventAttrib.colorType = NVTX_COLOR_ARGB;
		eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;
		
		eventAttrib.color = COLOR_RED;
		eventAttrib.message.ascii = "Ehrf";
		nvtxRangePushEx(&eventAttrib);
	#endif
	#pragma omp for
	for(size_t k = 0; k < *pb->num_configs; k++) {
		// hoz
		for(size_t n = 0; n < pb->num_sites; n++) {
			const size_t index = k*pb->num_sites + n;
			pb->vec2_density[index] = cabs2(vec2_el_coeff[index]);
			vec2_d_hoz[index] = -I*pb->omega*vec2_hoz[index] +I*pb->gamma*pb->vec2_density[index];
		}
		
		// el_coeff (without shift of expected interaction energy)
		double half_mean_energy_el_ph = 0.;
		for(size_t n_y = 0; n_y < pb->num_sites_y; n_y++) {
			for(size_t n_x = 0; n_x < pb->num_sites_x; n_x++) {
				const size_t index = k*pb->num_sites + n_y*pb->num_sites_x + n_x;
				
				double complex neighbor_el_coeff = 0.;
				// TODO PBC?
				if(n_x != 0) {
					neighbor_el_coeff += vec2_el_coeff[index - 1];
				}
				if(n_x+1 != pb->num_sites_x) {
					neighbor_el_coeff += vec2_el_coeff[index + 1];
				}
				if(n_y != 0) {
					neighbor_el_coeff += vec2_el_coeff[index - pb->num_sites_x];
				}
				if(n_y+1 != pb->num_sites_y) {
					neighbor_el_coeff += vec2_el_coeff[index + pb->num_sites_x];
				}
				vec2_d_el_coeff[index] = I*(  pb->hopping*neighbor_el_coeff 
				                            + pb->gamma*2*creal(vec2_hoz[index])*vec2_el_coeff[index]);
				half_mean_energy_el_ph += pb->vec2_density[index]*creal(vec2_hoz[index]);
			}
		}
		half_mean_energy_el_ph *= pb->gamma;
		for(size_t n = 0; n < pb->num_sites; n++) {
			const size_t index = k*pb->num_sites + n;
			
			vec2_d_el_coeff[index] -= I*half_mean_energy_el_ph*vec2_el_coeff[index];
		}
	}
	
	#if USE_CUDA_PROFILE_MARKS
		nvtxRangePop();
		eventAttrib.color = COLOR_GREEN;
		eventAttrib.message.ascii = "coupling";
		nvtxRangePushEx(&eventAttrib);
	#endif
	if(pb->gamma == 0. || EMULATE_MTE) {
		#pragma omp single
			memset(vec_d_coupling_coeff, 0, *pb->num_configs*sizeof(double complex));
	} else {
		// coupling coefficients
		// fill right hand side of its eqaution of motion and the overlap_matrix
		#if USE_CUDA_PROFILE_MARKS
			eventAttrib.color = COLOR_RED;
			eventAttrib.message.ascii = "fill overlap_matrix";
			nvtxRangePushEx(&eventAttrib);
		#endif
		MCE2_fill_overlap_matrix(pb, vec2_hoz, vec2_el_coeff, 0);
		#if USE_CUDA_PROFILE_MARKS
			nvtxRangePop();
			eventAttrib.color = COLOR_GREEN;
			eventAttrib.message.ascii = "fill eom_rhs";
			nvtxRangePushEx(&eventAttrib);
		#endif
		MCE2_fill_eom_rhs(pb, vec2_hoz, vec2_el_coeff, vec_coupling_coeff, vec_d_coupling_coeff);
		#if USE_CUDA_PROFILE_MARKS
			nvtxRangePop();
			eventAttrib.color = COLOR_BLUE;
			eventAttrib.message.ascii = "potrf";
			nvtxRangePushEx(&eventAttrib);
		#endif
		
		#if MCEv2_REFINE_SOLUTION
			// copy total overlap matrix and right hand side for later solution refinement
			#pragma omp single
			{
				#pragma omp task
					cblas_zcopy(*pb->num_configs, vec_d_coupling_coeff, 1, pb->vec2_buffer, 1);
				for(size_t i = 0; i < *pb->num_configs; i++) {
					#pragma omp task
						cblas_zcopy(*pb->num_configs, &pb->total_overlap_matrix[i*pb->max_num_configs], 1, &pb->phononic_overlap_matrix[i*pb->max_num_configs], 1);
				}
			}
		#endif
		
		// potrf = cholesky decomposition
		#if PRINT_FULL_DEBUG
			#pragma omp single
				printf("starting zpotrf\n");
		#endif
		#ifdef _OPENMP
			zpotrf_omp('U', *pb->num_configs, pb->total_overlap_matrix, pb->max_num_configs, pb->info);
		#else
			*pb->info =
			LAPACKE_zpotrf(LAPACK_COL_MAJOR, 'U', *pb->num_configs,
			               (lapack_complex_double*) pb->total_overlap_matrix, pb->max_num_configs);
		#endif
		
		// potrs = solving equation
		if(!*pb->info) {
			#if USE_CUDA_PROFILE_MARKS
				nvtxRangePop();
				eventAttrib.color = COLOR_GREEN;
				eventAttrib.message.ascii = "potrs";
				nvtxRangePushEx(&eventAttrib);
			#endif
			#pragma omp single
			{
				#if PRINT_FULL_DEBUG
					printf("starting zpotrs\n");
				#endif
				
				*pb->info =
				LAPACKE_zpotrs(LAPACK_COL_MAJOR, 'U', *pb->num_configs, 1,
				               (const lapack_complex_double*) pb->total_overlap_matrix, pb->max_num_configs,
				               (lapack_complex_double*) vec_d_coupling_coeff, pb->max_num_configs);
				
				#if MCEv2_REFINE_SOLUTION
					double berr, ferr;
					LAPACKE_zporfs(LAPACK_COL_MAJOR,'U', *pb->num_configs, 1,
					               (const lapack_complex_double*) pb->phononic_overlap_matrix, pb->max_num_configs,
					               (const lapack_complex_double*) pb->total_overlap_matrix, pb->max_num_configs,
					               (const lapack_complex_double*) pb->vec2_buffer, pb->max_num_configs,
					               (lapack_complex_double*) vec_d_coupling_coeff, pb->max_num_configs,
					               &ferr, &berr);
					//printf("zporfs: %g\t%g\n", ferr, berr);
				#endif
			}
		}
		#if USE_CUDA_PROFILE_MARKS
			nvtxRangePop();
		#endif
	}
	#if USE_CUDA_PROFILE_MARKS
		nvtxRangePop();
	#endif
	if(*pb->info) { // if solving failed -> abort
		print_lapack_solve_error(*pb->num_configs, pb->max_num_configs, 1,
		                         pb->total_overlap_matrix, vec_d_coupling_coeff,
		                         vec2_hoz, vec2_el_coeff, *pb->info);
	}
	#if PRINT_FULL_DEBUG
		#pragma omp single
			printf("ending derivative\n");
	#endif
}

#if ALLOW_CLONING_MCEv2
static inline void config_mark_to_delete(struct param_and_buffer* pb, const size_t config) {
	#pragma omp critical(mce2_config_mark_to_delete)
		pb->configs_to_delete[(*pb->num_configs_to_delete)++] = config;
	
	#if (PRINT_MCEv2_REDUCTION_DEBUG >=2)
		#pragma omp critical(print)
			printf("marking %lu for deletion\n", config);
	#endif
}


/*
THIS FUNCTION MUST BE WRAPPED IN A #pragma omp single

orders the configurations in a way such that the largest number of entries fullfills
pb->configs_to_delete[i] = *pb->num_configs-i-1
which corresponds to a fixed point regarding the swapping routine
*/
static inline void order_configs_to_delete(struct param_and_buffer* pb) {
	#if (PRINT_MCEv2_REDUCTION_DEBUG >= 3)
		printf("\nordering configs:");
	#endif
	for(size_t i = 0; i < *pb->num_configs_to_delete; ) {
		#if (PRINT_MCEv2_REDUCTION_DEBUG >= 3)
			for(size_t j = 0; j < *pb->num_configs_to_delete; j++) {
				printf("%lu  ", pb->configs_to_delete[j]);
			}
			printf("\n");
		#endif
		const size_t target_ind = *pb->num_configs - pb->configs_to_delete[i] - 1;
		if(target_ind < *pb->num_configs_to_delete && target_ind != i) {
			size_t swapper = pb->configs_to_delete[target_ind];
			pb->configs_to_delete[target_ind] = pb->configs_to_delete[i];
			pb->configs_to_delete[i] = swapper;
		} else {
			i++;
		}
	}
	#if (PRINT_MCEv2_REDUCTION_DEBUG >= 3)
		for(size_t j = 0; j < *pb->num_configs_to_delete; j++) {
			printf("%lu  ", pb->configs_to_delete[j]);
		}
		printf("\n");
	#endif
}


/*
before calling this function, order_configs_to_delete must be called, because that order is assumed

this functions swaps the configurations according to pb->configs_to_delete, such that these are at the end
of the configuration list.
the list to delete is updated as well as hoz, el_coeff, coupling and both overlap matrices
*/
static void swap_configurations_to_delete(struct param_and_buffer* pb,
                                          double complex* restrict vec2_hoz,
                                          double complex* restrict vec2_el_coeff,
                                          double complex* restrict vec_coupling_coeff) {
	#if (PRINT_MCEv2_REDUCTION_DEBUG >= 2)
		#pragma omp single
		{
			printf("swapping, initial delete list:\n");
			for(size_t j = 0; j < *pb->num_configs_to_delete; j++) {
				printf("%lu  ", pb->configs_to_delete[j]);
			}
			printf("\n");
		}
	#endif
	#pragma omp for schedule(static, 1)
	for(size_t i = 0; i < *pb->num_configs_to_delete; i++) {
		const size_t target_ind = *pb->num_configs - i - 1;
		
		// might already be in right spot, if not, swap all related buffers
		if(pb->configs_to_delete[i] != target_ind) {
			// state buffers
			cblas_zswap(pb->num_sites,
			            &vec2_hoz[target_ind*pb->num_sites], 1,
			            &vec2_hoz[pb->configs_to_delete[i]*pb->num_sites], 1);
			cblas_zswap(pb->num_sites,
			            &vec2_el_coeff[target_ind*pb->num_sites], 1,
			            &vec2_el_coeff[pb->configs_to_delete[i]*pb->num_sites], 1);
			const double complex swapper = vec_coupling_coeff[pb->configs_to_delete[i]];
			vec_coupling_coeff[pb->configs_to_delete[i]] = vec_coupling_coeff[target_ind];
			vec_coupling_coeff[target_ind] = swapper;
			
			// columns of overlap matrices
			cblas_zswap(*pb->num_configs,
			            &pb->total_overlap_matrix[target_ind*pb->max_num_configs], 1,
			            &pb->total_overlap_matrix[pb->configs_to_delete[i]*pb->max_num_configs], 1);
			cblas_zswap(*pb->num_configs,
			            &pb->phononic_overlap_matrix[target_ind*pb->max_num_configs], 1,
			            &pb->phononic_overlap_matrix[pb->configs_to_delete[i]*pb->max_num_configs], 1);
			
			// don't update num_configs_to_delete yet, because matrices rows are still missing
		}
	}
	
	#pragma omp for schedule(static, 1)
	for(size_t i = 0; i < *pb->num_configs_to_delete; i++) {
		const size_t target_ind = *pb->num_configs - i - 1;
		// might already be in right spot, if not, swap all related buffers
		if(pb->configs_to_delete[i] != target_ind) {
			// rows of overlap matrices
			cblas_zswap(*pb->num_configs,
			            &pb->total_overlap_matrix[target_ind], pb->max_num_configs,
			            &pb->total_overlap_matrix[pb->configs_to_delete[i]], pb->max_num_configs);
			cblas_zswap(*pb->num_configs,
			            &pb->phononic_overlap_matrix[target_ind], pb->max_num_configs,
			            &pb->phononic_overlap_matrix[pb->configs_to_delete[i]], pb->max_num_configs);
			
			pb->configs_to_delete[i] = target_ind;
		}
	}
	#if (PRINT_MCEv2_REDUCTION_DEBUG >= 2)
		#pragma omp single
		{
			printf("swapping, final delete list:\n");
			for(size_t j = 0; j < *pb->num_configs_to_delete; j++) {
				printf("%lu  ", pb->configs_to_delete[j]);
			}
			printf("\n");
		}
	#endif
}


// pb->total_overlap_matrix must be filled before calling this function
static void fill_nearest_neighbors(struct param_and_buffer* pb, const double max_overlap_abs2) {
	#pragma omp for
	for(size_t j = 0; j < *pb->num_configs-*pb->num_configs_to_delete; j++) {
		
		size_t max_weight_ind = j;
		double max_weight_val = 0.;
		size_t high_weight_ind = j; // shall become second largest element
		double high_weight_val = 0.;
		
		for(size_t k = 0; k < *pb->num_configs-*pb->num_configs_to_delete; k++) {
			if(k != j) {
				const double weight = cabs2(pb->total_overlap_matrix[k + j*pb->max_num_configs]);
				if(max_weight_ind == j || weight >= max_weight_val) {
					high_weight_ind = max_weight_ind;
					high_weight_val = max_weight_val;
					max_weight_ind = k;
					max_weight_val = weight;
				}
				else if(high_weight_ind == k || weight >= high_weight_val) {
					high_weight_ind = k;
					high_weight_val = weight;
				}
			}
		}
		if(max_weight_val >= max_overlap_abs2) {
			pb->nearest_neighbors[2*j+0] = max_weight_ind;
		} else {
			pb->nearest_neighbors[2*j+0] = pb->max_num_configs; // signal that the value is uncritical
		}
		if(high_weight_val >= max_overlap_abs2) {
			pb->nearest_neighbors[2*j+1] = high_weight_ind;
		} else {
			pb->nearest_neighbors[2*j+1] = pb->max_num_configs; // signal that the value is uncritical
		}
	}
	
	#if (PRINT_MCEv2_REDUCTION_DEBUG >= 2)
		#pragma omp single
		{
			printf("critical nearest_neighbors=\n");
			for(size_t j = 0; j < *pb->num_configs-*pb->num_configs_to_delete; j++) {
				if(pb->nearest_neighbors[2*j] != pb->max_num_configs) {
					const double max1_weight = cabs2(pb->total_overlap_matrix[pb->nearest_neighbors[2*j  ] + j*pb->max_num_configs]);
					const double max2_weight = cabs2(pb->total_overlap_matrix[pb->nearest_neighbors[2*j+1] + j*pb->max_num_configs]);
					printf("%2lu: %2lu %2lu, %g %g\n", j, pb->nearest_neighbors[2*j], pb->nearest_neighbors[2*j+1], max1_weight, max2_weight);
				}
			}
			printf("\n");
		}
	#endif
}

// rec_depth must be 0 when calling this function from outside this function
// returns whether an unsolved overlap has beend detected and a new check of all configuratiosn is required
static int check_config_max_overlap(struct param_and_buffer* pb, const size_t config, int rec_depth) {
	if(!pb->config_checked[config]) {
		pb->config_checked[config] = 1;
		// check subcritical overlap
		if(pb->nearest_neighbors[2*config + 0] == pb->max_num_configs) {
			return 0;
		} else if(pb->nearest_neighbors[2*pb->nearest_neighbors[2*config+0] + 0] == config) {
			// config is nearest neighbor of its nearest neighbor -> shortest link found
			pb->config_checked[pb->nearest_neighbors[2*config]] = 1;
			
			// delete that one with larger second largest overlap
			size_t config_to_delete;
			if(pb->nearest_neighbors[2*config + 1] == pb->max_num_configs) {
				config_to_delete = pb->nearest_neighbors[2*config];
			} else if(pb->nearest_neighbors[2*pb->nearest_neighbors[2*config] + 1] == pb->max_num_configs) {
				config_to_delete = config;
			} else {
				// second largest weights
				const double weight2_this  = cabs2(pb->total_overlap_matrix[
				                                      pb->nearest_neighbors[2*config+1]
				                                      + config*pb->max_num_configs]);
				const double weight2_other = cabs2(pb->total_overlap_matrix[
				                                      pb->nearest_neighbors[2*pb->nearest_neighbors[2*config]+1]
				                                      + pb->nearest_neighbors[2*config]*pb->max_num_configs]);
				if(weight2_this > weight2_other) {
					config_to_delete = config;
				} else {
					config_to_delete = pb->nearest_neighbors[2*config];
				}
			}
			config_mark_to_delete(pb, config_to_delete);
			// return at end of function is used
		} else {
			// shortest link is not symmetric... propagated through chain
			return check_config_max_overlap(pb, pb->nearest_neighbors[2*config], rec_depth+1);
		}
		
	//	check_config_max_overlap(pb, j);
	}
	return rec_depth >= 2;
}



void MCE2_reduce_configs(struct param_and_buffer* pb,
                         double complex* restrict vec2_hoz,
                         double complex* restrict vec2_el_coeff,
                         double complex* restrict vec_coupling_coeff){
	 // TODO how to scale with num_configs?
	//const double min_config_weight_abs2 = 1e-2/(*pb->num_configs);
	const double min_config_weight_abs2 = 1e-4/MAX(50, *pb->num_configs);
	const double max_overlap_abs2 = 0.98;
	#if PRINT_MCEv2_REDUCTION_DEBUG
		#pragma omp single
		{
			printf("attempt basis reduction, min_config_weight_abs2=%e", min_config_weight_abs2);
		}
	#endif
	// first: find all configurations with a too small weight
	#pragma omp single nowait // synchronize in following function call
		*pb->num_configs_to_delete = 0;
	MCE2_fill_overlap_matrix(pb, vec2_hoz, vec2_el_coeff, 1);
	
	#pragma omp single
	{
		#if (PRINT_MCEv2_REDUCTION_DEBUG>=2)
			printf("\nweights=");
		#endif
		for(size_t j = 0; j < *pb->num_configs; j++) {
			// check config weight
			const double config_weight_abs2 = cabs2(vec_coupling_coeff[j]);
			#if (PRINT_MCEv2_REDUCTION_DEBUG>=2)
				printf("%g  ", config_weight_abs2);
			#endif
			if(config_weight_abs2 < min_config_weight_abs2) {
				pb->configs_to_delete[(*pb->num_configs_to_delete)++] = j;
			}
		}
		order_configs_to_delete(pb);
		#if PRINT_MCEv2_REDUCTION_DEBUG
			printf("\ndelete configs unordered:\n");
			for(size_t j = 0; j < *pb->num_configs_to_delete; j++) {
				printf("%lu  ", pb->configs_to_delete[j]);
			}
			printf("\n");
		#endif
	}
	
	//const size_t num_configs_after_weigth_delete = *pb->num_configs - *pb->num_configs_to_delete;
	swap_configurations_to_delete(pb, vec2_hoz, vec2_el_coeff, vec_coupling_coeff);
	#if PRINT_MCEv2_REDUCTION_DEBUG
		#pragma omp single
		{
			printf("reduction: finished first swap\nsearching large overlap\n");
		}
	#endif
	
	// second: find configurations with too large overlap between each other in the remaining set
	if(*pb->num_configs-*pb->num_configs_to_delete <= 2) {
		#pragma omp single
		{
			if(   *pb->num_configs-*pb->num_configs_to_delete == 2
			   && cabs2(pb->total_overlap_matrix[1]) > max_overlap_abs2) {
				pb->configs_to_delete[(*pb->num_configs_to_delete)++] = 1;
			}
		}
		#if PRINT_MCEv2_REDUCTION_DEBUG
			#pragma omp single
			{
				printf("reduction: at most 2 configs left\n");
			}
		#endif
	} else {
		short need_to_check = 1;
		while(need_to_check) {
			#if PRINT_MCEv2_REDUCTION_DEBUG
				#pragma omp single
				{
					printf("reduction: checking for too close configs\n");
				}
			#endif
			need_to_check = 0;
			#pragma omp single nowait // synch in next function call
				memset(pb->config_checked, 0, *pb->num_configs*sizeof(short));
			fill_nearest_neighbors(pb, max_overlap_abs2);
			
			const size_t old_num_to_delete = *pb->num_configs_to_delete;
			#pragma omp for schedule(dynamic)
			for(size_t j = 0; j < *pb->num_configs - old_num_to_delete; j++) {
				if(!pb->config_checked[j]) {
					check_config_max_overlap(pb, j, 0);
					// TODO update need to check to allow multiple runs
				}
			}
			#if PRINT_MCEv2_REDUCTION_DEBUG
				#pragma omp single
				{
					printf("delete configs unordered:\n");
					for(size_t j = 0; j < *pb->num_configs_to_delete; j++) {
						printf("%lu  ", pb->configs_to_delete[j]);
					}
					printf("\n");
				}
			#endif
			#pragma omp single
			{
				order_configs_to_delete(pb);
			}
			swap_configurations_to_delete(pb, vec2_hoz, vec2_el_coeff, vec_coupling_coeff);
		}
	}
	
	#if PRINT_MCEv2_REDUCTION_DEBUG
		#pragma omp single
		{
			printf("delete configs   ordered:\n");
			for(size_t j = 0; j < *pb->num_configs_to_delete; j++) {
				printf("%lu  ", pb->configs_to_delete[j]);
			}
			printf("\n");
		}
	#endif
	
	
	// project unneeded configurations onto remaining subset
	if(*pb->num_configs_to_delete != 0) {
		#pragma omp single
		{
			*pb->num_configs -= *pb->num_configs_to_delete;
		}
		zpotrf_omp('U', *pb->num_configs, pb->total_overlap_matrix, pb->max_num_configs, pb->info);
		
		if(*pb->info) {
			#pragma omp single nowait
				fprintf(stderr, "lapack error in MCE2_reduce_configs, zpotrf returned %i\n", *pb->info);
			return;
		}
		// after swapping, pb->configs_to_delete is a decending list of
		// [*pb->num_configs-1, ..., *pb->num_configs-*pb->num_configs_to_delete-1)
		// and *pb->num_configs_to_delete configurations from the end are deleted
		// proect them onto reduced basis set and then forget them
		#pragma omp for
		for(size_t i = 0; i < *pb->num_configs_to_delete; i++) {
			if(!pb->info) {
				const size_t j = *pb->num_configs + i;
				int info =
				LAPACKE_zpotrs(LAPACK_COL_MAJOR, 'U', *pb->num_configs, 1,
				               (const lapack_complex_double*) pb->total_overlap_matrix, pb->max_num_configs,
				               (lapack_complex_double*) &pb->total_overlap_matrix[0 + j*pb->max_num_configs], pb->max_num_configs);
				
				#pragma omp critical(mce2_config_reduction)
				{
					if(info) {
						fprintf(stderr, "lapack error in MCE2_reduce_configs, zpotrs returned %i\n", info);
						*pb->info = info;
					} else {
						cblas_zaxpy(*pb->num_configs, &vec_coupling_coeff[j],
						            &pb->total_overlap_matrix[j*pb->max_num_configs], 1,
						            vec_coupling_coeff, 1);
					}
				}
			}
		}
	}
	
}


double norm2_force_diff_configs(const size_t length,
                                const double* const restrict a, size_t lda,
                                const double* const restrict b, size_t ldb) {
	double diff2 = 0.;
	for(size_t i =0; i < length; i++) {
		double diff_i = fma(-1., a[i*lda]*a[i*lda], b[i*lda]*b[i*lda]); // b^2 - a^2
		diff2 += diff_i*diff_i;
	}
	return diff2;
}


void MCE2_clone_configs(struct param_and_buffer* pb,
                        double complex* restrict vec2_hoz,
                        double complex* restrict vec2_el_coeff,
                        double complex* restrict vec_coupling_coeff) {
	// TODO paramaters only seem to work in special cases, scaling?
	const double min_config_weight = 1./MAX(50, *pb->num_configs);
	// minimal difference between the two selected eigenstates for cloning
	// force diff NOT scaled by gamma -> maximum achievable value is 1/2
	const double config_usage = *pb->num_configs/(double)pb->max_num_configs;
	const double breaking_weighted_force_diff_norm2 = 0.2*config_usage*config_usage;
	
	
	// phonon nieghborhood cutoff, radius of one phonon is 1 or 1/sqrt(2) TODO which one?
	// to convert to x space, multiply this value by sqrt(2)
	const double neighbor_radius = 0.1;
	// number of tolerated neighbors, curr_config is not counted
	const int max_num_neighbors = ((int) pb->num_sites) - 1;
	
	if(*pb->num_configs == pb->max_num_configs || pb->num_sites < 2 || breaking_weighted_force_diff_norm2 >= 0.5) {
		// no free space for clones, nothing to do here
		// cloning needs to split electronic wavefunction -> at least 2 sites
		// force difference need to exceeds breaking_weighted_force_diff_norm2, but is at the same time limited to [0, 1/2]
		return;
	}
	
	// TODO also check neighborhood to avoid similiar states?
	// would require to fill phononic_overlap_matrix here
	// update: phononic overlap_matrix is still filled from previous call to MCE2_reduce_configs
	// total_overlap_matrix is destroyed
	// update2: delete migth not be executed before calling this. stick to direct distance calculation
	
	#ifdef _OPENMP
		const int thread_ID = omp_get_thread_num();
	#else
		const int thread_ID = 0;
	#endif
	#if PRINT_MCEv2_CLONING_DEBUG
		printf("attempt cloning\n");
	#endif
	const size_t num_configs_old = *pb->num_configs;
	#pragma omp for schedule(dynamic)
	for(size_t curr_config = 0; curr_config < num_configs_old; curr_config++) {
		if(*pb->num_configs < pb->max_num_configs && cabs(vec_coupling_coeff[curr_config]) > min_config_weight) {
			// TODO also check derivative coupling?
			
			// count number of neighbors
			int num_neighbors = 0;
			// TODO how about concurrent changes by other cloning ? ignored for now by just checking old number
			for(size_t c = 0; c < num_configs_old && num_neighbors <= max_num_neighbors; c++ ) {
				if(c != curr_config) {
					vdSubI(pb->num_sites,
					       (double*) &vec2_hoz[c*pb->num_sites], 2,
					       (double*) &vec2_hoz[curr_config*pb->num_sites], 2,
					       (double*) &pb->vec2_buffer[thread_ID*pb->num_sites], 1);
					const double dist = cblas_dnrm2(pb->num_sites, (double*) &pb->vec2_buffer[thread_ID*pb->num_sites], 1);
					if(dist < neighbor_radius) {
						num_neighbors++;
					}
				}
			}
			
			#if PRINT_MCEv2_CLONING_DEBUG
			if(num_neighbors > max_num_neighbors) {
				printf("%i: %lu bad  neighbor\n",  thread_ID, curr_config);
			} else {
				printf("%i: %lu good neighbor, check forces\n",  thread_ID, curr_config);
			#else
			if(num_neighbors <= max_num_neighbors) {
			#endif
				// calculate eigenvectors of curr_config-th configuration into pb->BO_eigvecs (col major),
				// eigenvalues not needed (written to pb->vec_buffer)
				/*
				printf("old values:\nhoz: ");
				for(size_t n = 0; n < pb->num_sites; n++) {
					printf("%g+i%g  ",
					       creal(vec2_hoz[curr_config*pb->num_sites + n]),
					       cimag(vec2_hoz[curr_config*pb->num_sites + n]));
				}
				printf("\n");
				*/
				int info = MCE_fill_BO_basis(pb, -2*pb->gamma, (double*) vec2_hoz, 2, 2*curr_config*pb->num_sites, 0);
				
				#if (PRINT_MCEv2_CLONING_DEBUG >= 2)
					#pragma omp critical(print)
					{
						printf("BO_eigvecs=[\n");
						for(size_t i = 0; i < pb->num_sites; i++) {
							printf("[");
							for(size_t j = 0; j < pb->num_sites; j++) {
								printf("%g  ", pb->BO_eigvecs[i + j*pb->num_sites]);
							}
							printf("],\n");
						}
						printf("]\n\n");
					}
				#endif
				
				
				if(info) {
					fprintf(stderr, "info=%i error in MCE2_clone_configs after MCE_fill_BO_basis, curr_config=%lu\n",
						    info, curr_config);
				}
				
				double complex* const restrict vec_BO_coeff = &pb->vec2_buffer[thread_ID*pb->num_sites];
				double* const restrict vec_BO_prob = pb->tau;
				
				// project electronic state onto calculated BO-eigenbasis
				// real part
				cblas_dgemv(CblasColMajor, CblasConjTrans,
				            pb->num_sites, pb->num_sites, 1., pb->BO_eigvecs, pb->num_sites,
				            ((double*) &vec2_el_coeff[curr_config*pb->num_sites]) + 0, 2,
				            0., ((double*) vec_BO_coeff) + 0, 2);
				
				// imag part
				cblas_dgemv(CblasColMajor, CblasConjTrans,
				            pb->num_sites, pb->num_sites, 1., pb->BO_eigvecs, pb->num_sites,
				            ((double*) &vec2_el_coeff[curr_config*pb->num_sites]) + 1, 2,
				            0., ((double*) vec_BO_coeff) + 1, 2);
				
				
				// calculate BO occupation probability and search largest and second largest element
				size_t max_weight_ind = 1;
				size_t high_weight_ind = 1; // shall become second largest element
				vec_BO_prob[1] = 0.;
				for(size_t i = 0; i < pb->num_sites; i++) {
					vec_BO_prob[i] = cabs2(vec_BO_coeff[i]);
					if(vec_BO_prob[i] >= vec_BO_prob[max_weight_ind]) {
						high_weight_ind = max_weight_ind;
						max_weight_ind = i;
					}
					else if(vec_BO_prob[i] >= vec_BO_prob[high_weight_ind]) {
						high_weight_ind = i;
					}
				}
				/*
				printf("curr_config=%lu\nweight_ind, max=%lu, high=%lu\n", curr_config, max_weight_ind, high_weight_ind);
				printf("all:\n");
				for(size_t i = 0; i < pb->num_sites; i++) {
					printf("%g  ", vec_BO_prob[i]);
				}
				printf("\n");
				fflush(stdout);
				*/
				// TODO density of max occupied one might be buffered?
				
				size_t max_weighted_force_diff_ind = high_weight_ind;
				double max_weighted_force_diff_val = norm2_force_diff_configs(pb->num_sites,
				                                                              &pb->BO_eigvecs[max_weight_ind*pb->num_sites], 1,
				                                                              &pb->BO_eigvecs[high_weight_ind*pb->num_sites], 1);
				max_weighted_force_diff_val *= vec_BO_prob[high_weight_ind];
				
				for(size_t i = 0; i < pb->num_sites; i++) {
					// TODO vec_BO_prob[i]*sqrt(2) is upper bound? can save norm calculation
					if(i != max_weighted_force_diff_ind && i != max_weight_ind) {
						// calc force diff
						double force_diff = norm2_force_diff_configs(pb->num_sites,
						                                             &pb->BO_eigvecs[max_weight_ind*pb->num_sites], 1,
						                                             &pb->BO_eigvecs[i*pb->num_sites], 1);
						force_diff *= vec_BO_prob[i];
						
						if(force_diff > max_weighted_force_diff_val) {
							max_weighted_force_diff_ind = i;
							max_weighted_force_diff_val = force_diff;
						}
						
					}
				}
				max_weighted_force_diff_val *= vec_BO_prob[max_weight_ind];
				
				
				#if PRINT_MCEv2_CLONING_DEBUG
				if(max_weighted_force_diff_val <= breaking_weighted_force_diff_norm2) {
					printf("%i: %lu reject force_diff=%g\n",  thread_ID, curr_config, max_weighted_force_diff_val);
				}
				else {
					printf("%i: %lu accept force_diff=%g\n",  thread_ID, curr_config, max_weighted_force_diff_val);
				#else
				if(max_weighted_force_diff_val > breaking_weighted_force_diff_norm2) {
				#endif
					// clone, attempt to get a target number
					size_t target_config;
					#pragma omp critical(MCEv2_clone_target_config)
					{
						target_config = *pb->num_configs;
						if(target_config < pb->max_num_configs) {
							*pb->num_configs += 1;
						}
					}
					// check whether target number is valid
					if(target_config < pb->max_num_configs) {
						
						#if PRINT_MCEv2_CLONING_DEBUG
							printf("%i: %lu clone to %lu, main states=%lu, %lu\n",
							       thread_ID, curr_config, target_config,
							       max_weight_ind, max_weighted_force_diff_ind);
						#endif
						#if (PRINT_MCEv2_CLONING_DEBUG >= 2)
							#pragma omp critical(print)
							{
							printf("old values:\nhoz: ");
							for(size_t n = 0; n < pb->num_sites; n++) {
								printf("%g+i%g  ",
								       creal(vec2_hoz[curr_config*pb->num_sites + n]),
								       cimag(vec2_hoz[curr_config*pb->num_sites + n]));
							}
							printf("\nelc: ");
							for(size_t n = 0; n < pb->num_sites; n++) {
								printf("%g+i%g  ",
								       creal(vec2_el_coeff[curr_config*pb->num_sites + n]),
								       cimag(vec2_el_coeff[curr_config*pb->num_sites + n]));
							}
							printf("\ncoupling: %g+i%g\n",
							       creal(vec_coupling_coeff[curr_config]),
							       cimag(vec_coupling_coeff[curr_config]));
							}
						#endif
						// copy coupling_coeff and hoz
						vec_coupling_coeff[target_config] = vec_coupling_coeff[curr_config];
						cblas_zcopy(pb->num_sites,
						            &vec2_hoz[curr_config*pb->num_sites], 1,
						            &vec2_hoz[target_config*pb->num_sites], 1);
						
						// split el_coeff into two orthogonal vectors
						// highest occupied BO eigenstate...
						const size_t BO_group1_ref_ind = max_weight_ind;
						// ... and that one with largest weighted force deviation...
						const size_t BO_group2_ref_ind = max_weighted_force_diff_ind;
						// are used to define 2 groups of eigenvectors,
						// each group forms one of the two configurations after cloning
						// all eigenstates are assigened to the groups by largest scalar product of forces on nuclei
						// ignoring harmonic contribution, which is the same for all, and normalizing
						// -> electronic density (site local) scalar product
						
						double* const restrict ref1_density
							= ((double*) &pb->vec2_buffer2[thread_ID*pb->num_sites]) + 0;
						double* const restrict ref2_density
							= ((double*) &pb->vec2_buffer2[thread_ID*pb->num_sites]) + pb->num_sites;
						double* const restrict iterate_density = &pb->vec2_density[thread_ID*pb->num_sites];
						
						vdMul(pb->num_sites,
						      &pb->BO_eigvecs[BO_group1_ref_ind*pb->num_sites],
						      &pb->BO_eigvecs[BO_group1_ref_ind*pb->num_sites],
						      ref1_density);
						vdMul(pb->num_sites,
						      &pb->BO_eigvecs[BO_group2_ref_ind*pb->num_sites],
						      &pb->BO_eigvecs[BO_group2_ref_ind*pb->num_sites],
						      ref2_density);
						
						memset(&vec2_el_coeff[target_config*pb->num_sites], 0, pb->num_sites*sizeof(double complex));
						
						#if (PRINT_MCEv2_CLONING_DEBUG >= 2)
							double norm1 = cblas_dznrm2(pb->num_sites, &vec2_el_coeff[curr_config*pb->num_sites], 1);
							printf("norm before clone=%g\n", norm1);
						#endif
						
						for(size_t i = 0; i < pb->num_sites; i++) {
							int group; // 0 for group1, 1 for group2
							if(i == BO_group1_ref_ind) {
								group = 0;
							} else if(i == BO_group2_ref_ind) {
								group = 1;
							} else {
								vdMul(pb->num_sites,
								      &pb->BO_eigvecs[i*pb->num_sites],
								      &pb->BO_eigvecs[i*pb->num_sites],
								      iterate_density);
								const double dot1 = cblas_ddot(pb->num_sites, iterate_density, 1, ref1_density, 1);
								const double dot2 = cblas_ddot(pb->num_sites, iterate_density, 1, ref2_density, 1);
								
								group = dot2 > dot1;
							}
							
							// split curr_config by copying group2 to target_config
							if(group == 1) {
								// real part
								cblas_daxpy(pb->num_sites, creal(vec_BO_coeff[i]),
								            &pb->BO_eigvecs[i*pb->num_sites], 1,
								            ((double*) &vec2_el_coeff[target_config*pb->num_sites]) + 0, 2);
								// imag part
								cblas_daxpy(pb->num_sites, cimag(vec_BO_coeff[i]),
								            &pb->BO_eigvecs[i*pb->num_sites], 1,
								            ((double*) &vec2_el_coeff[target_config*pb->num_sites]) + 1, 2);
							}
						}
						// to complete split, substract target_config from curr_config
						cblas_zaxpy(pb->num_sites, &cplx_neg_one,
						            &vec2_el_coeff[target_config*pb->num_sites], 1,
						            &vec2_el_coeff[  curr_config*pb->num_sites], 1);
						
						#if (PRINT_MCEv2_CLONING_DEBUG >= 2)
							double norm_check = cblas_dznrm2(pb->num_sites, &vec2_el_coeff[curr_config*pb->num_sites], 1);
							printf("\nnew config1 norm = %g\n", norm_check);
							norm_check = cblas_dznrm2(pb->num_sites, &vec2_el_coeff[target_config*pb->num_sites], 1);
							printf("new config2 norm = %g\n", norm_check);
						#endif
						
						
						// TODO norm can be calculated easier by summing cabs2(BO coefficients)
						// renormalize both configs involved in cloning and scale coupling_coeff
						double norm = cblas_dznrm2(pb->num_sites, &vec2_el_coeff[target_config*pb->num_sites], 1);
						cblas_zdscal(pb->num_sites, 1./norm, &vec2_el_coeff[target_config*pb->num_sites], 1);
						vec_coupling_coeff[target_config] *= norm;
						
						// norm_c^2 = 1-norm_target^2 since orthogonal
						//norm = cblas_dznrm2(pb->num_sites, &vec2_el_coeff[curr_config*pb->num_sites], 1);
						norm = sqrt(1-norm*norm);
						cblas_zdscal(pb->num_sites, 1./norm, &vec2_el_coeff[curr_config*pb->num_sites], 1);
						vec_coupling_coeff[curr_config] *= norm;
						
						
						#if (PRINT_MCEv2_CLONING_DEBUG >= 2)
							#pragma omp critical(print)
							{
							printf("old config:\nhoz: ");
							for(size_t n = 0; n < pb->num_sites; n++) {
								printf("%g+i%g  ",
								       creal(vec2_hoz[curr_config*pb->num_sites + n]),
								       cimag(vec2_hoz[curr_config*pb->num_sites + n]));
							}
							printf("\nelc: ");
							for(size_t n = 0; n < pb->num_sites; n++) {
								printf("%g+i%g  ",
								       creal(vec2_el_coeff[curr_config*pb->num_sites + n]),
								       cimag(vec2_el_coeff[curr_config*pb->num_sites + n]));
							}
							printf("\nnew config:\nhoz: ");
							for(size_t n = 0; n < pb->num_sites; n++) {
								printf("%g+i%g  ",
								       creal(vec2_hoz[target_config*pb->num_sites + n]),
								       cimag(vec2_hoz[target_config*pb->num_sites + n]));
							}
							printf("\nelc: ");
								for(size_t n = 0; n < pb->num_sites; n++) {
								printf("%g+i%g  ", creal(vec2_el_coeff[target_config*pb->num_sites + n]),
								cimag(vec2_el_coeff[target_config*pb->num_sites + n]));
							}
							printf("\ncoupling: %g+i%g  %g+i%g\n",
							       creal(vec_coupling_coeff[curr_config]), cimag(vec_coupling_coeff[curr_config]),
							       creal(vec_coupling_coeff[target_config]), cimag(vec_coupling_coeff[target_config]));
							//printf("coupling sum: %g+i%g\n",
							//       creal(vec_coupling_coeff[c])+creal(vec_coupling_coeff[target_config]),
							//       cimag(vec_coupling_coeff[c])+cimag(vec_coupling_coeff[target_config]));
							
							double complex dot_check;
							cblas_zdotc_sub(pb->num_sites,
							                &vec2_el_coeff[curr_config*pb->num_sites], 1,
							                &vec2_el_coeff[target_config*pb->num_sites], 1,
							                &dot_check);
							printf("config1^* @ config2 = %g + i%g\n\n", creal(dot_check), cimag(dot_check));
							}
						#endif
					}
				}
			}
		}
	}
}
#endif
