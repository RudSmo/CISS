#include "sim_MCEv1_1e.h"
#include "sim_MCE_util.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include <mkl.h>

#include <math.h>

const size_t BO_occupation_MC_integral_min_num_iterations = 1024; // must be at least 2
const double BO_occupation_MC_intergal_tolerance = 0.01;

// for num_sites <= 2:
const size_t BO_occupation_quad_integral_num_points2 = 40;
const double BO_occupation_quad_integral_erfinv_scale2 = 3.;
// for num_sites > 2:
const size_t BO_occupation_quad_integral_num_pointsL = 30;
const double BO_occupation_quad_integral_erfinv_scaleL = 2.;

const double BO_occupation_quad_integral_erfinv_range = 0.98;

const double BO_occupation_quad_integral_min_weight = 1e-7;

// for malloc in sim_MCE.c
const size_t BO_occupation_quad_integral_num_pointsMAX = MAX(BO_occupation_quad_integral_num_points2,
                                                             BO_occupation_quad_integral_num_pointsL);

/*
calculates all weights(= electron densities within one configuration)
w_jn = |a_jn|^2/\sum_l |a_jl|^2
for a fixed configuration j
weigths for same site in different configurations is stored contiuous (column major)
*/
static void MCE1_fill_density(const size_t num_configs, const size_t max_num_configs,
                              const size_t num_sites,
                              double* restrict vec2_density,
                              const double complex* const restrict vec2_el_coeff) {
	
	#pragma omp for
	for(size_t j=0; j < num_configs; j++) {
		double norm2_site = 0.;
		for(size_t n = 0; n < num_sites; n++) {
			const double complex coeff = vec2_el_coeff[n*max_num_configs + j];
			vec2_density[n*max_num_configs + j] = creal(coeff)*creal(coeff) + cimag(coeff)*cimag(coeff);
			norm2_site += vec2_density[n*max_num_configs + j];
		}
		cblas_dscal(num_sites, 1./norm2_site, &vec2_density[j], max_num_configs);
	}
}


/*
computes the overlap matrix of the phononic states between configurations
only the upper triangle and main diagonal are filled (CblasUpper)
CblasColMajor is used
the diagonal is by definition always 1 (overlap of a gaussian with itself)
*/
static void MCE1_fill_overlap_matrix(double complex* const restrict phononic_overlap_matrix,
                                   const size_t num_configs, const size_t max_num_configs,
                                   const size_t num_sites,
                                   const double complex* const restrict vec2_hoz,
                                   double complex* const restrict vec2_buffer) {
	// buffer[k] = \sum_n |z_kn|^2  (buffer is here treated as double buffer)
	#pragma omp for
	for(size_t k = 0; k < num_configs; k++) {
		double complex z = vec2_hoz[k];
		((double*) vec2_buffer)[k] = creal(z)*creal(z) + cimag(z)*cimag(z);
		for(size_t n = 1; n < num_sites; n++) {
			z = vec2_hoz[n*max_num_configs+k];
			((double*) vec2_buffer)[k] += creal(z)*creal(z) + cimag(z)*cimag(z);
		}
	}
	
	#pragma omp for schedule(dynamic)
	for(size_t j = 0; j < num_configs; j++) {
		for(size_t k = 0; k < j; k++) { // only upper triangle excluding main diagonal
			double complex exponent;
			// exponent = \sum_l conj(z_kl)*z_jl
			cblas_zdotc_sub(num_sites,
			                &vec2_hoz[k], max_num_configs,
			                &vec2_hoz[j], max_num_configs,
			                &exponent);
			
			// exponent -= .5*\sum_l (|z_kl|^2 + |z_jl|^2) (sums stored in vec2_buffer)
			exponent -= .5*(((double*) vec2_buffer)[k] + ((double*) vec2_buffer)[j]);
			#if EMULATE_MTE
				phononic_overlap_matrix[k + j*max_num_configs] = 0; //cexp(exponent);
			#else
				phononic_overlap_matrix[k + j*max_num_configs] = cexp(exponent);
			#endif
		}
		phononic_overlap_matrix[j*max_num_configs + j] = 1.; // main diag = 1
	}
}

/*
fills the right and side (rhs) of the equation of motion(eom), such that
R @ \dot a = rhs
the derivative can be obtained by solving the system of equations

rhs is a column major matrix with multiple right hand sides
(each column corresponds to one site n)
vec2_density and (phononic) phononic_overlap_matrix must be filled before calling this function
vec2_buffer must be set to 0 before calling this function
vec2_conf_overlap must be set to 0 before calling this function
*/
static void MCE1_fill_eom_rhs(double complex* const restrict vec2_eom_rhs,
                              const size_t num_configs,
                              const size_t max_num_configs,
                              const size_t num_sites_x,
                              const size_t num_sites_y,
                              const double hopping,
                              const double omega,
                              const double gamma,
                              const double complex* const restrict phononic_overlap_matrix,
                              const double complex* const restrict vec2_hoz,
                              const double complex* const restrict vec2_el_coeff,
                              const double* const restrict vec2_density,
                              double complex* const restrict vec2_buffer,
                              double complex* const restrict vec2_buffer2,
                              double complex* const restrict vec2_conf_overlap) {
	const size_t num_sites= num_sites_x*num_sites_y;
	
	const double complex imag_hopping = hopping*I;
	const double complex neg_imag_gamma = -gamma*I;
	const double complex imag_gamma = gamma*I;
	//const double complex neg_imag = -I;
	
	#if USE_CUDA_PROFILE_MARKS
		nvtxEventAttributes_t eventAttrib = {0};
		eventAttrib.version = NVTX_VERSION;
		eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
		eventAttrib.colorType = NVTX_COLOR_ARGB;
		eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;
		
		eventAttrib.color = COLOR_RED;
		eventAttrib.message.ascii = "1";
		nvtxRangePushEx(&eventAttrib);
	#endif
	
	// kinetic part (includes setting vec2_eom_rhs to zero)
	//&el_coeff[y, x, c] == &el_coeff[(y*num_sites_x + x)*max_num_configs + c]
	#pragma omp for collapse(2) schedule(dynamic)
	for(size_t y = 0; y < num_sites_y; y++) {
		for(size_t x = 0; x < num_sites_x; x++) {
			if(x != num_sites_x - 1) {
				cblas_zaxpy(num_configs, &cplx_one,
				            &vec2_el_coeff[(y*num_sites_x + x+1)*max_num_configs], 1,
				            &vec2_buffer  [(y*num_sites_x + x  )*max_num_configs], 1);
			}
			if(x != 0) {
				cblas_zaxpy(num_configs, &cplx_one,
				            &vec2_el_coeff[(y*num_sites_x + x-1)*max_num_configs], 1,
				            &vec2_buffer  [(y*num_sites_x + x  )*max_num_configs], 1);
			}
			if(y != num_sites_y - 1) {
				cblas_zaxpy(num_configs, &cplx_one,
				            &vec2_el_coeff[((y+1)*num_sites_x + x)*max_num_configs], 1,
				            &vec2_buffer  [( y   *num_sites_x + x)*max_num_configs], 1);
			}
			if(y != 0) {
				cblas_zaxpy(num_configs, &cplx_one,
				            &vec2_el_coeff[((y-1)*num_sites_x + x)*max_num_configs], 1,
				            &vec2_buffer  [( y   *num_sites_x + x)*max_num_configs], 1);
			}
		}
	}
	
	#pragma omp for schedule(static)
	for(size_t n = 0; n < num_sites; n++) {
		cblas_zhemv(CblasColMajor, CblasUpper,
		            num_configs, &imag_hopping,
		            phononic_overlap_matrix, max_num_configs,
		            &vec2_buffer[n*max_num_configs], 1,
		            &cplx_zero, &vec2_eom_rhs[n*max_num_configs], 1);
	}
	
	#if USE_CUDA_PROFILE_MARKS
		nvtxRangePop();
		eventAttrib.color = COLOR_GREEN;
		eventAttrib.message.ascii = "2";
		nvtxRangePushEx(&eventAttrib);
	#endif
	
	
	#ifdef _OPENMP
		const size_t thread_ID = omp_get_thread_num();
	#else
		const size_t thread_ID = 0;
	#endif
	
	// coupling part
	#pragma omp for schedule(static) nowait
	for(size_t n = 0; n < num_sites; n++) {
		// vec2_buffer_jn = a_jn*z_jn
		vzMul(num_configs,
		      (MKL_Complex16*) &vec2_el_coeff[n*max_num_configs],
		      (MKL_Complex16*) &vec2_hoz[n*max_num_configs],
		      (MKL_Complex16*) &vec2_buffer[thread_ID*max_num_configs]);
		
		// vec2_eom_rhs -= gamma * R @ vec2_buffer
		cblas_zhemv(CblasColMajor, CblasUpper,
		            num_configs, &imag_gamma,
		            phononic_overlap_matrix, max_num_configs,
		            &vec2_buffer[thread_ID*max_num_configs], 1,
		            &cplx_one, &vec2_eom_rhs[n*max_num_configs], 1);
		
		// vec2_buffer_n = - gamma * R @ a_.n
		cblas_zhemv(CblasColMajor, CblasUpper,
		            num_configs, &imag_gamma,
		            phononic_overlap_matrix, max_num_configs,
		            &vec2_el_coeff[n*max_num_configs], 1,
		            &cplx_zero, &vec2_buffer[thread_ID*max_num_configs], 1);
		
		// vec2_buffer2 = vec2_buffer * conj(vec2_hoz)
		vzMulByConj(num_configs,
		            (MKL_Complex16*) &vec2_buffer[thread_ID*max_num_configs],
		            (MKL_Complex16*) &vec2_hoz[n*max_num_configs],
		            (MKL_Complex16*) &vec2_buffer2[thread_ID*max_num_configs]);
		
		// vec2_eom_rhs += vec2_buffer2
		cblas_daxpy(num_configs*2, 1., (double*) &vec2_buffer2[thread_ID*max_num_configs], 1,
		            (double*) &vec2_eom_rhs[n*max_num_configs], 1);
	}
	
	
	
	#if USE_CUDA_PROFILE_MARKS
		nvtxRangePop();
		eventAttrib.color = COLOR_BLUE;
		eventAttrib.message.ascii = "3";
		nvtxRangePushEx(&eventAttrib);
	#endif
	
	#pragma omp for schedule(static)
	for(size_t j = 0; j < num_configs; j++) {
		double complex* const restrict overlap_column = &vec2_conf_overlap[j*max_num_configs];
		memset(overlap_column, 0, max_num_configs*sizeof(double complex));
		
		// vec2_conf_overlap[k, j] = sum_l z_kl^* *d_jl, k hidden in daxpy
		for(size_t l = 0; l < num_sites; l++) {
			// real part
			cblas_daxpy(num_configs, vec2_density[j + l*max_num_configs],
			            ((double*) &vec2_hoz[l*max_num_configs]) + 0, 2,
			            ((double*) overlap_column) + 0, 2);
			
			// imag part, negative for comjugate
			cblas_daxpy(num_configs, -vec2_density[j + l*max_num_configs],
			            ((double*) &vec2_hoz[l*max_num_configs]) + 1, 2,
			            ((double*) overlap_column) + 1, 2);
		}
		for(size_t k = 0; k < num_configs; k++) {
			const double complex R_kj = k<j
			                            ? phononic_overlap_matrix[k + j*max_num_configs]
			                            : conj(phononic_overlap_matrix[j + k*max_num_configs]);
			overlap_column[k] *= R_kj;
		}
	}
	
	#pragma omp for
	for(size_t n = 0; n < num_sites; n++) {
		cblas_zgemv(CblasColMajor, CblasNoTrans,
		            num_configs, num_configs, &neg_imag_gamma,
		            vec2_conf_overlap, max_num_configs,
		            &vec2_el_coeff[n*max_num_configs], 1,
		            &cplx_one, &vec2_eom_rhs[n*max_num_configs], 1);
	}
	
	
	
	#if USE_CUDA_PROFILE_MARKS
		nvtxRangePop();
		eventAttrib.color = COLOR_GREEN;
		eventAttrib.message.ascii = "4";
		nvtxRangePushEx(&eventAttrib);
	#endif
	
	#pragma omp for
	for(size_t j = 0; j < num_configs; j++) {
		// buffer_jl = \sum_l density_jl * i* imag(z_jl)
		vec2_buffer[j] = 0.;
		for(size_t l = 0; l < num_sites; l++) {
			vec2_buffer[j] +=  I*vec2_density[l*max_num_configs + j]
			                  *cimag(vec2_hoz[l*max_num_configs + j]);
		}
	}
	
	
	#pragma omp for
	for(size_t n = 0; n < num_sites; n++) {
		// buffer2_jn = a_jn*buffer_j0
		vzMul(num_configs,
		      (MKL_Complex16*) &vec2_el_coeff[n*max_num_configs],
		      (MKL_Complex16*) vec2_buffer,
		      (MKL_Complex16*) &vec2_buffer2[n*max_num_configs]);
		
		// eom_rhs_.n += R @ buffer2
		cblas_zhemv(CblasColMajor, CblasUpper,
		            num_configs, &neg_imag_gamma,
		            phononic_overlap_matrix, max_num_configs,
		            &vec2_buffer2[n*max_num_configs], 1,
		            &cplx_one, &vec2_eom_rhs[n*max_num_configs], 1);
	}
	
	#if USE_CUDA_PROFILE_MARKS
		nvtxRangePop();
	#endif
}

__attribute__((unused))
static void MCE1_calc_BO_occupation_by_mid_each(struct param_and_buffer* const pb,
                                                const double norm2,
                                                double* const restrict vec2_obs_BO_nel,
                                                const double complex* const restrict vec2_hoz,
                                                const double complex* const restrict vec2_el_coeff) {
	int info;
	#pragma omp for
	for(size_t k = 0; k < *pb->num_configs; k++) {
		#if PRINT_STORE_OBS_DEBUG
			printf("\nconfiguration k=%lu\nRe(z)=", k);
			for(size_t n = 0; n < pb->num_sites; n++) {
				printf("%g ", creal(vec2_hoz[k + n*pb->max_num_configs]));
			}
			printf("\n");
		#endif
		
		info = MCE_fill_BO_basis(pb, -2*pb->gamma, (double*) vec2_hoz, 2*pb->max_num_configs, 2*k, 1);
		
		if(info) {
			fprintf(stderr, "LAPACKE_dsteqr returned %i in observable storage of BO eigenstates!\n"
			                "BO coefficients set to 0.\n", info);
			cblas_dscal(pb->num_sites, 0., ((double*) &pb->vec2_buffer2[k]) + 0, 2*pb->max_num_configs);
			cblas_dscal(pb->num_sites, 0., ((double*) &pb->vec2_buffer2[k]) + 1, 2*pb->max_num_configs);
		} else {
			// project to eigenbasis
			// real part
			cblas_dgemv(CblasColMajor, CblasTrans, pb->num_sites, pb->num_sites,
			            1., pb->BO_eigvecs, pb->num_sites, ((double*) &vec2_el_coeff[k]) + 0, 2*pb->max_num_configs,
			            0., ((double*) &pb->vec2_buffer2[k]) + 0, 2*pb->max_num_configs);
			// imag part
			cblas_dgemv(CblasColMajor, CblasTrans, pb->num_sites, pb->num_sites,
			            1., pb->BO_eigvecs, pb->num_sites, ((double*) &vec2_el_coeff[k]) + 1, 2*pb->max_num_configs,
			            0., ((double*) &pb->vec2_buffer2[k] + 1), 2*pb->max_num_configs);
		}
		#if PRINT_STORE_OBS_DEBUG
			printf("eigenvectors:\n");
			for(size_t n = 0; n < pb->num_sites; n++) {
				for(size_t m = 0; m < pb->num_sites; m++) {
					printf("%g ", creal(pb->BO_eigvecs[n + m*pb->num_sites]));
				}
				printf("\n");
			}
		#endif
		
	} // for loop, projection to basis
	
	
	
	#if PRINT_STORE_OBS_DEBUG
		printf("\nlocal states:\n");
		for(size_t k = 0; k < *pb->num_configs; k++) {
			for(size_t n = 0; n < pb->num_sites; n++) {
				printf("%g+i%g ",
				       creal(vec2_el_coeff[k + n*pb->max_num_configs]),
				       cimag(vec2_el_coeff[k + n*pb->max_num_configs]));
			}
			printf("\n");
		}
		
		printf("\nprojected states:\n");
		for(size_t k = 0; k < *pb->num_configs; k++) {
			for(size_t n = 0; n < pb->num_sites; n++) {
				printf("%g+i%g ",
				       creal(pb->vec2_buffer2[k + n*pb->max_num_configs]),
				       cimag(pb->vec2_buffer2[k + n*pb->max_num_configs]));
			}
			printf("\n");
		}
			
		printf("overlap matrix:\n");
		MCE1_print_overlap_matrix(pb->phononic_overlap_matrix, *pb->num_configs, pb->max_num_configs);
		printf("occupation number:\n");
	#endif
	
	// vec2_buffer2 contains the coefficients of BO_states
	// calculate occupation of each eigenstate
	#pragma omp for schedule(static)
	for(size_t l = 0; l < pb->num_sites; l++) {
		// buffer = R @ e_jl
		cblas_zhemv(CblasColMajor, CblasUpper, *pb->num_configs,
		            &cplx_one, pb->phononic_overlap_matrix, pb->max_num_configs,
		            &pb->vec2_buffer2[l*pb->max_num_configs], 1,
		            &cplx_zero, &pb->vec2_buffer[l*pb->max_num_configs], 1);
		
		double complex z;
		// z = conj(a_kl)^T @ R @ e_jl
		cblas_zdotc_sub(*pb->num_configs,
		                &pb->vec2_buffer2[l*pb->max_num_configs], 1,
		                &pb->vec2_buffer [l*pb->max_num_configs], 1,
		                &z);
		vec2_obs_BO_nel[l] = creal(z);
		#if PRINT_STORE_OBS_DEBUG
			printf("%g+i%g ", creal(z), cimag(z));
		#endif
	}
	#if PRINT_STORE_OBS_DEBUG
		printf("\n");
	#endif
	
	#pragma omp single
		cblas_dscal (pb->num_sites, 1./norm2, vec2_obs_BO_nel, 1);
}


__attribute__((unused))
static void MCE1_calc_BO_occupation_by_mid_pair(struct param_and_buffer* const pb,
                                                const double norm2,
                                                double* const restrict vec2_obs_BO_nel,
                                                const double complex* const restrict vec2_hoz,
                                                const double complex* const restrict vec2_el_coeff) {
	int info;
	#ifdef _OPENMP
		const size_t thread_ID = omp_get_thread_num();
	#else
		const size_t thread_ID = 0;
	#endif
	double* restrict occupation_number_local = &pb->vec2_density[thread_ID*pb->num_sites];
	// clear buffer for occupation number on each site, collected over threads after loop
	memset(occupation_number_local, 0, pb->num_sites*sizeof(double));
	
	#pragma omp for schedule(dynamic) nowait
	for(int j = *pb->num_configs-1; j >= 0; j--) {
		for(size_t k = 0; k <= j; k++) {
			#if CLUSTER_MODE
				for(size_t n = 0; n < pb->num_sites; n++) {
					((double*) &pb->vec2_buffer2[thread_ID*pb->num_sites])[n]
					    =   creal(vec2_hoz[k + n*pb->max_num_configs])
					      + creal(vec2_hoz[j + n*pb->max_num_configs]);
				}
			#else
				vdAddI(pb->num_sites,
				       (double*) &vec2_hoz[k], 2*pb->max_num_configs,
				       (double*) &vec2_hoz[j], 2*pb->max_num_configs,
				       (double*) &pb->vec2_buffer2[thread_ID*pb->num_sites], 1);
			#endif
			
			// no factor 2 for gamma, because that is contained in sum
			info = MCE_fill_BO_basis(pb, -pb->gamma, (double*) &pb->vec2_buffer2[thread_ID*pb->num_sites], 1, 0, 1);
			if(info) {
				fprintf(stderr, "calculation of mid_pair BO occupation failed with k=%lu, j=%i, info=%i\n", k, j, info);
			}
			// project to eigenbasis, config k -> buffer
			// real part
			cblas_dgemv(CblasColMajor, CblasTrans, pb->num_sites, pb->num_sites,
			            1., pb->BO_eigvecs, pb->num_sites, ((double*) &vec2_el_coeff[k]) + 0, 2*pb->max_num_configs,
			            0., ((double*) &pb->vec2_buffer[thread_ID*pb->num_sites]) + 0, 2);
			// imag part
			cblas_dgemv(CblasColMajor, CblasTrans, pb->num_sites, pb->num_sites,
			            1., pb->BO_eigvecs, pb->num_sites, ((double*) &vec2_el_coeff[k]) + 1, 2*pb->max_num_configs,
			            0., ((double*) &pb->vec2_buffer[thread_ID*pb->num_sites] + 1), 2);
			
			// project to eigenbasis, config j -> buffer2
			// real part
			cblas_dgemv(CblasColMajor, CblasTrans, pb->num_sites, pb->num_sites,
			            1., pb->BO_eigvecs, pb->num_sites, ((double*) &vec2_el_coeff[j]) + 0, 2*pb->max_num_configs,
			            0., ((double*) &pb->vec2_buffer2[thread_ID*pb->num_sites]) + 0, 2);
			// imag part
			cblas_dgemv(CblasColMajor, CblasTrans, pb->num_sites, pb->num_sites,
			            1., pb->BO_eigvecs, pb->num_sites, ((double*) &vec2_el_coeff[j]) + 1, 2*pb->max_num_configs,
			            0., ((double*) &pb->vec2_buffer2[thread_ID*pb->num_sites] + 1), 2);
			
			
			// electronic overlap
			for(size_t n = 0; n < pb->num_sites; n++) {
				double complex electronic_overlap =  conj(pb->vec2_buffer [thread_ID*pb->num_sites+n])
				                                    *     pb->vec2_buffer2[thread_ID*pb->num_sites+n];
				
				if(k == j) {
					occupation_number_local[n] += creal(electronic_overlap); // phononic_overlap == 1
				} else {
					occupation_number_local[n] += 2*creal(electronic_overlap*pb->phononic_overlap_matrix[k+j*pb->max_num_configs]);
				}
			}
		}
	}
	
	#pragma omp critical
	{
		cblas_daxpy(pb->num_sites, 1./norm2,
		            occupation_number_local, 1,
		            vec2_obs_BO_nel, 1);
	}
	#pragma omp barrier
}


__attribute__((unused))
static void MCE1_calc_BO_occupation_by_MC_integral(struct param_and_buffer* const pb,
                                                   const double norm2,
                                                   double* const restrict vec2_obs_BO_nel,
                                                   const double complex* const restrict vec2_hoz,
                                                   const double complex* const restrict vec2_el_coeff) {
	int info;
	#ifdef _OPENMP
		const size_t thread_ID = omp_get_thread_num();
	#else
		const size_t thread_ID = 0;
	#endif
	double* restrict occupation_number_local = &pb->vec2_density[thread_ID*pb->num_sites];
	
	const double scaled_var_tolerance = BO_occupation_MC_intergal_tolerance*BO_occupation_MC_intergal_tolerance/((*pb->num_configs*(*pb->num_configs+1))/2);
	
	#pragma omp for schedule(dynamic) nowait
	for(int j = *pb->num_configs-1; j >= 0; j--) {
		for(size_t k = 0; k <= j; k++) {
			#if CLUSTER_MODE
				for(size_t n = 0; n < pb->num_sites; n++) {
					pb->vec_buffer_local[n]
					    =   creal(vec2_hoz[k + n*pb->max_num_configs])
					      + creal(vec2_hoz[j + n*pb->max_num_configs]);
				}
			#else
				vdAddI(pb->num_sites,
				       (double*) &vec2_hoz[k], 2*pb->max_num_configs,
				       (double*) &vec2_hoz[j], 2*pb->max_num_configs,
				       pb->vec_buffer_local, 1);
			#endif
			
			double complex exponent_space_indenpendent = 0.;
			for(size_t n = 0; n < pb->num_sites; n++) {
				const double temp = creal(vec2_hoz[k + n*pb->max_num_configs])-creal(vec2_hoz[j + n*pb->max_num_configs]);
				exponent_space_indenpendent += -.5*temp*temp;
				exponent_space_indenpendent += I*(  creal(vec2_hoz[k + n*pb->max_num_configs])*cimag(vec2_hoz[k + n*pb->max_num_configs])
				                                  - creal(vec2_hoz[j + n*pb->max_num_configs])*cimag(vec2_hoz[j + n*pb->max_num_configs]));
			}
			
			// reset sum and sum of squared occupation number
			memset(occupation_number_local, 0, pb->num_sites*sizeof(double));
			memset(pb->vec_buffer_local_sq_occ, 0, pb->num_sites*sizeof(double));
			double current_var = 2*scaled_var_tolerance; // anything larger than scaled_var_tolerance
			size_t num_iterations = 0;
			while((num_iterations < BO_occupation_MC_integral_min_num_iterations || current_var > scaled_var_tolerance) && !sigterm_received) { // TODO ERROR_ESTIMATE_TOO_LARGE
				num_iterations++;
				vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, pb->rng_stream, pb->num_sites, pb->vec_buffer_local2, 0. , 1.);
				double complex exponent = exponent_space_indenpendent;
				
				for(size_t n = 0; n < pb->num_sites; n++) {
					// shift space to center point
					pb->vec_buffer_local2[n] += pb->vec_buffer_local[n];
					
					exponent += I*(cimag(vec2_hoz[j + n*pb->max_num_configs]) - cimag(vec2_hoz[k + n*pb->max_num_configs]))
					             *pb->vec_buffer_local2[n];
					
				}
				exponent = cexp(exponent);
				
				// calc eigenstates
				// no factor 2 for gamma, because space is rescaled by that factor
				info = MCE_fill_BO_basis(pb, -pb->gamma, pb->vec_buffer_local2, 1, 0, 1);
				if(info) {
					fprintf(stderr, "calculation of MC integral BO occupation failed with k=%lu, j=%i, info=%i\n", k, j, info);
				} else {
					// project to eigenbasis, config k -> buffer
					// real part
					cblas_dgemv(CblasColMajor, CblasTrans, pb->num_sites, pb->num_sites,
					            1., pb->BO_eigvecs, pb->num_sites, ((double*) &vec2_el_coeff[k]) + 0, 2*pb->max_num_configs,
					            0., ((double*) &pb->vec2_buffer[thread_ID*pb->num_sites]) + 0, 2);
					// imag part
					cblas_dgemv(CblasColMajor, CblasTrans, pb->num_sites, pb->num_sites,
					            1., pb->BO_eigvecs, pb->num_sites, ((double*) &vec2_el_coeff[k]) + 1, 2*pb->max_num_configs,
					            0., ((double*) &pb->vec2_buffer[thread_ID*pb->num_sites] + 1), 2);
					
					// project to eigenbasis, config j -> buffer2
					// real part
					cblas_dgemv(CblasColMajor, CblasTrans, pb->num_sites, pb->num_sites,
					            1., pb->BO_eigvecs, pb->num_sites, ((double*) &vec2_el_coeff[j]) + 0, 2*pb->max_num_configs,
					            0., ((double*) &pb->vec2_buffer2[thread_ID*pb->num_sites]) + 0, 2);
					// imag part
					cblas_dgemv(CblasColMajor, CblasTrans, pb->num_sites, pb->num_sites,
					            1., pb->BO_eigvecs, pb->num_sites, ((double*) &vec2_el_coeff[j]) + 1, 2*pb->max_num_configs,
					            0., ((double*) &pb->vec2_buffer2[thread_ID*pb->num_sites] + 1), 2);
					
					
					current_var = 0;
					
					// electronic overlap
					for(size_t n = 0; n < pb->num_sites; n++) {
						double complex electronic_overlap =  conj(pb->vec2_buffer [thread_ID*pb->num_sites+n])
						                                    *     pb->vec2_buffer2[thread_ID*pb->num_sites+n];
						//printf("eo: %g+i%g\n", creal(electronic_overlap), cimag(electronic_overlap)); 
						const double integrand = creal(electronic_overlap*exponent)
						                         *(k == j? 1. : 2.)/norm2;
						occupation_number_local[n] += integrand;
						pb->vec_buffer_local_sq_occ[n] += integrand*integrand;
						
						// estimate current variance of mean as (<x^2> - <x>^2)/N
						const double current_mean = occupation_number_local[n]/num_iterations;
						const double current_site_var = (  pb->vec_buffer_local_sq_occ[n]/num_iterations
						                                 - current_mean*current_mean)/num_iterations;
						if(current_site_var > current_var) {
							current_var = current_site_var;
						}
					}
				}
			}
			cblas_dscal(pb->num_sites, 1./num_iterations, occupation_number_local, 1);
			
			#pragma omp critical(MCE1_calc_BO_occupation_by_MC_integral_collect)
			{
				cblas_daxpy(pb->num_sites, 1.,
				            occupation_number_local, 1,
				            vec2_obs_BO_nel, 1);
			}
		}
	}
	#pragma omp barrier
}

// x_grid, pb->vec_buffer_local2 are assumed to be disjoint buffers
// &pb->vec_buffer_local2[dim*pb->num_sites] (+pb->num_sites slots) contain the function output
// pb->vec_buffer_local_sq_occ is overwritten with positions in the grid
__attribute__((unused))
static void MCE1_calc_BO_occupation_by_quad_integral_rec_iterator(const size_t dim,
                                                                  const size_t j, const size_t k,
                                                                  struct param_and_buffer* const pb,
                                                                  const size_t BO_occupation_quad_integral_num_points, // size of x_grid
                                                                  const double* const restrict x_grid,
                                                                  const double complex* const restrict vec2_hoz,
                                                                  const double complex* const restrict vec2_el_coeff,
                                                                  const double complex exponent_space_indenpendent) {
	double* output = &pb->vec_buffer_local2[dim*pb->num_sites];
	memset(output, 0, pb->num_sites*sizeof(double));
	double* x_current = pb->vec_buffer_local_sq_occ;
	
	#ifdef _OPENMP
		const size_t thread_ID = omp_get_thread_num();
	#else
		const size_t thread_ID = 0;
	#endif
	
	for(size_t i = 0; i < BO_occupation_quad_integral_num_points && !sigterm_received; i++) {
		x_current[dim] =   creal(vec2_hoz[k + dim*pb->max_num_configs])
		                 + creal(vec2_hoz[j + dim*pb->max_num_configs])
		                 + x_grid[i];
		double dx;
		if(i == 0) {
			dx = x_grid[1] - x_grid[0];
		} else if(i == BO_occupation_quad_integral_num_points - 1) {
			dx = x_grid[i] - x_grid[i-1];
		} else {
			dx = .5*(x_grid[i+1] - x_grid[i-1]);
		}
		
		if(dim == pb->num_sites-1) {
			// actual calculation of the integrand
			double complex exponent = exponent_space_indenpendent;
				
			for(size_t n = 0; n < pb->num_sites; n++) {
				const double temp = x_current[n] - creal(vec2_hoz[k + n*pb->max_num_configs])
				                                 - creal(vec2_hoz[j + n*pb->max_num_configs]);
				exponent += I*(cimag(vec2_hoz[j + n*pb->max_num_configs]) - cimag(vec2_hoz[k + n*pb->max_num_configs]))
				             *x_current[n]
				            - .5*temp*temp;
			}
			exponent = cexp(exponent);
			
			// calc eigenstates
			// no factor 2 for gamma, because space is rescaled by that factor
			int info = MCE_fill_BO_basis(pb, -pb->gamma, x_current, 1, 0, 1);
			if(info) {
				fprintf(stderr, "calculation of MC integral BO occupation failed with k=%lu, j=%lu, info=%i\n", k, j, info);
			} else {
				// project to eigenbasis, config k -> buffer
				// real part
				cblas_dgemv(CblasColMajor, CblasTrans, pb->num_sites, pb->num_sites,
				            1., pb->BO_eigvecs, pb->num_sites, ((double*) &vec2_el_coeff[k]) + 0, 2*pb->max_num_configs,
				            0., ((double*) &pb->vec2_buffer[thread_ID*pb->num_sites]) + 0, 2);
				// imag part
				cblas_dgemv(CblasColMajor, CblasTrans, pb->num_sites, pb->num_sites,
				            1., pb->BO_eigvecs, pb->num_sites, ((double*) &vec2_el_coeff[k]) + 1, 2*pb->max_num_configs,
				            0., ((double*) &pb->vec2_buffer[thread_ID*pb->num_sites] + 1), 2);
				
				// project to eigenbasis, config j -> buffer2
				// real part
				cblas_dgemv(CblasColMajor, CblasTrans, pb->num_sites, pb->num_sites,
				            1., pb->BO_eigvecs, pb->num_sites, ((double*) &vec2_el_coeff[j]) + 0, 2*pb->max_num_configs,
				            0., ((double*) &pb->vec2_buffer2[thread_ID*pb->num_sites]) + 0, 2);
				// imag part
				cblas_dgemv(CblasColMajor, CblasTrans, pb->num_sites, pb->num_sites,
				            1., pb->BO_eigvecs, pb->num_sites, ((double*) &vec2_el_coeff[j]) + 1, 2*pb->max_num_configs,
				            0., ((double*) &pb->vec2_buffer2[thread_ID*pb->num_sites] + 1), 2);
				
				
				// electronic overlap
				for(size_t n = 0; n < pb->num_sites; n++) {
					double complex electronic_overlap =  conj(pb->vec2_buffer [thread_ID*pb->num_sites+n])
					                                    *     pb->vec2_buffer2[thread_ID*pb->num_sites+n];
					//printf("eo: %g+i%g\n", creal(electronic_overlap), cimag(electronic_overlap)); 
					const double integrand = creal(electronic_overlap*exponent)
					                         *(k == j? 1. : 2.);
					output[n] += dx*integrand;
				}
			}
		} else { // recursion
			MCE1_calc_BO_occupation_by_quad_integral_rec_iterator(dim + 1, j, k, pb,
			                                                      BO_occupation_quad_integral_num_points, x_grid,
			                                                      vec2_hoz, vec2_el_coeff,
			                                                      exponent_space_indenpendent);
			cblas_daxpy(pb->num_sites, dx,
			            &pb->vec_buffer_local2[(dim+1)*pb->num_sites], 1,
			            &pb->vec_buffer_local2[ dim   *pb->num_sites], 1);
		}
		
	}
}

__attribute__((unused))
static void MCE1_calc_BO_occupation_by_quad_integral(struct param_and_buffer* const pb,
                                                     const double norm2,
                                                     double* const restrict vec2_obs_BO_nel,
                                                     const double complex* const restrict vec2_hoz,
                                                     const double complex* const restrict vec2_el_coeff) {
	#ifdef _OPENMP
		const size_t thread_ID = omp_get_thread_num();
	#else
		const size_t thread_ID = 0;
	#endif
	double* restrict occupation_number_local = &pb->vec2_density[thread_ID*pb->num_sites];
	memset(occupation_number_local, 0, pb->num_sites*sizeof(double));
	
	const size_t BO_occupation_quad_integral_num_points = pb->num_sites <= 2
	                                                      ? BO_occupation_quad_integral_num_points2
	                                                      : BO_occupation_quad_integral_num_pointsL;
	const size_t BO_occupation_quad_integral_erfinv_scale = pb->num_sites <= 2
	                                                      ? BO_occupation_quad_integral_erfinv_scale2
	                                                      : BO_occupation_quad_integral_erfinv_scaleL;
	
	// pb->vec_buffer_local = x_grid for quad integral, sampled as in ehrenfest_init.py: scale*linspace(-range, range, num_points)
	for(size_t i = 0; i < BO_occupation_quad_integral_num_points; i++) {
		pb->vec_buffer_local2[i] = BO_occupation_quad_integral_erfinv_range*(-1.+2.*i/(BO_occupation_quad_integral_num_points-1.));
	}
	vdErfInv(BO_occupation_quad_integral_num_points, pb->vec_buffer_local2, pb->vec_buffer_local);
	cblas_dscal(BO_occupation_quad_integral_num_points, BO_occupation_quad_integral_erfinv_scale, pb->vec_buffer_local, 1);
	
	/*
	#pragma omp single
	{
		printf("x_grid:\n");
		for(size_t i = 0; i < BO_occupation_quad_integral_num_points; i++) {
			printf("%g ", pb->vec_buffer_local[i]);
		}
		printf("\n");
	}
	*/
	#if PRINT_BO_PROGRESS
		#pragma omp single nowait
		{
			printf("evaluating BO integral...\n");
		}
	#endif
	
	const double integral_factor = pow(2*M_PI, -.5*pb->num_sites)/norm2;
	
	#pragma omp for schedule(dynamic) nowait
	for(int j = *pb->num_configs-1; j >= 0; j--) {
		for(size_t k = 0; k <= j; k++) {
			double complex exponent_space_indenpendent = 0.;
			for(size_t n = 0; n < pb->num_sites; n++) {
				const double temp = creal(vec2_hoz[k + n*pb->max_num_configs])-creal(vec2_hoz[j + n*pb->max_num_configs]);
				exponent_space_indenpendent += -.5*temp*temp;
				exponent_space_indenpendent += I*(  creal(vec2_hoz[k + n*pb->max_num_configs])*cimag(vec2_hoz[k + n*pb->max_num_configs])
				                                  - creal(vec2_hoz[j + n*pb->max_num_configs])*cimag(vec2_hoz[j + n*pb->max_num_configs]));
			}
			const double conf_j_norm = cblas_dznrm2(pb->num_sites, &vec2_el_coeff[j], pb->max_num_configs);
			const double conf_k_norm = cblas_dznrm2(pb->num_sites, &vec2_el_coeff[k], pb->max_num_configs);
			const double max_int_weight = conf_j_norm*conf_k_norm*exp(creal(exponent_space_indenpendent));
			//printf("%e\n", max_int_weight);
			
			if(max_int_weight > BO_occupation_quad_integral_min_weight) {
				MCE1_calc_BO_occupation_by_quad_integral_rec_iterator(0, j, k, pb,
				                                                      BO_occupation_quad_integral_num_points, pb->vec_buffer_local,
				                                                      vec2_hoz, vec2_el_coeff, exponent_space_indenpendent);
				
				cblas_daxpy(pb->num_sites, integral_factor, pb->vec_buffer_local2, 1, occupation_number_local, 1);
			}
		}
		#if PRINT_BO_PROGRESS
			if(0 == j%MAX(1, *pb->num_configs/20)) {
				#pragma omp critical(print_bo_progress)
				{
					printf("%lu/%lu\n", *pb->num_configs-j, *pb->num_configs);
					fflush(stdout);
				}
			}
		#endif
	}
	
	#pragma omp critical(MCE1_calc_BO_occupation_by_MC_integral_collect)
	{
		cblas_daxpy(pb->num_sites, 1.,
		            occupation_number_local, 1,
		            vec2_obs_BO_nel, 1);
	}
	#pragma omp barrier
}

/*
calculates and stores the observables
 - x_l, p_l (stored as x+ip = y)
 - number of phonons (nph)
 - electron density (nel)
 - kinetic energy of each electron (unscaled)
 - electron-phonon interaction per site (unscaled)

unscaled means that the factor of hopping or gamma was not mutliplied in

all pointers are taken from the passed struct param_and_buffer* const pb, no pointers needs to be shifted in time.
The thread local buffers BO_eigvecs and tau from pb are referenced iff pb->vec2_obs_BO_nel != NULL && time % obs_BO_nel_save_dist == 0
obs_BO_nel_save_dist can be used to reduce the number of evaluations of the BO states.

returns the norm of the wavefunction
*/
double MCE1_store_observables(const size_t time, const size_t obs_BO_nel_save_dist,
                              struct param_and_buffer* const pb,
                              const double complex* const restrict vec2_hoz,
                              const double complex* const restrict vec2_el_coeff) {
	#if PRINT_STORE_OBS_DEBUG
		printf("storing observables, time=%lu\n", time);
	#endif
	double complex* const restrict vec2_obs_y = &pb->vec2_obs_y[time*pb->num_sites];
	double* const restrict vec2_obs_nph = &pb->vec2_obs_nph[time*pb->num_sites];
	double* const restrict vec2_obs_nel = &pb->vec2_obs_nel[time*pb->num_sites];
	double* const restrict vec_obs_el_ekin = &pb->vec_obs_el_ekin[time];
	double* const restrict vec2_obs_el_ph_inter = &pb->vec2_obs_el_ph_inter[time*pb->num_sites];
	double* const restrict vec2_obs_BO_nel = pb->vec2_obs_BO_nel != NULL && (time % obs_BO_nel_save_dist == 0)
	                                         ? &pb->vec2_obs_BO_nel[(time/obs_BO_nel_save_dist)*pb->num_sites]
	                                         : NULL;
	
	
	MCE1_fill_overlap_matrix(pb->phononic_overlap_matrix, *pb->num_configs, pb->max_num_configs, pb->num_sites,
	                       vec2_hoz, pb->vec2_buffer);
	
	#pragma omp for schedule(static) nowait
	for(size_t n = 0; n < pb->num_sites; n++) {
		// buffer = R @ a_n
		cblas_zhemv(CblasColMajor, CblasUpper, *pb->num_configs,
		            &cplx_one, pb->phononic_overlap_matrix, pb->max_num_configs,
		            &vec2_el_coeff[n*pb->max_num_configs], 1,
		            &cplx_zero, &pb->vec2_buffer[n*pb->max_num_configs], 1);
	}
	
	// norm and electron density
	double norm2 = 0.;  // norm2 is thread local, collect at end of function
	#pragma omp for schedule(static)
	for(size_t l = 0; l < pb->num_sites; l++) {
		double complex z;
		// z = conj(a_kl)^T @ R @ a_jl
		cblas_zdotc_sub(*pb->num_configs,
		                &vec2_el_coeff[l*pb->max_num_configs], 1,
		                &pb->vec2_buffer[l*pb->max_num_configs], 1,
		                &z);
		norm2 += creal(z);
		vec2_obs_nel[l] = creal(z);
	}
	
	// kinetic energy of electron (not scaled with t0)
	// R @ a_j_l = buffer
	double e_kin_thread = 0.0;
	#pragma omp for collapse(2)
	for(size_t y = 0; y < pb->num_sites_y; y++) {
		for(size_t x = 0; x < pb->num_sites_x; x++) {
			double complex z;
			
			if(x != pb->num_sites_x - 1) {
				cblas_zdotc_sub(*pb->num_configs,
				                &  vec2_el_coeff[(y*pb->num_sites_x + x+1)*pb->max_num_configs], 1,
				                &pb->vec2_buffer[(y*pb->num_sites_x + x  )*pb->max_num_configs], 1,
				                &z);
				e_kin_thread -= creal(z);
				
				cblas_zdotc_sub(*pb->num_configs,
				                &  vec2_el_coeff[(y*pb->num_sites_x + x  )*pb->max_num_configs], 1,
				                &pb->vec2_buffer[(y*pb->num_sites_x + x+1)*pb->max_num_configs], 1,
				                &z);
				e_kin_thread -= creal(z);
			}
			
			if(y != pb->num_sites_y - 1) {
				cblas_zdotc_sub(*pb->num_configs,
				                &  vec2_el_coeff[((y+1)*pb->num_sites_x + x)*pb->max_num_configs], 1,
				                &pb->vec2_buffer[( y   *pb->num_sites_x + x)*pb->max_num_configs], 1,
				                &z);
				e_kin_thread -= creal(z);
				
				cblas_zdotc_sub(*pb->num_configs,
				                &  vec2_el_coeff[( y   *pb->num_sites_x + x)*pb->max_num_configs], 1,
				                &pb->vec2_buffer[((y+1)*pb->num_sites_x + x)*pb->max_num_configs], 1,
				                &z);
				e_kin_thread -= creal(z);
			}
		}
	}
	
	// electron-phonon interaction (not scaled with gamma)
	#pragma omp for schedule(static) nowait
	for(size_t l = 0; l < pb->num_sites; l++) {
		
		vzMul(*pb->num_configs,
		      (MKL_Complex16*) &vec2_el_coeff[l*pb->max_num_configs],
		      (MKL_Complex16*) &vec2_hoz[l*pb->max_num_configs],
		      (MKL_Complex16*) &pb->vec2_buffer2[l*pb->max_num_configs]);
		
		double complex z;
		
		// z = conj(a_kl*z_kl = buffer2) @ (R @ a_jl = buffer)
		cblas_zdotc_sub(*pb->num_configs,
		                &pb->vec2_buffer2[l*pb->max_num_configs], 1,
		                &pb->vec2_buffer[l*pb->max_num_configs], 1,
		                &z);
		
		vec2_obs_el_ph_inter[l] = -creal(z);
		
		// buffer = R @ (z_jl * a_jl)
		cblas_zhemv(CblasColMajor, CblasUpper, *pb->num_configs,
		            &cplx_one, pb->phononic_overlap_matrix, pb->max_num_configs,
		            &pb->vec2_buffer2[l*pb->max_num_configs], 1,
		            &cplx_zero, &pb->vec2_buffer[l*pb->max_num_configs], 1);
		
		// z = conj(a_kl) @ (R @ (a_jl*z_jl) = buffer)
		cblas_zdotc_sub(*pb->num_configs,
		                &  vec2_el_coeff[l*pb->max_num_configs], 1,
		                &pb->vec2_buffer[l*pb->max_num_configs], 1,
		                &z);
		
		vec2_obs_el_ph_inter[l] -= creal(z);
	}
	
	// number of phonons and y_l = (x_l, p_l)
	#pragma omp for schedule(static)
	for(size_t l = 0; l < pb->num_sites; l++) {
		for(size_t n = 0; n < pb->num_sites; n++) {
			double complex z;
			// number of phonons (corresponds to phonon energy/omega)
			// buffer_l = a_n*z_l
			vzMul(*pb->num_configs,
			      (MKL_Complex16*) &vec2_hoz[l*pb->max_num_configs],
			      (MKL_Complex16*) &vec2_el_coeff[n*pb->max_num_configs],
			      (MKL_Complex16*) &pb->vec2_buffer[l*pb->max_num_configs]);
			
			// buffer2_l = R @ (buffer_l = a_n*z_l)
			cblas_zhemv(CblasColMajor, CblasUpper, *pb->num_configs,
			            &cplx_one, pb->phononic_overlap_matrix, pb->max_num_configs,
			            &pb->vec2_buffer[l*pb->max_num_configs], 1,
			            &cplx_zero, &pb->vec2_buffer2[l*pb->max_num_configs], 1);
			
			// z = conj(buffer_l = a_n*z_l) @ (R @ (a_n*z_l) = buffer2_l)
			cblas_zdotc_sub(*pb->num_configs,
			                &pb->vec2_buffer [l*pb->max_num_configs], 1,
			                &pb->vec2_buffer2[l*pb->max_num_configs], 1,
			                &z);
			
			vec2_obs_nph[l] += creal(z);
			
			// obs_y
			cblas_zdotc_sub(*pb->num_configs,
			                &  vec2_el_coeff [n*pb->max_num_configs], 1,
			                &pb->vec2_buffer2[l*pb->max_num_configs], 1,
			                &z);
			vec2_obs_y[l] += z;
		}
	}
	
	if(pb->fullStorage_num_configs) {
		#pragma omp single nowait
			pb->fullStorage_num_configs[time] = *pb->num_configs;
	}
	if(pb->fullStorage_el_coeff){
		#pragma omp single nowait
			memcpy(&pb->fullStorage_el_coeff[time*pb->num_sites*pb->max_num_configs],
			       vec2_el_coeff, pb->num_sites*pb->max_num_configs*sizeof(double complex));
	}
	if(pb->fullStorage_hoz) {
		#pragma omp single nowait
			memcpy(&pb->fullStorage_hoz[time*pb->num_sites*pb->max_num_configs],
			       vec2_hoz, pb->num_sites*pb->max_num_configs*sizeof(double complex));
	}
	
	
	
	#pragma omp single 
	{
		((double*) pb->vec2_buffer)[0] = 0.; // store norm2 in buffer
		//printf("ID    e_kin_thread,  sum_e_kin   norm2_thread    norm2_sum\n");
	}
	// collect norm and kinetic energy
	#pragma omp critical(MCE1_store_observables_collect_ekin_norm)
	{
		*vec_obs_el_ekin += e_kin_thread;
		((double*) pb->vec2_buffer)[0] += norm2;
	}
	#pragma omp barrier
	
	norm2 = ((double*) pb->vec2_buffer)[0];
	
	// divide all observables by total squared norm
	#pragma omp sections
	{
		#pragma omp section
			cblas_zdscal(pb->num_sites, SQRT2/norm2, vec2_obs_y,   1);
		#pragma omp section
			cblas_dscal (pb->num_sites,    1./norm2, vec2_obs_nph, 1);
		#pragma omp section
			cblas_dscal (pb->num_sites,    1./norm2, vec2_obs_nel, 1);
		
		#pragma omp section
			cblas_dscal (pb->num_sites,    1./norm2, vec2_obs_el_ph_inter, 1);
		#pragma omp section
			*vec_obs_el_ekin /= norm2;
	}
	
	
	if(vec2_obs_BO_nel != NULL) {
		#if PRINT_STORE_OBS_DEBUG
			printf("storing Born Oppenheimer occupation number\n");
		#endif
		#pragma omp single
		{
			int info = MCE_fill_BO_basis(pb, -sqrt(2.)*pb->gamma, (double*) vec2_obs_y, 2, 0, 0);
			memcpy(pb->vec2_conf_overlap, pb->BO_eigvecs, pb->num_sites*pb->num_sites*sizeof(double));
			if(info) {
				fprintf(stderr, "calculating reference basis of BO eigenstates failed, info=%i\n", info);
				*pb->info = info;
			}
		}
		const int sigterm_received_before_BO_occ = sigterm_received;
		
		if(!*pb->info) {
			#if BO_OCCUPATION_CALC_METHOD == BO_MID_EACH
				MCE1_calc_BO_occupation_by_mid_each(pb, norm2, vec2_obs_BO_nel, vec2_hoz, vec2_el_coeff);
			#elif BO_OCCUPATION_CALC_METHOD == BO_MID_PAIR
				MCE1_calc_BO_occupation_by_mid_pair(pb, norm2, vec2_obs_BO_nel, vec2_hoz, vec2_el_coeff);
			#elif BO_OCCUPATION_CALC_METHOD == BO_MC_INTEGRAL
				MCE1_calc_BO_occupation_by_MC_integral(pb, norm2, vec2_obs_BO_nel, vec2_hoz, vec2_el_coeff);
			#elif BO_OCCUPATION_CALC_METHOD == BO_QUAD_INTEGRAL
				MCE1_calc_BO_occupation_by_quad_integral(pb, norm2, vec2_obs_BO_nel, vec2_hoz, vec2_el_coeff);
			#else
				#error no method for BO eigenstates selected
			#endif
		}
		
		if(sigterm_received != sigterm_received_before_BO_occ) {
			// incomplete calculation of BO occupation number -> ste to 0
			memset(vec2_obs_BO_nel, 0, pb->num_sites*sizeof(double));
		}
		
	} // end of occupation of BO eigenstates
	
	return norm2;
}


/*
calculates the derivatives of phononic (hoz) and electronic corrdinates
returns info from lapack solving in *pb->info, abort simulation if not 0
*/
void MCE1_joint_derivative_for_rkf(struct param_and_buffer* const pb,
                                   const double complex* restrict vec2_hoz_old,
                                   const double complex* restrict vec2_el_coeff,
                                   double complex* restrict vec2_d_hoz,
                                   double complex* restrict vec2_d_el_coeff) {
	#if PRINT_FULL_DEBUG
		#pragma omp single nowait
			printf("starting derivative\n");
	#endif
	
	#pragma omp single nowait
	{
		*pb->info = 0;
	}
	
	#if USE_CUDA_PROFILE_MARKS
		nvtxEventAttributes_t eventAttrib = {0};
		eventAttrib.version = NVTX_VERSION;
		eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
		eventAttrib.colorType = NVTX_COLOR_ARGB;
		eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;
		
		eventAttrib.color = COLOR_GREEN;
		eventAttrib.message.ascii = "derivative";
		nvtxRangePushEx(&eventAttrib);
		
		eventAttrib.color = COLOR_RED;
		eventAttrib.message.ascii = "fill R";
		nvtxRangePushEx(&eventAttrib);
	#endif
	
	MCE1_fill_overlap_matrix(pb->phononic_overlap_matrix, *pb->num_configs, pb->max_num_configs, pb->num_sites,
	                         vec2_hoz_old, pb->vec2_buffer);
	
	#if USE_CUDA_PROFILE_MARKS
		nvtxRangePop();
		eventAttrib.color = COLOR_GREEN;
		eventAttrib.message.ascii = "memcpy+dot z";
		nvtxRangePushEx(&eventAttrib);
	#endif
	
	
	#if USE_CUDA
		#if PRINT_FULL_DEBUG
			#pragma omp single nowait
				printf("starting gpu zpotrf\n");
		#endif
		struct timespec time_start, time_end;
		#pragma omp master
		{
			if(*pb->runtime_cuda == 0. || *pb->runtime_cpu==0.) {
				clock_gettime(CLOCK_MONOTONIC_RAW, &time_start);
			}
		}
	
		if(*pb->runtime_cuda <= *pb->runtime_cpu) {
			#pragma omp sections nowait
			{
				#pragma omp section
				{
					pb->cuda_status = cudaMemcpyAsync(pb->R_overlap_gpu, pb->phononic_overlap_matrix,
					                                  sizeof(double complex)*pb->max_num_configs*pb->max_num_configs,
					                                  cudaMemcpyHostToDevice, 0);
					pb->cusolver_status = cusolverDnZpotrf(pb->cusolver_handle, CUBLAS_FILL_MODE_UPPER,
					                                       pb->max_num_configs, pb->R_overlap_gpu,
					                                       pb->max_num_configs, pb->work_gpu, pb->Lwork_gpu, pb->info_gpu);
				}
				#pragma omp section
				{
					memset(pb->vec2_buffer, 0, pb->max_num_configs*pb->num_sites*sizeof(double complex));
				}
			}
		} else {
	#endif
			#pragma omp single nowait
			{
				memset(pb->vec2_buffer, 0, pb->max_num_configs*pb->num_sites*sizeof(double complex));
			}
	#if USE_CUDA
		}
	#endif
	
	
	int info;
	MCE1_fill_density(*pb->num_configs, pb->max_num_configs, pb->num_sites, pb->vec2_density, vec2_el_coeff);
	
	// uncoupled movement of hoz
	#pragma omp for
	for(size_t n = 0; n < pb->num_sites; n++) {
		// dx/dt
		cblas_dcopy(*pb->num_configs,
		            (double*) &vec2_hoz_old[n*pb->max_num_configs] + 1, 2,
		            (double*) &vec2_d_hoz[n*pb->max_num_configs], 2);
		cblas_dscal(*pb->num_configs, pb->omega,
		            (double*) &vec2_d_hoz[n*pb->max_num_configs], 2);
		// dp/dt
		cblas_dcopy(*pb->num_configs,
		            (double*) &vec2_hoz_old[n*pb->max_num_configs], 2,
		            (double*) &vec2_d_hoz[n*pb->max_num_configs] + 1, 2);
		cblas_dscal(*pb->num_configs, -pb->omega,
		            (double*) &vec2_d_hoz[n*pb->max_num_configs] + 1, 2);
		
		// coupling effect on hoz
		cblas_daxpy(*pb->num_configs,
		            pb->gamma, &pb->vec2_density[n*pb->max_num_configs], 1,
		            (double*) &vec2_d_hoz[n*pb->max_num_configs] + 1, 2);
	}
	
	#if USE_CUDA_PROFILE_MARKS
		nvtxRangePop();
		eventAttrib.color = COLOR_BLUE;
		eventAttrib.message.ascii = "eom_rhs";
		nvtxRangePushEx(&eventAttrib);
	#endif
	
	// electronic coefficients
	MCE1_fill_eom_rhs(vec2_d_el_coeff,
	                  *pb->num_configs, pb->max_num_configs, pb->num_sites_x, pb->num_sites_y,
	                  pb->hopping, pb->omega, pb->gamma,
	                  pb->phononic_overlap_matrix, vec2_hoz_old, vec2_el_coeff, pb->vec2_density,
	                  pb->vec2_buffer, pb->vec2_buffer2, pb->vec2_conf_overlap);
	
	#if PRINT_FULL_DEBUG
	#pragma omp single
	{
		printf("\neom_rhs=\n");
		MCE1_print_vec2(vec2_d_el_coeff, *pb->num_configs, pb->max_num_configs, pb->num_sites);
		printf("\n");
	}
	#endif
	
	#if USE_CUDA_PROFILE_MARKS
		nvtxRangePop();
		eventAttrib.color = COLOR_RED;
		eventAttrib.message.ascii = "potrf";
		nvtxRangePushEx(&eventAttrib);
	#endif
	// right hand side of the equation of motion is completed
	// solve R \dot a = vec2_eom_rhs
	
	#if USE_CUDA
		if(*pb->runtime_cuda <= *pb->runtime_cpu) {
			#if PRINT_FULL_DEBUG
				#pragma omp single
					printf("cuda synchronize\n");
			#endif
			
			#pragma omp master
			{
				cudaDeviceSynchronize();
				if(pb->cuda_status != cudaSuccess || pb->cusolver_status != CUSOLVER_STATUS_SUCCESS) {
					*pb->info = CUDA_INFO_RETURN;
				} else {
					pb->cuda_status = cudaMemcpy(pb->phononic_overlap_matrix, pb->R_overlap_gpu,
					                             sizeof(complex double)**pb->num_configs*pb->max_num_configs,
					                             cudaMemcpyDeviceToHost);
					if(pb->cuda_status != cudaSuccess) {
						*pb->info = CUDA_INFO_RETURN;
					}
				}
			}
			#if PRINT_FULL_DEBUG
				#pragma omp single nowait
					printf("cuda zpotrf results copied back\n");
			#endif
			
			#pragma omp barrier
		} else {
	#endif
	#if PRINT_FULL_DEBUG
		#pragma omp single nowait
			printf("starting cpu zpotrf\n");
	#endif
	#ifdef _OPENMP
		zpotrf_omp('U', *pb->num_configs, pb->phononic_overlap_matrix, pb->max_num_configs, pb->info);
	#else
		#pragma omp single  //not needed, because this is the no-openmp case
			*pb->info =
			LAPACKE_zpotrf(LAPACK_COL_MAJOR, 'U', *pb->num_configs,
			               (lapack_complex_double*) pb->phononic_overlap_matrix, pb->max_num_configs);
	#endif
	
	#if USE_CUDA
		}
		#pragma omp master
		{
			if(*pb->runtime_cuda == 0. || *pb->runtime_cpu==0.) {
				clock_gettime(CLOCK_MONOTONIC_RAW, &time_end);
				const double runtime = (double) (time_end.tv_sec - time_start.tv_sec) + (time_end.tv_nsec - time_start.tv_nsec) / 1e9;
				if(*pb->runtime_cuda == 0.) {
					*pb->runtime_cuda = runtime;
				} else {
					*pb->runtime_cpu = runtime;
				}
			}
		}
	#endif
	
	#if USE_CUDA_PROFILE_MARKS
		nvtxRangePop();
		eventAttrib.color = COLOR_GREEN;
		eventAttrib.message.ascii = "potrs";
		nvtxRangePushEx(&eventAttrib);
	#endif
	
	if(!*pb->info) {
		#if PRINT_FULL_DEBUG
			#pragma omp single nowait
				printf("starting zpotrs\n");
		#endif
		// split rhs to the different threads
		#ifdef _OPENMP
			div_t site_splitting = div(pb->num_sites, omp_get_num_threads());
			const int thread_num_sites = site_splitting.quot + (omp_get_thread_num() < site_splitting.rem);
			const int thread_site_offset = omp_get_thread_num()*site_splitting.quot
			                               + MIN(omp_get_thread_num(), site_splitting.rem);
		#else
			const int thread_num_sites = pb->num_sites;
			const int thread_site_offset = 0;
		#endif
		
		info =
		LAPACKE_zpotrs(LAPACK_COL_MAJOR, 'U', *pb->num_configs, thread_num_sites,
		               (const lapack_complex_double*) pb->phononic_overlap_matrix, pb->max_num_configs,
		               (lapack_complex_double*) &vec2_d_el_coeff[thread_site_offset*pb->max_num_configs], pb->max_num_configs);
		// TODO solution could be refined and error estimated by zporfs, but needs copy of phononic_overlap_matrix
		#pragma omp critical(joint_derivative_info_collect)
		{
			if(info) {
				*pb->info = info;
			}
		}
	}
	#if USE_CUDA_PROFILE_MARKS
		nvtxRangePop();
		nvtxRangePop();
	#endif
	
	#pragma omp barrier
	
	if(*pb->info) { // if solving failed -> abort
		print_lapack_solve_error(*pb->num_configs, pb->max_num_configs, pb->num_sites,
		                         pb->phononic_overlap_matrix, vec2_d_el_coeff,
		                         vec2_d_hoz, vec2_el_coeff, *pb->info);
	}
}


/*
attempts to clone some of the states in the current configuration.
any configuration is checked and just cloned, if some criteria is fullfilled TODO which one?
the maximum number of configurations given by pb->max_num_configs is never exceeded
the number of configs is changed in the given pb struct
*/
#if ALLOW_CLONING_MCEv1
void MCE1_clone_configs(struct param_and_buffer* pb,
                        double complex* restrict vec2_hoz,
                        double complex* restrict vec2_el_coeff) {
	
	// if gamma == 0, no cloning is attempted, because the systems are decoupled anyway
	if(*pb->num_configs == pb->max_num_configs || pb->gamma == 0.) {
		return;
	}
	
	
	double* const restrict norm2_contribution_config = (double*) pb->vec2_density;
	#pragma omp single nowait // barrier inside of MCE1_fill_phononic_overlap_matrix
		memset(norm2_contribution_config, 0, *pb->num_configs*sizeof(double));
	
	const size_t old_num_configs = *pb->num_configs;
	const double min_norm_square_weight = MIN(5./old_num_configs, 0.05); // TODO should depend on num_configs?
	const double min_force_diff_norm = 0.3; // max value of compared value is sqrt(2) TODO check that
	const double hoz_noise = 1e-4;
	const double min_distance = 0.3;
	const double max_overlap = exp(-min_distance*min_distance/2.);
	/*
	distance d_kj betwwen configs k and is related to R_kj by
	cabs(R_kj) = exp(-d_kj^2/2) or d_kj = sqrt(-2*ln(cabs(R_kj)))
	
	some values:
	cabs(R_kj)  | d_kj
	------------|--------
	       0.99 | 0.142
	       0.98 | 0.201
	       0.97 | 0.247
	       0.96 | 0.286
	       0.95 | 0.320
	       0.90 | 0.459
	*/
	
	MCE1_fill_overlap_matrix(pb->phononic_overlap_matrix, old_num_configs, pb->max_num_configs, pb->num_sites,
	                       vec2_hoz, pb->vec2_buffer);
	#if PRINT_CLONE && 0
		printf("\n\noverlap:\n");
		MCE1_print_overlap_matrix(pb->phononic_overlap_matrix, old_num_configs, pb->max_num_configs);
	#endif
	
	#pragma omp for
	for(size_t n = 0; n < pb->num_sites; n++) {
		// buffer = R @ a_n
		cblas_zhemv(CblasColMajor, CblasUpper, old_num_configs,
		            &cplx_one, pb->phononic_overlap_matrix, pb->max_num_configs,
		            &vec2_el_coeff[n*pb->max_num_configs], 1,
		            &cplx_zero, &pb->vec2_buffer[n*pb->max_num_configs], 1);
		vzMulByConj(old_num_configs,
		            (MKL_Complex16*) &pb->vec2_buffer[n*pb->max_num_configs],
		            (MKL_Complex16*) &vec2_el_coeff[n*pb->max_num_configs],
		            (MKL_Complex16*) &pb->vec2_buffer2[n*pb->max_num_configs]);
		
		// accumulate real part, because that sums over all configs to the squared norm of psi
		
		cblas_daxpy(old_num_configs,
		            1., (double*) &pb->vec2_buffer2[n*pb->max_num_configs], 2,
		            norm2_contribution_config, 1);
	}
	
	double norm_square_sum = 0.;
	#pragma omp for
	for(size_t k = 0; k < old_num_configs; k++) {
		for(size_t n = 0; n < pb->num_sites; n++) {
			const double nrm2 = pb->vec2_buffer2[n*pb->max_num_configs+k];
			norm2_contribution_config[k] += nrm2;
			norm_square_sum += nrm2;
		}
	}
	
	#pragma omp single
	{
		#if PRINT_CLONE
			printf("contribution to squared norm:\n");
			for(size_t k = 0; k < old_num_configs; k++) {
				printf("%lu  %g\n", k, norm2_contribution_config[k]);
			}
		#endif
		// collect sum of squared norm
		((double*) pb->vec2_buffer2)[0] = 0.;
	}
	#pragma omp critical(clone_norm2_collection)
	{
		((double*) pb->vec2_buffer2)[0] += norm_square_sum;
	}
	
	// check distance to neighbors (via element in overlap matrix) and maybe skip that configuration
	// distance d_jk corresponds to cabs(R_jk) = exp(-d_jk^2/2)
	#pragma omp for schedule(static, 1)
	for(size_t k = 0; k < old_num_configs; k++) {
		for(size_t i = k+1; i < old_num_configs; i++) {
			if(cabs(pb->phononic_overlap_matrix[k + i*pb->max_num_configs]) > max_overlap) {
				norm2_contribution_config[k] = -fabs(norm2_contribution_config[k]);
				norm2_contribution_config[i] = -fabs(norm2_contribution_config[i]);
			}
		}
	}
	// logic belongs for the loop, but here one barrier is avoided
	norm_square_sum = ((double*) pb->vec2_buffer2)[0];
	
	// argsort the norms in decending order with insertion sort, neglecting to low weights
	#pragma omp single
	{
		pb->order[0] = 0;
		for (int i = 1; i < old_num_configs; i++) {
			int key = pb->order[i] = i; // initial order is also filled on-the-fly
			int j = i - 1;
			for( ; j >= 0 && norm2_contribution_config[pb->order[j]] < norm2_contribution_config[key]; j--) {
				pb->order[j + 1] = pb->order[j];
			}
			pb->order[j + 1] = key;
		}
		#if PRINT_CLONE
			printf("checking for cloning\n%g\n", min_norm_square_weight*norm_square_sum);
			
			for(size_t k = 0; k < old_num_configs; k++) {
				printf("%i\t%g\n", pb->order[k], norm2_contribution_config[pb->order[k]]);
			}
			printf("\n");
		#endif
	}
	
	
	#pragma omp for// schedule(static, 1) // FIXME
	for(size_t k = 0 ;   k < old_num_configs; k++) {
		// these are reasons to abort, but they are imcompatible with OpenMP loops
		// so just skip those cases (after one continue, can't recover)
		if(   norm2_contribution_config[pb->order[k]] <= min_norm_square_weight*norm_square_sum
		   || *pb->num_configs >= pb->max_num_configs) {
			continue;
		}
		
		// find separted buffers. vec2_buffer2 is complex, therefore space for diag AND offdiag
		// and buffer is size max(num_threads, max_num_configs)*num_sites. same for vec2_buffer (later as force_diff)
		const int thread_ID = omp_get_thread_num();
		double* const restrict H_diag = (double*) &pb->vec2_buffer2[thread_ID*pb->num_sites];
		double* const restrict H_offdiag = &H_diag[pb->num_sites];
		
		int info;
		if(pb->num_sites_y == 1) {
			cblas_dcopy(pb->num_sites,
			            (double*) &vec2_hoz[pb->order[k]], pb->max_num_configs*2, // *2 for only real part (x/sqrt(2))
			            H_diag, 1);
			cblas_dscal(pb->num_sites, -2.*pb->gamma, H_diag, 1);
			for(size_t i = 0; i < pb->num_sites-1; i++) {
				H_offdiag[i] = -pb->hopping;
			}
			// diagonalize
			info
			= LAPACKE_dsteqr(LAPACK_COL_MAJOR, 'I', pb->num_sites,
			                 H_diag, H_offdiag, pb->BO_eigvecs, pb->num_sites);
		} else {
			memset(pb->BO_eigvecs, 0, pb->num_sites*pb->num_sites*sizeof(double));
			
			// diagonal part from coupling
			cblas_dcopy(pb->num_sites, (double*) &vec2_hoz[pb->order[k]], pb->max_num_configs*2,
			            pb->BO_eigvecs, pb->num_sites+1);
			cblas_dscal(pb->num_sites, -2.*pb->gamma, pb->BO_eigvecs, pb->num_sites+1);
			
			// offdiagonal part from hopping
			for(size_t y = 0; y < pb->num_sites_y; y++) {
				for(size_t x = 0; x < pb->num_sites_x; x++) {
					const size_t site = y*pb->num_sites_x + x;
					if(x != pb->num_sites_x-1) {
						// from right (left by symmetry)
						pb->BO_eigvecs[site + (site+1)*pb->num_sites] = -pb->hopping;
						//pb->BO_eigvecs[site+1 + site*num_sites] = -pb->hopping;
					}
					// from above (below by symmetry)
					if(y != pb->num_sites_y-1) {
						pb->BO_eigvecs[site + (site+pb->num_sites_x)*pb->num_sites] = -pb->hopping;
						//pb->BO_eigvecs[site+num_sites_x + site*num_sites] = -pb->hopping;
					}
				}
			}
			// diagonalize
			LAPACKE_dsytrd(LAPACK_COL_MAJOR, CblasUpper, pb->num_sites, pb->BO_eigvecs, pb->num_sites, H_diag, H_offdiag, pb->tau);
			LAPACKE_dorgtr(LAPACK_COL_MAJOR, CblasUpper, pb->num_sites, pb->BO_eigvecs, pb->num_sites, pb->tau);
			info
			= LAPACKE_dsteqr(LAPACK_COL_MAJOR, 'V', pb->num_sites,
			                 H_diag, H_offdiag, pb->BO_eigvecs, pb->num_sites);
		}
		
		if(info) {
			fprintf(stderr, "LAPACKE_dsteqr returned %i in cloning!\n", info);
			continue;
		}
		/*
		printf("eigenvalues:\n");
		for(size_t n = 0; n < pb->num_sites; n++) {
			printf("%g\t", H_diag[n]);
		}
		
		printf("\n");
		MCE1_print_vec2(&vec2_hoz[pb->order[k]], 1, pb->max_num_configs,pb->num_sites);
		*/
		
		const double norm_config = cblas_dznrm2(pb->num_sites, &vec2_el_coeff[pb->order[k]], pb->max_num_configs);
		
		// project to eigenbasis and normalize configuration wavefunction (and discard eigenvalues)
		// real part
		cblas_dgemv(CblasColMajor, CblasTrans, pb->num_sites, pb->num_sites,
		            1./norm_config, pb->BO_eigvecs, pb->num_sites, (double*) &vec2_el_coeff[pb->order[k]] + 0, 2*pb->max_num_configs,
		            0., (double*) &pb->vec2_buffer2[thread_ID*pb->num_sites] + 0, 2);
		// imag part
		cblas_dgemv(CblasColMajor, CblasTrans, pb->num_sites, pb->num_sites,
		            1./norm_config, pb->BO_eigvecs, pb->num_sites, (double*) &vec2_el_coeff[pb->order[k]] + 1, 2*pb->max_num_configs,
		            0., (double*) &pb->vec2_buffer2[thread_ID*pb->num_sites] + 1, 2);
		
		//printf("\nconfiguration %lu, norm=%g:\n", k, cblas_dznrm2(pb->num_sites, pb->vec2_buffer2, 1));
		//MCE1_print_vec2(pb->vec2_buffer2, pb->num_sites, pb->num_sites, 1);
		
		// search highest occupied eigenstate
		size_t index_max_eigen_weight = 0;
		double max_eigen_weight = cabs(pb->vec2_buffer2[0]);
		for(size_t n = 1; n < pb->num_sites; n++) {
			const double weight = cabs(pb->vec2_buffer2[n]);
			if(weight > max_eigen_weight) {
				max_eigen_weight = weight;
				index_max_eigen_weight = n;
			}
		}
		
		double* const force_diff = &((double*) pb->vec2_buffer)[thread_ID*pb->num_sites];
		
		double* const highest_weighted_eigenstate = &pb->BO_eigvecs[index_max_eigen_weight*pb->num_sites];
		for(size_t n = 0; n < pb->num_sites; n++) {
			const double force_eigen = highest_weighted_eigenstate[n]*highest_weighted_eigenstate[n];
			// one could also look at force acting on remaining state, but then cases like (0.99, 0.01) can be splitted although useless
			//const complex double other_z = vec2_el_coeff[pb->order[k] + n*pb->max_num_configs]/norm_config - pb->vec2_buffer2[index_max_eigen_weight]*highest_weighted_eigenstate[n];
			//const double force_other = (creal(other_z)*creal(other_z) + cimag(other_z)*cimag(other_z))/(1.-max_eigen_weight*max_eigen_weight);
			
			const complex double all_z = vec2_el_coeff[pb->order[k] + n*pb->max_num_configs]/norm_config;
			const double force_all = creal(all_z)*creal(all_z) + cimag(all_z)*cimag(all_z);
			
			force_diff[n] = force_eigen - force_all;
		
		}
		
		#if PRINT_CLONE
		#pragma omp critical(clone_print_force_diffs)
		{
		printf("\nconfiguration %i, norm2(force_diff)=%g\n", order[k], cblas_dnrm2(pb->num_sites, force_diff, 1));
		printf("q=%g, max eigenstate= %lu\n",
		       (creal(vec2_hoz[pb->order[k]+1*pb->max_num_configs])-creal(vec2_hoz[pb->order[k]+0*pb->max_num_configs]))*pb->omega/pb->gamma,
		       index_max_eigen_weight);
		
		for(size_t n = 0; n < pb->num_sites; n++) {
			printf("%g\t%g\t%g\n",
			       cabs(pb->vec2_buffer2[n])*cabs(pb->vec2_buffer2[n]),
			       //creal(pb->vec2_buffer2[n]),
			       //cimag(pb->vec2_buffer2[n]),
			       pb->BO_eigvecs[n*pb->num_sites+0],
			       pb->BO_eigvecs[n*pb->num_sites+1]);
		}
		}
		#endif
		
		// last condition. If fullfilled, cloning happnes (if no other thread took the last storage space already
		if(cblas_dnrm2(pb->num_sites, force_diff, 1) > min_force_diff_norm) {
			// synchronize cloning
			// clone_index must tbe larger than the index of the old config, which shall be cloned
			// and therefore can't be 0
			size_t clone_index = 0;
			#pragma omp critical(clone_get_index)
			{
				if(*pb->num_configs < pb->max_num_configs) {
					clone_index = (*pb->num_configs)++;
					printf("%i: clone %i to %lu\n", omp_get_thread_num(), pb->order[k], clone_index);
				} else {
					printf("%i: cloning %i to %lu refused\n", omp_get_thread_num(), pb->order[k], clone_index);
				}
			}
			if(clone_index) {
				/*
				printf("\n\ncloning config %i\nhoz:\n", pb->order[k]);
				MCE1_print_vec2(vec2_hoz, pb->num_configs, pb->max_num_configs, pb->num_sites);
				printf("\nel_coeff:\n");
				MCE1_print_vec2(vec2_el_coeff, pb->num_configs, pb->max_num_configs, pb->num_sites);
				*/
				
				// split el_coeff into most important eigenvector and rest
				const double complex eigen_weight = pb->vec2_buffer2[index_max_eigen_weight]*norm_config;
				// copy eigenstate to new config
				// real
				cblas_dcopy(pb->num_sites,
					        &pb->BO_eigvecs[index_max_eigen_weight*pb->num_sites], 1,
					        (double*) &vec2_el_coeff[clone_index] + 0, 2*pb->max_num_configs);
				cblas_dscal(pb->num_sites, creal(eigen_weight),
					        (double*) &vec2_el_coeff[clone_index] + 0, 2*pb->max_num_configs);
				// imag
				cblas_dcopy(pb->num_sites,
					        &pb->BO_eigvecs[index_max_eigen_weight*pb->num_sites], 1,
					        (double*) &vec2_el_coeff[clone_index] + 1, 2*pb->max_num_configs);
				cblas_dscal(pb->num_sites, cimag(eigen_weight),
					        (double*) &vec2_el_coeff[clone_index] + 1, 2*pb->max_num_configs);
				
				// substract new configuration from old configuration
				cblas_zaxpy(pb->num_sites,
					        &cplx_neg_one, &vec2_el_coeff[clone_index], pb->max_num_configs,
					        &vec2_el_coeff[pb->order[k]], pb->max_num_configs);
				
				// copy hoz
				cblas_zcopy(pb->num_sites,
					        &vec2_hoz[pb->order[k]], pb->max_num_configs,
					        &vec2_hoz[clone_index], pb->max_num_configs);
				
				// add some noise to the cloned hoz to make R non-singular
				// TODO get better random numbers
				for(size_t n = 0; n < pb->num_sites; n++) {
					vec2_hoz[pb->order[k] + n*pb->max_num_configs] += ((double) random()/RAND_MAX - .5)*hoz_noise + I*((double) random()/RAND_MAX - .5)*hoz_noise;
				}
				/*
				printf("\nafter cloning:\nhoz:\n");
				MCE1_print_vec2(vec2_hoz, pb->num_configs, pb->max_num_configs, pb->num_sites);
				//MCE1_print_vec2(&vec2_hoz[pb->num_configs-1], 1, pb->max_num_configs, pb->num_sites);
				printf("\nel_coeff:\n");
				MCE1_print_vec2(vec2_el_coeff, pb->num_configs, pb->max_num_configs, pb->num_sites);
				*/
			}
		}
	}
}
#endif
