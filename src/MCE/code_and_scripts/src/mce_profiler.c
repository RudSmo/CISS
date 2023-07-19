#include "sim_MCEv1_1e.h"
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[]) {
	
	if(argc != 5) {
		printf("required arguments: num_sites_x num_sites_y num_configs=max_num_configs num_times\n");
		return 1;
	}
	size_t num_sites_x = atol(argv[1]);
	size_t num_sites_y = atol(argv[2]);
	size_t num_configs = atol(argv[3]);
	
	size_t num_times = atol(argv[4]);
	double dt = 0.0005;
	
	double hopping = 1.;
	double omega = 1.;
	double gamma = 1.;
	
	FILE* f = fopen("data_export/raw_start.dat", "rb");
	if(!f) {
		printf("file not found\n");
		exit(1);
	}
	
	double complex* restrict vec2_hoz = malloc(sizeof(double complex)*num_configs*num_sites_x*num_sites_y*2);
	double complex* restrict vec2_el_coeff = vec2_hoz + num_configs*num_sites_x*num_sites_y;
	
	
	// TODO load vec2_hoz and vec2_el_coeff from file
	size_t length_input = 2*num_configs*num_sites_x*num_sites_y;
	size_t read_bytes = fread(vec2_hoz, sizeof(double complex), length_input, f);
	fclose(f);
	f = NULL;
	
	if(read_bytes != length_input) {
		printf("read failed\n");
		exit(2);
	}
	
	double complex* restrict vec2_obs_y = malloc(sizeof(double complex)*num_sites_x*num_sites_y*num_times);
	double* restrict vec2_obs_nph = malloc(sizeof(double)*num_sites_x*num_sites_y*num_times);
	double* restrict vec2_obs_nel = malloc(sizeof(double)*num_sites_x*num_sites_y*num_times);
	double* restrict vec_obs_el_ekin = malloc(sizeof(double)*num_times);
	double* restrict vec2_obs_el_ph_inter = malloc(sizeof(double)*num_sites_x*num_sites_y*num_times);
	
	
	MCE1_simulate(num_configs, num_configs,
	              num_sites_x, num_sites_y,
	              1, num_times, dt,
	              hopping, omega, gamma,
	              vec2_hoz,  // hoz = harmonic oscillator z
	              vec2_el_coeff,
	              vec2_obs_y,
	              vec2_obs_nph,
	              vec2_obs_nel,
	              vec_obs_el_ekin,
	              vec2_obs_el_ph_inter,
	              1);
	
	free(vec2_hoz);
	free(vec2_obs_y);
	free(vec2_obs_nph);
	free(vec2_obs_nel);
	free(vec_obs_el_ekin);
	free(vec2_obs_el_ph_inter);
	
	
	return 0;
}
