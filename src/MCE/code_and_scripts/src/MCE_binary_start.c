#include "sim_MCEv1_1e.h"
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>


#define NAME_BUFFER_SIZE 256
#define OUTPUT_DIR "MCE_binary_results/"
#define INPUT_DIR "MCE_binary_start/"

void save_observable(char* io_file_name, const int io_file_name_size,
                     const char* output_path, const char* observable_name,
                     void* observable, size_t itemsize, size_t num_items) {
	
	int char_written
	= snprintf(io_file_name, io_file_name_size, "%s%s.bin", output_path, observable_name);
	
	if(char_written > 0 && char_written < io_file_name_size) {
		FILE* f = fopen(io_file_name, "wb");
		fwrite(observable, itemsize, num_items, f);
		fclose(f);
	} else {
		fprintf(stderr, "filename buffer too small!\n");
	}
}

// TODO allow MCEv2
int main(int argc, char* argv[]) {
	
	if(argc != 12) {
		fprintf(stderr, "required arguments: ph_init el_init num_sites_x num_sites_y num_configs max_num_configs dt num_times hopping omega gamma\n");
		return 1;
	}
	
	size_t num_sites_x     = atol(argv[ 3]);
	size_t num_sites_y     = atol(argv[ 4]);
	size_t num_configs     = atol(argv[ 5]);
	size_t max_num_configs = atol(argv[ 6]);
	
	if(max_num_configs < num_configs) {
		fprintf(stderr, "max_num_configs < num_configs is not allowed, falling back to max_num_configs = num_configs\n");
		max_num_configs = num_configs;
	}
	
	double dt              = atof(argv[ 7]);
	size_t num_times       = atol(argv[ 8]);
	
	double hopping         = atof(argv[ 9]);
	double omega           = atof(argv[10]);
	double gamma           = atof(argv[11]);
	
	char name[NAME_BUFFER_SIZE];
	char io_file_name[NAME_BUFFER_SIZE] = OUTPUT_DIR;
	
	
	// set up all filenames and folders
	if(mkdir(OUTPUT_DIR, S_IRWXU) != 0 && errno != EEXIST) {
		fprintf(stderr, "creating the parent output directory failed!\n");
		return 1;
	}
	
	int char_written // TODO update name
	= snprintf(name, NAME_BUFFER_SIZE, "MCE_%s_%s_Lx=%lu_Ly=%lu_hopping=%g_omega=%g_gamma=%g_numconf=%lu",
	          argv[1], argv[2], num_sites_x, num_sites_y, hopping, omega, gamma, num_configs);
	if(char_written <= 0 || char_written >= NAME_BUFFER_SIZE-1) {
		fprintf(stderr, "filename buffer couldn't be build!\n");
		return 1;
	}
	
	char_written
	= snprintf(io_file_name, NAME_BUFFER_SIZE, INPUT_DIR"%s.bin", name);
	if(char_written <= 0 || char_written >= NAME_BUFFER_SIZE) {
		fprintf(stderr, "filename buffer couldn't be build!\n");
		return 1;
	}
	
	
	char output_path[NAME_BUFFER_SIZE];
	// create output directory
	char_written
	= snprintf(output_path, NAME_BUFFER_SIZE, OUTPUT_DIR "%s-%lu/", name, max_num_configs);
	if(char_written <= 0 || char_written >= NAME_BUFFER_SIZE) {
		fprintf(stderr, "output_path couldn't be build!\n");
		return 1;
	}
	if(mkdir(output_path, S_IRWXU) != 0 && errno != EEXIST) {
		fprintf(stderr, "creating the output directory failed!\n");
		return 1;
	}
	
	printf("output_path=%s\n\n", output_path);
	
	// read the input data
	FILE* f = fopen(io_file_name, "rb");
	if(!f) {
		fprintf(stderr, "file \"%s\" not found\n", io_file_name);
		return 1;
	}
	
	double complex* restrict vec2_hoz = malloc(sizeof(double complex)*max_num_configs*num_sites_x*num_sites_y);
	double complex* restrict vec2_el_coeff = malloc(sizeof(double complex)*max_num_configs*num_sites_x*num_sites_y);
	
	size_t read_bytes = 0;
	// load vec2_hoz from file
	for(size_t n = 0; n < num_sites_x*num_sites_y; n++) {
		read_bytes += fread(&vec2_hoz[n*max_num_configs], sizeof(double complex), num_configs, f);
	}
	// load vec2_el_coeff from file (directly follwing vec2_hoz values)
	for(size_t n = 0; n < num_sites_x*num_sites_y; n++) {
		read_bytes += fread(&vec2_el_coeff[n*max_num_configs], sizeof(double complex), num_configs, f);
	}
	
	fclose(f);
	f = NULL;
	
	if(read_bytes != 2*num_configs*num_sites_x*num_sites_y) {
		fprintf(stderr, "reading input file failed\n");
		free(vec2_hoz);
		free(vec2_el_coeff);
		return 2;
	}
	
	double complex* restrict vec2_obs_y = malloc(sizeof(double complex)*num_sites_x*num_sites_y*num_times);
	double* restrict vec2_obs_nph = malloc(sizeof(double)*num_sites_x*num_sites_y*num_times);
	double* restrict vec2_obs_nel = malloc(sizeof(double)*num_sites_x*num_sites_y*num_times);
	double* restrict vec_obs_el_ekin = malloc(sizeof(double)*num_times);
	double* restrict vec2_obs_el_ph_inter = malloc(sizeof(double)*num_sites_x*num_sites_y*num_times);
	
	
	
	
	MCE1_simulate(num_configs, max_num_configs,
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
	
	
	// TODO write all observables
	save_observable(io_file_name, NAME_BUFFER_SIZE, output_path, "obs_y",
	                vec2_obs_y, sizeof(double complex), num_sites_x*num_sites_y*num_times);
	save_observable(io_file_name, NAME_BUFFER_SIZE, output_path, "obs_nph",
	                vec2_obs_nph, sizeof(double), num_sites_x*num_sites_y*num_times);
	save_observable(io_file_name, NAME_BUFFER_SIZE, output_path, "obs_nel",
	                vec2_obs_nel, sizeof(double), num_sites_x*num_sites_y*num_times);
	save_observable(io_file_name, NAME_BUFFER_SIZE, output_path, "obs_el_ekin",
	                vec_obs_el_ekin, sizeof(double), num_times);
	save_observable(io_file_name, NAME_BUFFER_SIZE, output_path, "obs_el_ph_inter",
	                vec2_obs_el_ph_inter, sizeof(double), num_sites_x*num_sites_y*num_times);
	
	// overwrite vec_obs_el_ekin by time and save that too
	for(size_t t = 0; t < num_times; t++) {
		vec_obs_el_ekin[t] = dt*t;
	}
	save_observable(io_file_name, NAME_BUFFER_SIZE, output_path, "time",
	                vec_obs_el_ekin, sizeof(double), num_times);
	
	
	abort_simulation:
	
	
	free(vec2_obs_y);
	free(vec2_obs_nph);
	free(vec2_obs_nel);
	free(vec_obs_el_ekin);
	free(vec2_obs_el_ph_inter);
	
	free(vec2_hoz);
	free(vec2_el_coeff);
	
	return 0;
}
