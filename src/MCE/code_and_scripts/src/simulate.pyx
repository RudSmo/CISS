#!python
#cython: language_level=3

#cimport numpy as np
import numpy as np
import ctypes
import signal


cdef extern from "sim_MTE.h":
	void MTE_simulate(size_t num_samples,
	                  size_t num_sites_x,
	                  size_t num_sites_y,
	                  size_t num_electrons,
	                  size_t num_times,
	                  double dt,
	                  double tol_energy_error_per_site,
	                  const double* y_initial,
	                  const complex* psi_initial,
	                  double hopping,
	                  double omega,
	                  double gamma,
	                  double* y_mean_a_t,
	                  double* density_a_t,
	                  double* y_squared_a_t,
	                  int printProgress,
	                  double* y_full_a_t,
	                  double* el_kinetic_energy_a_t,
	                  double* energy_interaction_a_t);

cdef extern from "sim_MCE.h":
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
	                 complex* vec2_hoz,
	                 complex* vec2_el_coeff,
	                 complex* vec_coupling_coeff,
	                 complex* vec2_obs_y,
	                 double* vec2_obs_nph,
	                 double* vec2_obs_nel,
	                 double* vec_obs_el_ekin,
	                 double* vec2_obs_el_ph_inter,
	                 double* vec2_obs_BO_nel,
	                 int obs_BO_nel_save_dist,
	                 size_t* fullStorage_num_configs,
	                 complex* fullStorage_hoz,
	                 complex* fullStorage_el_coeff,
	                 int printProgress);


def simulate_system_MTE(y_initial,
                        psi_initial,
                        double hopping, double omega, double gamma,
                        num_times, time_step, tol_energy_error_per_site=1e-6,
                        returnNumPhonons=False,
                        returnYSquared=False,
                        printProgress=True,
                        fullSampleReturn=False,
                        returnElEkin=True,
                        returnInteractionEnergy=True):
	cdef size_t num_samples = y_initial.shape[0]
	cdef size_t num_sites_y = y_initial.shape[1]
	cdef size_t num_sites_x = y_initial.shape[2]
	cdef size_t num_electrons = psi_initial.shape[1]
	
	y_all_times = np.zeros((num_times, num_sites_y, num_sites_x, 2), dtype=np.float)
	density_all_times = np.zeros((num_times, num_sites_y, num_sites_x), dtype=np.float)
	
	cdef int printProg = printProgress
	cdef double* y_squared_all_times_p = <double*> 0u
	y_squared_all_times = None
	if returnNumPhonons or returnYSquared:
		y_squared_all_times = np.zeros((num_times, num_sites_y, num_sites_x, 2), dtype=np.float)
		y_squared_all_times_p = <double*> <size_t> y_squared_all_times.__array_interface__['data'][0]
	
	cdef double* full_y_data_pointer = <double*> 0u
	full_y_data = None
	if fullSampleReturn:
		full_y_data = np.empty((num_times, num_samples, num_sites_y, num_sites_x, 2))
		full_y_data_pointer = <double*> <size_t> full_y_data.__array_interface__['data'][0]
	
	cdef double* el_kinetic_energy_pointer = <double*> 0u
	el_kinetic_energy = None
	if returnElEkin:
		el_kinetic_energy = np.empty((num_times, num_electrons))
		el_kinetic_energy_pointer = <double*> <size_t> el_kinetic_energy.__array_interface__['data'][0]
	
	cdef double* energy_interaction_pointer = <double*> 0u
	energy_interaction = None
	if returnInteractionEnergy:
		energy_interaction = np.empty((num_times, num_sites_y, num_sites_x))
		energy_interaction_pointer = <double*> <size_t> energy_interaction.__array_interface__['data'][0]
	
	
	# call the actual simulation code
	MTE_simulate(num_samples, num_sites_x, num_sites_y, num_electrons,
	             num_times, time_step, tol_energy_error_per_site,
	             <double*> <size_t> y_initial.__array_interface__['data'][0],
	             <complex*> <size_t> psi_initial.__array_interface__['data'][0],
	             hopping, omega,
	             gamma,
	             <double*> <size_t> y_all_times.__array_interface__['data'][0],
	             <double*> <size_t> density_all_times.__array_interface__['data'][0],
	             y_squared_all_times_p, printProg,
	             full_y_data_pointer,
	             el_kinetic_energy_pointer,
	             energy_interaction_pointer);
	
	returnValue = {}
	returnValue["y_mean"] = y_all_times
	returnValue["density"] = density_all_times
	returnValue["y_squared"] = y_squared_all_times
	returnValue["y_full"] = full_y_data
	returnValue["el_kinetic_energy"] = el_kinetic_energy
	returnValue["energy_interaction"] = energy_interaction
	
	return returnValue


# fullStorage always returns, but if it wasn't request, that value is None (instead of np.ndarray)
# same for BO_electron_density
def simulate_system_MCE(hoz_initial,
                        el_coeff_initial,
                        coupling_coeff_initial,
                        double hopping, double omega, double gamma,
                        size_t num_times, double dt,
                        size_t num_configs,
                        int allow_config_delete,
                        int allow_config_clone,
                        printProgress=True,
                        fullStorage=False,
                        calc_BO_electron_density=False,
                        BO_nel_save_dist=1):
	
	use_MCEv1 = coupling_coeff_initial is None
	
	if use_MCEv1:
		assert(hoz_initial.shape[0] == el_coeff_initial.shape[1])
		assert(hoz_initial.shape[1] == el_coeff_initial.shape[2])
		assert(hoz_initial.shape[2] == el_coeff_initial.shape[3])
	else:
		assert(hoz_initial.shape[0] == el_coeff_initial.shape[0])
		assert(hoz_initial.shape[1] == el_coeff_initial.shape[2	])
		assert(hoz_initial.shape[2] == el_coeff_initial.shape[3])
	
	coupling_coeff_initial_p = <complex*> <size_t> 0u
	if not use_MCEv1:
		coupling_coeff_initial_p = <complex*> <size_t> coupling_coeff_initial.__array_interface__['data'][0]
	
	# TODO assert, that for every configuration at least one coefficient is not 0
	cdef size_t num_sites_x = hoz_initial.shape[1] if use_MCEv1 else hoz_initial.shape[2]
	cdef size_t num_sites_y = hoz_initial.shape[0] if use_MCEv1 else hoz_initial.shape[1]
	cdef size_t max_num_configs = hoz_initial.shape[2] if use_MCEv1 else hoz_initial.shape[0]
	cdef size_t num_electrons = el_coeff_initial.shape[0] if use_MCEv1 else el_coeff_initial.shape[1]
	
	obs_y   = np.zeros((num_times, num_sites_y, num_sites_x, 2), dtype=np.float)
	obs_nph = np.zeros((num_times, num_sites_y, num_sites_x), dtype=np.float)
	obs_nel = np.zeros((num_times, num_sites_y, num_sites_x), dtype=np.float)
	obs_el_ekin = np.zeros((num_times, num_electrons), dtype=np.float)
	obs_el_ph_inter = np.zeros((num_times, num_sites_y, num_sites_x), dtype=np.float)
	
	
	fullStorage_num_configs = None
	fullStorage_num_configs_p = <size_t*> <size_t> 0u
	
	if fullStorage or allow_config_delete or allow_config_clone:
		fullStorage_num_configs = np.empty(num_times, dtype=np.int64)
		fullStorage_num_configs_p = <size_t*> <size_t> fullStorage_num_configs.__array_interface__['data'][0]
	
	fullStorage_el_coeff = None
	fullStorage_el_coeff_p = <complex*> <size_t> 0u
	
	fullStorage_hoz = None
	fullStorage_hoz_p = <complex*> <size_t> 0u
	
	if fullStorage:
		fullStorage_el_coeff   = np.empty((num_times, num_sites_y, num_sites_x, max_num_configs), dtype=np.complex)
		fullStorage_el_coeff_p = <complex*> <size_t> fullStorage_el_coeff.__array_interface__['data'][0]
		
		fullStorage_hoz   = np.empty((num_times, num_sites_y, num_sites_x, max_num_configs), dtype=np.complex)
		fullStorage_hoz_p = <complex*> <size_t> fullStorage_hoz.__array_interface__['data'][0]
	
	
	BO_nel = None
	BO_nel_p = <double*> <size_t> 0u
	if calc_BO_electron_density:
		if type(BO_nel_save_dist) != int or BO_nel_save_dist <= 0:
			raise TypeError("BO_nel_save_dist must be a positive integer")
		BO_nel = np.zeros((1+(num_times-1)//BO_nel_save_dist, num_sites_y*num_sites_x), dtype=np.float)
		BO_nel_p = <double*> <size_t> BO_nel.__array_interface__['data'][0]
		BO_nel_num_times = BO_nel.shape[0]
	
	original_sigint_handler = signal.getsignal(signal.SIGINT)
	
	failed = \
	MCE_simulate(num_configs, max_num_configs,
	             allow_config_delete, allow_config_clone,
	             num_sites_x, num_sites_y, num_electrons, num_times,
	             dt, hopping, omega, gamma,
	             <complex*> <size_t> hoz_initial.__array_interface__['data'][0],
	             <complex*> <size_t> el_coeff_initial.__array_interface__['data'][0],
	             <complex*> <size_t> coupling_coeff_initial_p,
	             <complex*> <size_t> obs_y.__array_interface__['data'][0],
	             <double*> <size_t> obs_nph.__array_interface__['data'][0],
	             <double*> <size_t> obs_nel.__array_interface__['data'][0],
	             <double*> <size_t> obs_el_ekin.__array_interface__['data'][0],
	             <double*> <size_t> obs_el_ph_inter.__array_interface__['data'][0],
	             BO_nel_p, BO_nel_save_dist,
	             fullStorage_num_configs_p,
	             fullStorage_hoz_p,
	             fullStorage_el_coeff_p,
	             <int> printProgress)
	
	signal.signal(signal.SIGINT, original_sigint_handler)
	
	if failed:
		print("MCE simulation failed with return value %i!"%failed)
	
	return obs_y, obs_nph, obs_nel, obs_el_ekin, obs_el_ph_inter,\
	       fullStorage_num_configs, fullStorage_el_coeff, fullStorage_hoz, BO_nel
