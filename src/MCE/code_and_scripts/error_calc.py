import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import simps

import data_manager

def error(sim_time, sim_obs, ref_time, ref_observable):
	interp_ref = interp1d(ref_time, ref_observable, axis=0, kind="cubic", assume_sorted=True, copy=False)
	diff = sim_obs - interp_ref(sim_time)
	return np.sqrt(np.mean(simps(diff**2, dx=dt, axis=0)))		# mean over sites


initial_state = "BCDW"
psi_initial="local"
hopping = 1.
omega, gamma = .1, 0.4
t_max = 6.5
symmetrize_start = False

output_file = open("error_ehrenfest/%s_omega=%g_gamma=%g.dat"%(initial_state, omega, gamma), "a")


dt = 0.01
num_sites = 2

print("#sites\tdt\tnum_samples\terror_nel\terror_nph\terror_x\terror_p\terror_ekin", file=output_file)

for num_sites in [2, 5, 7, 9, 11, 13]:
	print("sites=", num_sites)
	ref = data_manager.RefData(initial_state, psi_initial, num_sites, omega, gamma)
	ref.set_cut_time(t_max)

	for num_samples in [400, 1000, 4000]:
		print("samples=", num_samples)
		
		for i in range(10):
			sim = data_manager.Simulation(omega, gamma, t_max, dt,
						                  num_samples, num_sites=num_sites,
						                  symmetric_y_init=symmetrize_start,
						                  num_electrons=num_sites//2,
						                  printProgress=False,
						                  initial_state=initial_state,
						                  psi_initial=psi_initial,
						                  calcYSquared=True)
			
			error_x = error(sim.time(), sim.x(), ref.time(), ref.x())
			error_p = error(sim.time(), sim.p(), ref.time(), ref.p())
			error_nel = error(sim.time(), sim.nel(), ref.time(), ref.nel())
			error_nph = error(sim.time(), sim.nph(), ref.time(), ref.nph())
			if num_sites == 2:
				error_ekin = np.nan
			else:
				error_ekin = error(sim.time(), sim.energy_kinetic(), ref.time(), ref.energy_kinetic())
			
			print(num_sites, num_samples, dt, error_nel, error_nph, error_x, error_p, error_ekin, file=output_file)

output_file.close()
