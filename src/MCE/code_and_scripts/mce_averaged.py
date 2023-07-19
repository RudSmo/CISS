import numpy as np
import matplotlib.pyplot as plt
import os
import time as systime

from src.ehrenfest_init import initial_state_MCE
from simulate import simulate_system_MCE

import data_manager

hopping = 1.
omega = 2
gamma = 4
num_sites_x = 2
num_sites_y = 1
num_electrons = 1

random_hoz = True	# very evil, if disabled
symmetrize_start = False

useDisorderRegularization = False  # False is better -> use eigenbasis to invert
regularize_strength = 0e-6
cutoff_eigval_overlap = 0e-3
qss = True
num_configs = 40
max_num_configs = num_configs
num_reps = 20
realization_alpha=0.2

regularize_diag = False


save_dt = 0.05
#dt = min(0.05, np.inf if omega*gamma == 0. else np.sqrt(1e-8/omega/gamma**2))
#print(dt)
#quit()
dt = 0.1  # 0.0001
save_step = 1 #int(np.round(save_dt/dt))
t_max = 50. #6*2*np.pi

initial_state = "BCDW"
psi_initial = "local"

save_plots = False
figsize = (8, 4)
title = initial_state + " " + psi_initial + ": $\\omega=%g$, $\\gamma=%g,$ conf %i$\\rightarrow$%i"%(omega, gamma, num_configs, max_num_configs)

try:
	ref = data_manager.RefData(initial_state, psi_initial, num_sites_x, omega, gamma)
	ref.set_cut_time(t_max)
	t_max = ref.time()[-1]
except Exception as e:
	print("no reference Data found:")
	print(repr(e))
	ref = None

energy_weight_power = 1
min_energy_diff = 1e-4

num_times = int(np.round(t_max/dt)) + 1  # +1 for t=0
time = np.linspace(0, t_max, num_times)
weight_sum = np.zeros((num_times, 1, 1))

weights_full = np.zeros((num_reps, num_times))

obs_nel_m  = np.zeros((num_times, num_sites_y, num_sites_x))
obs_nel_wm = np.zeros((num_times, num_sites_y, num_sites_x))

obs_nph_m  = np.zeros((num_times, num_sites_y, num_sites_x))
obs_nph_wm = np.zeros((num_times, num_sites_y, num_sites_x))

obs_nel_full = np.zeros((num_reps, num_times, num_sites_y, num_sites_x))
obs_nph_full = np.zeros((num_reps, num_times, num_sites_y, num_sites_x))

fig_nel, ax_nel = plt.subplots(1, 1)
fig_nel.suptitle("nel "+title)
ax_nel.set_xlabel("time $tt_0$")
ax_nel.set_ylabel("nel")


fig_nph, ax_nph = plt.subplots(1, 1)
fig_nph.suptitle("nph "+title)
ax_nph.set_xlabel("time $tt_0$")
ax_nph.set_ylabel("nph")

fig_etot, ax_etot = plt.subplots(1, 1)

time_start = systime.time()

for rep in range(num_reps):
	# initial conditions
	hoz, el_coeff = initial_state_MCE(num_sites_x, num_sites_y, num_configs, max_num_configs, num_electrons, 
	                                  symmetrize_start, initial_state, psi_initial,
	                                  hopping, omega, gamma,
	                                  regularize_diag, regularize_strength,
	                                  plotOverlapMatrix=False, use_QSS=qss, random_hoz=random_hoz,
	                                  cutoff_eigval_overlap=cutoff_eigval_overlap,
	                                  useDisorderRegularization=useDisorderRegularization)
	
	# simulation
	obs_y, obs_nph, obs_nel, obs_el_ekin, obs_el_ph_inter = simulate_system_MCE(hoz, el_coeff, hopping, omega, gamma, num_times, dt, num_configs, printProgress=True)
	
	e_tot =   omega * obs_nph.sum(axis=(1, 2))\
	        + gamma * obs_el_ph_inter.sum(axis=(1, 2))\
	        + hopping * obs_el_ekin.sum(axis=1)
	#plt.title(r"$\Delta E_{tot}$ " + title)
	e_tot -= e_tot[0]
	
	e_tot_maxed = np.abs(e_tot)
	max_e_tot = 0.
	for i in range(1, e_tot_maxed.size):
		if e_tot_maxed[i] > max_e_tot:
			max_e_tot = e_tot_maxed[i]
		else:
			e_tot_maxed[i] = max_e_tot
	
	weight = (np.where(e_tot_maxed > min_energy_diff, e_tot_maxed, min_energy_diff)**(-energy_weight_power)).reshape(num_times, 1, 1)
	weights_full[rep] = weight.flatten()
	
	weight_sum += weight
	
	ax_nel.set_prop_cycle(None)
	ax_nel.plot(time, obs_nel[:, 0, :], alpha=realization_alpha)
	
	obs_nel_m  += obs_nel
	obs_nel_wm += weight*obs_nel
	
	ax_nph.set_prop_cycle(None)
	ax_nph.plot(time, obs_nph[:, 0, :], alpha=realization_alpha)
	
	obs_nph_m  += obs_nph
	obs_nph_wm += weight*obs_nph
	
	
	ax_etot.plot(time, e_tot)

print("runtime:", systime.time()-time_start)

obs_nel_m  /= num_reps
obs_nel_wm /= weight_sum
obs_nph_m  /= num_reps
obs_nph_wm /= weight_sum


plt.figure()
plt.title("weights")
weights_full /= np.sum(weights_full, axis=0, keepdims=True)
plt.plot(time, weights_full.T)

ax_nel.set_prop_cycle(None)
ax_nel.plot(time, obs_nel_m[:, 0, :])[0].set_label("mean")
#ax_nel.set_prop_cycle(None)
ax_nel.plot(time, obs_nel_wm[:, 0, :], '--')[0].set_label("weighted mean")
if ref is not None:
	ax_nel.plot(ref.time(), ref.nel(), '-.')[0].set_label("ref")
ax_nel.legend()

ax_nph.set_prop_cycle(None)
ax_nph.plot(time, obs_nph_m[:, 0, :])[0].set_label("mean")
#ax_nph.set_prop_cycle(None)
ax_nph.plot(time, obs_nph_wm[:, 0, :], '--')[0].set_label("weighted mean")
if ref is not None:
	ax_nph.plot(ref.time(), ref.nph(), '-.')[0].set_label("ref")
ax_nph.legend()


plt.figure()
plt.title("nel error")
plt.plot(time, obs_nel_m[:, 0, :]-ref.nel())
plt.gca().set_prop_cycle(None)
plt.plot(time, obs_nel_wm[:, 0, :]-ref.nel(), '--')


plt.figure()
plt.title("nph error")
plt.plot(time, obs_nph_m[:, 0, :]-ref.nph())
plt.gca().set_prop_cycle(None)
plt.plot(time, obs_nph_wm[:, 0, :]-ref.nph(), '--')

plt.show()

