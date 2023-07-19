import numpy as np
import matplotlib.pyplot as plt
import data_manager
import os

# latex style
use_tex = False
if use_tex:
	from matplotlib import rc
	rc('text', usetex=use_tex)
	rc('font', family='serif')

hopping = 0.5
omega = 1
gamma = 2.

initial_state = "dressed"
psi_initial = "local"
psi_state = (0, 0)

num_configs = 200
num_configs1           = num_configs
num_configs2_local     = num_configs
num_configs2_localFlip = num_configs
num_configs2_eigen     = num_configs
num_configs2_eigenFlip = num_configs

sim_name = "MCE2localFlip"  # MCE1local MCE2local MCE2localFlip MCE2eigen MCE2eigenFlip
extra_subfolder = "0"
save_simulation_data = True        # whether qpq file sshall be stored
overwrite_simulation_data = False  # whether regular simulation data shall be overwritten.
                                   # auto switched to True if save_simulation_data and no simualtion found
save_config_placements = True
load_sampled_state = True
printProgress = True
save_plots = True
num_contour_lines = 28

shown_times = np.arange(0., 301., 10.)

sim_t_max = shown_times.max()
dt = shown_times[1]
save_path = "dimer_tunnel_plot/"

if sim_name == "MCE1local":
	used_num_configs = num_configs1
if sim_name == "MCE2local":
	used_num_configs = num_configs2_local
if sim_name == "MCE2localFlip":
	used_num_configs = num_configs2_localFlip
if sim_name == "MCE2eigen":
	used_num_configs = num_configs2_eigen
if sim_name == "MCE2eigenFlip":
	used_num_configs = num_configs2_eigenFlip

# create a fake simulation with no simulation time or storage to construct name and search for data
fakeSim = data_manager.Simulation_MCE(int(sim_name[3]), hopping, omega, gamma,
                                      2, 1,
                                      0, dt,
                                      used_num_configs,
                                      initial_state, psi_initial, psi_state,
                                      used_num_configs,
                                      allow_config_delete=False,
                                      allow_config_clone=False,
                                      printProgress=False,
                                      load_sampled_state=load_sampled_state,
                                      store_results=False,
                                      local_psi_init_select_eigenstate="eigen" in sim_name,
                                      flip_initial_state_for_half_configs="Flip" in sim_name,
                                      extra_subfolder=extra_subfolder)

data_path = fakeSim.get_storage_path(extra_subfolder=extra_subfolder)
save_path += fakeSim.get_sim_name().replace(" ", "")
print(f"{data_path=}")
print(f"{save_path=}")
if save_plots:
	os.system("mkdir -p " + save_path.replace("(", "\(").replace(")", "\)"))

overwrite_simulation_data = overwrite_simulation_data or (save_simulation_data and not os.path.isdir(data_path))

qpq_filename_template  = data_path + "fullQPQ_wt=%g.txt"

try:
	# TODO load data from files
	q  = np.empty((len(shown_times), num_configs))
	pq = np.empty((len(shown_times), num_configs))
	for t in range(len(shown_times)):
		data = np.loadtxt(qpq_filename_template%shown_times[t])
		q[t] = data[:, 0]
		pq[t] = data[:, 1]
	time = shown_times
except OSError:
	sim = data_manager.Simulation_MCE(int(sim_name[3]), hopping, omega, gamma,
	                                  2, 1,
	                                  sim_t_max, dt,
	                                  used_num_configs,
	                                  initial_state, psi_initial, psi_state,
	                                  used_num_configs,
	                                  allow_config_delete=False,
	                                  allow_config_clone=False,
	                                  printProgress=printProgress,
	                                  load_sampled_state=load_sampled_state,
	                                  store_results=overwrite_simulation_data,
	                                  local_psi_init_select_eigenstate="eigen" in sim_name,
	                                  flip_initial_state_for_half_configs="Flip" in sim_name,
	                                  fullStorage=True,
	                                  extra_subfolder=extra_subfolder)
	
	# not scaled with sqrt(2) because that cancels: q = (x1-x2)/sqrt(2) = (hoz1*sqrt(2) - hoz2*sqrt(2))/sqrt(2)
	q  = sim.fullStorage_hoz[:, 0, 0, :].real - sim.fullStorage_hoz[:, 0, 1, :].real
	pq = sim.fullStorage_hoz[:, 0, 0, :].imag - sim.fullStorage_hoz[:, 0, 1, :].imag
	time = sim.time()
	
	# save the data which will be plotted
	if save_simulation_data:
		for show_time in shown_times:
			show_time_ind = np.argmax(time >= show_time)
			data = np.stack((q[show_time_ind], pq[show_time_ind])).T
			filename = qpq_filename_template%time[show_time_ind]
			np.savetxt(filename, data, header="q, pq for all configurations")
	




if np.all(q[1] == 0):
	quit()

def BO_surfaces(q):
	
	half_q2 = .5*q**2
	e_diff = np.sqrt(q**2*gamma**2 + hopping**2)
	
	return np.stack((half_q2 - e_diff, half_q2+e_diff))


q_range  = max(q.max(), -q.min())
pq_range = max(pq.max(), -pq.min())
q_range = 10
pq_range = 10
print(f"{q_range=}")
print(f"{pq_range=}")

label_q  = r"$q/\sqrt{\frac{\hbar}{m\omega_0}}$"
label_pq = r"$p_q/\sqrt{\hbar m\omega_0}$"

q_grid = np.linspace(-q_range, q_range, num=100)
energy_q = BO_surfaces(q_grid)[0]

pq_grid = np.linspace(-pq_range, pq_range, num=100)
energy_pq = .5*pq_grid**2

energy_full = energy_q.reshape(q_grid.size, 1) + energy_pq.reshape(1, pq_grid.size)


for show_time in shown_times:
	show_time_ind = np.argmax(time >= show_time)
	
	fig, ax = plt.subplots(1, 1, figsize=(6/2, 3))
	#fig.suptitle(r"$\omega_0t = %g$"%show_time)
	cm = ax.contour(q_grid, pq_grid, energy_full.T, num_contour_lines)
	#cb = fig.colorbar(cm)
	
	if "Flip" in sim_name:
		ax.plot(q[show_time_ind][:num_configs//2], pq[show_time_ind][:num_configs//2], '.')
		ax.plot(q[show_time_ind][num_configs//2:], pq[show_time_ind][num_configs//2:], '+')
	else:
		ax.plot(q[show_time_ind], pq[show_time_ind], '.')
	
	ax.set_xlabel(label_q)
	ax.set_xlim((-q_range, q_range))
	ax.set_ylabel(label_pq)
	ax.set_ylim((-pq_range, pq_range))
	#cb.set_label(r"$E_0/\hbar\omega_0$")
	fig.tight_layout(pad=0.7)
	if save_plots:
		fig.savefig(save_path + "/wt=%g.pdf"%show_time)

if not save_plots:
	plt.show()
