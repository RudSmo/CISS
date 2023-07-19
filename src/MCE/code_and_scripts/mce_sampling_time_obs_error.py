import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import data_manager
import os

# latex style
use_tex = True
if use_tex:
	from matplotlib import rc
	rc('text', usetex=use_tex)
	rc('font', family='serif')


version = 1
hopping = 2.5
omega = 1
gamma = 2.5
num_sites_x = 2
num_sites_y = 1

num_configs = 60

num_repetitions = 20
num_disorder_scalings_ph = 5
num_disorder_scalings_el = 9

local_psi_init_select_eigenstate = False

seeds = np.random.default_rng().integers(low=0, high=(1<<32) - 1, size=num_repetitions)
#scale_phs = np.exp(np.linspace(0., np.log(1e-1), num_disorder_scalings_ph, endpoint=True))
scale_phs = np.linspace(1., 5e-1, num_disorder_scalings_ph, endpoint=True)
if local_psi_init_select_eigenstate:
	num_disorder_scalings_el = 1
	scale_els = np.zeros(num_disorder_scalings_el)
else:
	#scale_els = np.exp(np.linspace(np.log(1.), np.log(1e-1), num_disorder_scalings_el, endpoint=True))
	scale_els = np.linspace(1., 0, num_disorder_scalings_el, endpoint=True)

print(f"{scale_phs=}")
print(f"{scale_els=}")


initial_state = "bare"
psi_initial = "local"
psi_state = None

initial_site = (num_sites_x-1)//2
if psi_initial == "local":
	psi_state = (0, initial_site)
	
if initial_state == "bare":
	initial_phonon_number = 0
else:
	initial_phonon_number = (gamma/omega)**2

ref = data_manager.RefData(initial_state, psi_initial, psi_state, num_sites_x, hopping, omega, gamma)
dt = ref.time()[1]		# seems to be 0.04
print(f"{dt=}")

shown_time_indices = np.array([0, 200, 400])		# indices of shown time
shown_times = dt*shown_time_indices
t_max = max(20, max(shown_time_indices)*dt)

phonon_dressed_displaced_sites = None

load_sampled_state = False
force_new_simulation = False
printMoreProgress = False

save_data = True
show_plots = True
save_plots = True

site = (num_sites_x-1)//2

save_path = "convergence_plots/mce_sampling_time/"


obs_full = {}
obs_full[f"nel{site}"] = np.empty((num_repetitions, num_disorder_scalings_ph, num_disorder_scalings_el, len(shown_time_indices)))
obs_full[f"nph{site}"] = np.empty((num_repetitions, num_disorder_scalings_ph, num_disorder_scalings_el, len(shown_time_indices)))
obs_full["failed_sims"]= np.zeros((num_repetitions, num_disorder_scalings_ph, num_disorder_scalings_el), dtype=np.int)

max_num_configs = num_configs
print(f"{num_configs=}")
num_runs = num_repetitions*num_disorder_scalings_ph*num_disorder_scalings_el
for i in range(num_repetitions):
	seed = seeds[i]
	print("repetition", i)
	for j in range(num_disorder_scalings_ph):
		scale_ph = scale_phs[j]
		print("scale_ph", scale_ph)
		for k in range(num_disorder_scalings_el):
			scale_el = scale_els[k]
			print("scale_el", scale_el)
			subfolder_name = f"sampling_test_{scale_ph=}_{scale_el=}_eigen={local_psi_init_select_eigenstate}/{i}"
			
			sim = data_manager.Simulation_MCE(version, hopping, omega, gamma,
			                                  num_sites_x, num_sites_y,
			                                  t_max, dt,
			                                  num_configs,
			                                  initial_state, psi_initial, psi_state,
			                                  max_num_configs, printProgress=printMoreProgress,
			                                  phonon_dressed_displaced_sites=phonon_dressed_displaced_sites,
			                                  force_new_simulation=force_new_simulation,
			                                  load_sampled_state=load_sampled_state,
			                                  store_results=save_data,
			                                  extra_subfolder=subfolder_name,
			                                  local_psi_init_select_eigenstate=local_psi_init_select_eigenstate,
			                                  seed=seed, allowCompression=False,
			                                  scale_ph=scale_ph, scale_el=scale_el,
			                                  flip_initial_state_for_half_configs=False,
			                                  calcCondition = True)
			
			sim_name = sim.get_sim_name()
			
			obs_full[f"nel{site}"][i, j, k, :] = sim.nel(site)[shown_time_indices]
			obs_full[f"nph{site}"][i, j, k, :] = sim.nph(site)[shown_time_indices]
			
			if sim.nel()[-1].sum() == 0:
				obs_full["failed_sims"][i, j, k] += 1
			
			del sim
			
			print(f"{k + num_disorder_scalings_el*(j + num_disorder_scalings_ph*i)+1}/{num_runs}")

# calculate mean absolute deviations
obs_meanAbsDev = {}
obs_meanAbsDev[f"nel{site}"] = np.mean(np.abs(obs_full[f"nel{site}"]-ref.nel(site)[shown_time_indices].reshape(1,1,1,-1)), axis=0)
obs_meanAbsDev[f"nph{site}"] = np.mean(np.abs(obs_full[f"nph{site}"]-ref.nph(site)[shown_time_indices].reshape(1,1,1,-1)), axis=0)

# 
if psi_initial == "local":
	obs_meanAbsDev[f"nel{site}"][:, scale_els == 0, shown_time_indices==0] = 0


for i in range(len(shown_time_indices)):
	plt.figure()
	plt.title(r"nel $\omega_0 t = %g$"%(shown_times[i]))
	plt.pcolormesh(scale_phs, scale_els, obs_meanAbsDev[f"nel{site}"][:, :, i].T, shading="nearest", norm=matplotlib.colors.LogNorm())
	
	#plt.xscale("log")
	#plt.yscale("log")
	#ax_lim_fact = 2
	#plt.xlim((scale_phs.min()/ax_lim_fact, scale_phs.max()*ax_lim_fact))
	#plt.ylim((scale_els.min()/ax_lim_fact, scale_els.max()*ax_lim_fact))
	
	plt.xlabel(r"$\sigma_{ph}$")
	plt.ylabel(r"$\sigma_{el}$")
	
	cb = plt.colorbar()
	cb.set_label(r"$\langle |\langle n_{el, %i} \rangle_{MCE} - \langle n_{el, %i} \rangle_{ED}| \rangle_r$"%(site+1, site+1))
	
	
	plt.figure()
	plt.title(r"nph $\omega_0 t = %g$"%(shown_times[i]))
	plt.pcolormesh(scale_phs, scale_els, obs_meanAbsDev[f"nph{site}"][:, :, i].T, shading="nearest", norm=matplotlib.colors.LogNorm())
	
	plt.xlabel(r"$\sigma_{ph}$")
	plt.ylabel(r"$\sigma_{el}$")
	
	cb = plt.colorbar()
	cb.set_label(r"$\langle |\langle N_{ph, %i} \rangle_{MCE} - \langle N_{ph, %i} \rangle_{ED}| \rangle_r$"%(site+1, site+1))

plt.figure()
plt.title("number of failed simulations")
plt.pcolormesh(scale_phs, scale_els, obs_meanAbsDev["failed_sims"][:, :].T, shading="nearest", vmin=0)
plt.xlabel(r"$\sigma_{ph}$")
plt.ylabel(r"$\sigma_{el}$")

if show_plots:
	plt.show()

