import numpy as np
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

num_configs_list = [20, 40, 60] # 100, 200]
num_configs_label_list = [str(nc) for nc in num_configs_list]

num_repetitions = 1
num_disorder_scalings = 20

seeds = np.random.default_rng().integers(low=0, high=(1<<32) - 1, size=num_repetitions)
scale_phs = np.exp(np.linspace(0., np.log(1e-3), num_disorder_scalings, endpoint=False))

local_psi_init_select_eigenstate = False

dt = 0.01
t_max = 0
save_step = 1

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

phonon_dressed_displaced_sites = None

line_fmt = '.-'

load_sampled_state = False
force_new_simulation = False

save_data = True
show_plots = True
save_plots = True

save_path = "convergence_plots/mce_sampling/"


norms        = np.empty((len(num_configs_list), num_repetitions, num_disorder_scalings))
nph_init     = np.empty((len(num_configs_list), num_repetitions, num_disorder_scalings))
nel_init     = np.empty((len(num_configs_list), num_repetitions, num_disorder_scalings))
nel_offinit  = np.empty((len(num_configs_list), num_repetitions, num_disorder_scalings))
condition_tot = np.empty((len(num_configs_list), num_repetitions, num_disorder_scalings))
condition_ph = np.empty((len(num_configs_list), num_repetitions, num_disorder_scalings))

for k in range(len(num_configs_list)):
	num_configs = num_configs_list[k]
	max_num_configs = num_configs
	print(f"{num_configs=}")
	for i in range(num_repetitions):
		seed = seeds[i]
		for j in range(num_disorder_scalings):
			scale_ph = scale_phs[j]
			if local_psi_init_select_eigenstate:
				scale_el = 0
			else:
				scale_el = scale_ph/2 #num_sites_x
			subfolder_name = f"sampling_test_{scale_ph=}_{scale_el=}_eigen={local_psi_init_select_eigenstate}/{i}"
			
			sim = data_manager.Simulation_MCE(1, hopping, omega, gamma,
			                                  num_sites_x, num_sites_y,
			                                  t_max, dt,
			                                  num_configs,
			                                  initial_state, psi_initial, psi_state,
			                                  max_num_configs, printProgress=False,
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
			
			norms        [k, i, j] = np.abs(1 - sim.finalInitParams["total_norm"])
			nel_init     [k, i, j] = np.abs(1 - sim.nel(initial_site))
			nel_offinit  [k, i, j] = sim.nel(initial_site+1)
			nph_init     [k, i, j] = np.abs(initial_phonon_number - sim.nph(initial_site))
			condition_tot[k, i, j] = sim.finalInitParams["condition_tot"]
			condition_ph [k, i, j] = sim.finalInitParams["condition_ph"]
		print(f"{i+1}/{num_repetitions}")


axlabel = {"norm" : r"$|1-\langle \psi|\psi\rangle|$",
           "condP": r"$\kappa(R)$",
           "condT": r"$\kappa(S)$",
           "nel"  : r"$|\langle \psi_{init}|\hat n_{el, %i}|\psi_{init}\rangle - \langle \psi|\hat n_{el, %i}|\psi\rangle|$"%(initial_site+1, initial_site+1),
           #"nph"  : r"$|\langle \psi_{init}|\hat N_{ph, %i}|\psi_{init}\rangle - \langle \psi|\hat N_{ph, %i}|\psi\rangle|$"%(initial_site+1, initial_site+1),\
           "nph"  : r"$\Delta \langle \hat N_{ph, %i}\rangle$"%(initial_site+1),\
           "scale_ph":r"$\sigma_{ph}$"}

figsize = (3, 3)

fig_norm , ax_norm  = plt.subplots(1, 1, figsize=figsize)
fig_condP, ax_condP = plt.subplots(1, 1, figsize=figsize)
fig_condT, ax_condT = plt.subplots(1, 1, figsize=figsize)
fig_nph  , ax_nph   = plt.subplots(1, 1, figsize=figsize)
fig_nel  , ax_nel   = plt.subplots(1, 1, figsize=figsize)

axs = [ax_norm, ax_condP, ax_condT, ax_nph, ax_nel]

for r in range(num_repetitions):
	for ax in axs:
		ax.set_prop_cycle(None)
	
	ax_norm .plot(scale_phs    , norms        [:, r].T, '.-')
	ax_condP.plot(norms[:, r].T    , condition_tot[:, r].T, line_fmt)
	ax_condT.plot(norms[:, r].T    , condition_ph [:, r].T, line_fmt)
	ax_nph  .plot(norms[:, r].T, nph_init     [:, r].T, line_fmt)
	ax_nel  .plot(norms[:, r].T, nel_init     [:, r].T, line_fmt)
	

for ax in axs:
	ax.set_xscale("log")
	ax.set_yscale("log")
	
	ax.set_xlabel(axlabel["norm"])

ax_norm.set_xlabel(axlabel["scale_ph"])
ax_norm.set_ylabel(axlabel["norm"])
ax_norm.legend(num_configs_label_list, title=r"$N_c=$")

ax_condP.set_ylabel(axlabel["condP"])
ax_condT.set_ylabel(axlabel["condT"])
ax_nph  .set_ylabel(axlabel["nph"])
ax_nel  .set_ylabel(axlabel["nel"])

fig_norm.tight_layout()
fig_condP.tight_layout()
fig_condT.tight_layout()
fig_nph  .tight_layout()
fig_nel  .tight_layout()

if save_plots:
	sim_name = sim_name[:sim_name.rfind("_")]
	save_path += sim_name + f"_eigen={local_psi_init_select_eigenstate}/"
	os.system("mkdir -p " + save_path.replace("(", "\(").replace(")", "\)"))
	
	fig_norm .savefig(save_path + "norm_scale.pdf")
	fig_condP.savefig(save_path + "condP_norm.pdf")
	fig_condT.savefig(save_path + "condT_norm.pdf")
	fig_nph  .savefig(save_path + "nph_norm.pdf")
	fig_nel  .savefig(save_path + "nel_norm.pdf")

