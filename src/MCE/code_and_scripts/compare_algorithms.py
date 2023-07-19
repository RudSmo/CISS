from cycler import cycler
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


hopping = 0.5
omega = 1
gamma = 2
num_sites_x = 2
num_sites_y = 1

initial_state = "dressed"

#hopping, gamma, initial_state = 0.5, 2. , "bare"
#hopping, gamma, initial_state = 0.5, 2. , "dressed"
#hopping, gamma, initial_state = 2.5, 2.5, "bare"
#hopping, gamma, initial_state = 2.5, 2.5, "dressed"
#hopping, gamma, initial_state = 10 , 4  , "bare"
#hopping, gamma, initial_state = 10 , 4  , "dressed"


psi_initial = "local"
psi_state = (0, (num_sites_x-1)//2)

# int, iterable or "all"
show_site = psi_state[1]
show_site = "all"
MCE1pancake = False

if num_sites_x == 2:
	num_samples  = 4000
	num_configs  = 200
else:
	if hopping == 0.5:
		num_samples  = 10000
		num_configs  =  3000
	else:
		num_samples  =  4000
		num_configs  =  6000
num_configs1       = num_configs
num_configs2_local = num_configs
num_configs2_localFlip = num_configs
num_configs2_eigen = num_configs//2
num_configs2_eigenFlip = num_configs
num_configs_AIMC_init = num_configs
num_configs_AIMC_max = num_configs_AIMC_init*2

local_psi_init_select_eigenstate = None	# None for default values

if hopping == 10 and (initial_state == "bare" or num_sites_x == 2):
	show_t_max = 10
else:
	show_t_max = 20
	show_t_max = 320 #100*np.pi
sim_t_max = show_t_max
dt = 0.01 if sim_t_max < 100 else 0.1
#show_t_max = sim_t_max

print(f"{dt=}, {sim_t_max=}, {show_t_max=}")

save_data = True
force_new_simulation = False
show_plots = True
save_plots = False
save_path = "compare_plots/"	# each simulation creates a subfolder
shared_color_params = {"shading" : "nearest", "cmap" : "viridis"}
#colormap_xp = "bwr"

min_nel_cut = 1e-4
min_nph_cut = 1e-4

printProgress = True


figwide_width = {1:.7*6, 2:6, 3:6, 4:6}
figwide_height = 3
figheight_cb = 0.7

num_plots = 2
figsize = (figwide_width[num_plots]/num_plots, figwide_height)
tight_layout_pad = 0.7

# order in this dict also determines drawing oder. last in drawn on top
# either MTE or MCEv1 must be included to generate path name
show_simulation = {
                   "AIMC"      : False,
                   "MCE2eigenFlip":False,
                   "MCE2eigen" : False,
                   "MCE2localFlip" : True,
                   "MCE2local" : False,
                   "EHRENFEST" : hopping == 10 or num_sites_x != 2,
                   "MTE"       : False,
                   "MCE1"      : True,
                   "EXACT"     : True,
                   "long_exact_tunneling" : num_sites_x == 2 and initial_state == "dressed"\
                                            and hopping == 0.5 and omega == 1 and gamma == 2,
                   }


line_colors = {
               "AIMC"      : 'tab:gray',
               "MCE2eigenFlip": 'k',
               "MCE2eigen" : 'tab:brown',
               "MCE2localFlip": 'tab:red',
               "MCE2local" : 'tab:purple',
               "EHRENFEST" : 'tab:red',
               "MCE1"      : 'tab:orange',
               "EXACT"     : 'tab:blue',
               "MTE"       : 'tab:green',
               }

line_styles = {
               "AIMC"      : '-',
               "MCE2eigenFlip" : '-.',
               "MCE2eigen" : '-',
               "MCE2localFlip" : '-',
               "MCE2local" : '-',
               "EHRENFEST" : '-',
               "MCE1"      : '-',
               "EXACT"     : ':',
               "MTE"       : '--',
               }

# set linestyle as default propcycle
#default_cycler = (cycler(color=[line_colors[sim_name] for sim_name in show_simulation.keys() if show_simulation[sim_name]]) +
#                  cycler(linestyle=[line_styles[sim_name] for sim_name in show_simulation.keys() if show_simulation[sim_name]]))
#plt.rc('axes', prop_cycle=default_cycler)

sim_labels = {"EXACT"     : "ED" if num_sites_x <= 3 else "DMRG",
              "EHRENFEST" : "Ehrenfest",
              "MTE"       : "MTE",
              "MCE1"      : "MCEv1",
              "MCE2local" : "MCEv2 local",
              "MCE2localFlip" : "MCEv2 local Flip",
              "MCE2eigen" : "MCEv2 eigen",
              "MCE2eigenFlip" : "MCEv2 eigen Flip",
              "AIMC"      : "AIMC"}
time_label = r"$\omega_0 t$"

# more labels, site depend labels follow later in the loop over sites
num_sigmas = 2
ax_mean_label = {}
ax_std_label = {}
ax_mixed_label = {}

etot_label = r"$E_{tot}/\hbar\omega_0$"
etot_mean_label = r"$\langle E_{tot}\rangle/\hbar\omega_0$"

ax_mean_label ["OCDW"]   = r"$\langle\mathcal O_{\rm CDW}\rangle$"
ax_std_label  ["OCDW"]   = r"$\sigma_{\mathcal O_{\rm CDW}}$"
ax_mixed_label["OCDW"]   = r"$\langle\mathcal O_{\rm CDW}\rangle \pm %g\sigma_{\mathcal O_{\rm CDW}}$"%(num_sigmas)

ax_mean_label ["e_kin"]   = r"$\langle E_{kin}\rangle/\hbar\omega_0$"
ax_std_label  ["e_kin"]   = r"$\langle E_{kin}\rangle/\hbar\omega_0$"
ax_mixed_label["e_kin"]   = r"$\left(\langle E_{kin}\rangle\pm %g\sigma_{\langle E_{kin}\rangle}\right)/\hbar\omega_0$"\
                            %(num_sigmas)

used_sim_names = [sim_name for sim_name in show_simulation.keys() if show_simulation[sim_name]]

simulation = dict()

for sim_name in used_sim_names:
	if sim_name == "EXACT":
		print("loading exact data")
		try:
			simulation["EXACT"] = data_manager.RefData(initial_state, psi_initial, psi_state, num_sites_x, num_sites_y, hopping, omega, gamma)
		except Exception:
			print("#######failed loading exact data######")
			show_simulation["EXACT"] = False
	elif sim_name == "EHRENFEST":
		print("running Ehrenfest")
		simulation["EHRENFEST"] = data_manager.Simulation_MTE(hopping, omega, gamma,
		                                                      num_sites_x, num_sites_y,
		                                                      sim_t_max, dt,
		                                                      1,	#num_samples
		                                                      initial_state, psi_initial, psi_state,
		                                                      phonon_dressed_displaced_sites=None,
		                                                      calcYSquared=True,
		                                                      fullStorage=False,
		                                                      printProgress=printProgress,
		                                                      symmetric_y_init=False,
		                                                      force_new_simulation=force_new_simulation,
		                                                      load_sampled_state=False,
		                                                      store_results=save_data)
	elif sim_name == "MTE":
		print("running MTE")
		simulation["MTE"] = data_manager.Simulation_MTE(hopping, omega, gamma,
		                                                num_sites_x, num_sites_y,
		                                                sim_t_max, dt,
		                                                num_samples,
		                                                initial_state, psi_initial, psi_state,
		                                                phonon_dressed_displaced_sites=None,
		                                                calcYSquared=True,
		                                                fullStorage=False,
		                                                printProgress=printProgress,
		                                                symmetric_y_init=False,
		                                                force_new_simulation=force_new_simulation,
		                                                load_sampled_state=False,
		                                                store_results=save_data)
		shared_path = save_path
		full_sim_name = simulation["MTE"].get_sim_name()
		cropped_sim_name = full_sim_name[:full_sim_name.rfind("_numtraj=")]
		shared_path += cropped_sim_name[4:].replace(" ", "")
	elif sim_name == "MCE1":
		print("running MCE1")
		simulation["MCE1"] = data_manager.Simulation_MCE(1, hopping, omega, gamma,
		                                                 num_sites_x, num_sites_y,
		                                                 sim_t_max, dt,
		                                                 num_configs1,
		                                                 initial_state, psi_initial, psi_state,
		                                                 num_configs1,
		                                                 allow_config_delete=False,
		                                                 allow_config_clone=False,
		                                                 printProgress=printProgress,
		                                                 force_new_simulation=force_new_simulation,
		                                                 load_sampled_state=False,
		                                                 store_results=save_data,
		                                                 local_psi_init_select_eigenstate=False,
		                                                 pancakeSampling=MCE1pancake)
		shared_path = save_path
		full_sim_name = simulation["MCE1"].get_sim_name()
		cropped_sim_name = full_sim_name[:full_sim_name.rfind("_numconf=")]
		shared_path += cropped_sim_name[6:].replace(" ", "")
	elif "MCE2" in sim_name:
		print("running", sim_name)
		if sim_name == "MCE2local":
			used_num_configs = num_configs2_local
		if sim_name == "MCE2localFlip":
			used_num_configs = num_configs2_localFlip
		if sim_name == "MCE2eigen":
			used_num_configs = num_configs2_eigen
		if sim_name == "MCE2eigenFlip":
			used_num_configs = num_configs2_eigenFlip
		simulation[sim_name] = data_manager.Simulation_MCE(2, hopping, omega, gamma,
		                                                   num_sites_x, num_sites_y,
		                                                   sim_t_max, dt,
		                                                   used_num_configs,
		                                                   initial_state, psi_initial, psi_state,
		                                                   used_num_configs,
		                                                   allow_config_delete=False,
		                                                   allow_config_clone=False,
		                                                   printProgress=printProgress,
		                                                   force_new_simulation=force_new_simulation,
		                                                   load_sampled_state=False,
		                                                   store_results=save_data,
		                                                   local_psi_init_select_eigenstate="eigen" in sim_name,
		                                                   flip_initial_state_for_half_configs="Flip" in sim_name)
	elif sim_name == "AIMC":
		print("running AIMC")
		simulation["AIMC"] = data_manager.Simulation_MCE(2, hopping, omega, gamma,
		                                                 num_sites_x, num_sites_y,
		                                                 sim_t_max, dt,
		                                                 num_configs_AIMC_init,
		                                                 initial_state, psi_initial, psi_state,
		                                                 num_configs_AIMC_max,
		                                                 allow_config_delete=False,
		                                                 allow_config_clone=True,
		                                                 printProgress=printProgress,
		                                                 force_new_simulation=force_new_simulation,
		                                                 load_sampled_state=False,
		                                                 store_results=save_data,
		                                                 local_psi_init_select_eigenstate=True)

#if "EXACT" in used_sim_names and show_simulation["EXACT"]:
used_sim_names = [sim_name for sim_name in show_simulation.keys() if show_simulation[sim_name]]

if show_simulation["EHRENFEST"]:
	shared_path += "_EHRF"
if show_simulation["MTE"]:
	shared_path += f"_MTE={num_samples}"
if show_simulation["MCE1"]:
	shared_path += f"_MCE1={num_configs1}"
	# pancake is already included in simualtion name, no need to add it here
if show_simulation["MCE2local"]:
	shared_path += f"_MCE2local={num_configs2_local}"
if show_simulation["MCE2localFlip"]:
	shared_path += f"_MCE2localFlip={num_configs2_local}"
if show_simulation["MCE2eigen"]:
	shared_path += f"_MCE2eigen={num_configs2_eigen}"
if show_simulation["MCE2eigenFlip"]:
	shared_path += f"_MCE2eigenFlip={num_configs2_eigen}"
if show_simulation["AIMC"]:
	shared_path += f"_AIMC={num_configs_AIMC_init}-{num_configs_AIMC_max}"
if show_simulation["long_exact_tunneling"]:
	shared_path+= "_long"

shared_path += "/"

if save_plots:
	os.system("mkdir -p " + shared_path.replace("(", r"\(").replace(")", r"\)"))
	print("saving to:", shared_path)


if show_site == "all":
	sites_to_draw = np.arange(num_sites_x)
else:
	sites_to_draw = np.array(show_site, ndmin=1)

for sim_name in used_sim_names:
	if sim_name != "long_exact_tunneling":
		sim = simulation[sim_name]
		sim.set_cut_time(show_t_max)

plot_inset = "EHRENFEST" in used_sim_names and num_sites_x == 2 and hopping == 10

for site in sites_to_draw:
	ax_mean_label ["nel"] = r"$\langle \hat{n}_{el,%i}\rangle_r$"%(site+1)
	ax_std_label  ["nel"] = r"$\sigma_{\langle \hat n_{el,%i}\rangle}$"%(site+1)
	ax_mixed_label["nel"] = r"$\langle\hat n_{el,%i}\rangle_r \pm %g\sigma_{\langle\hat n_{el,%i}\rangle}$"%(site+1, num_sigmas, site+1)
	
	ax_mean_label ["nph"] = r"$\langle\hat N_{ph,%i}\rangle_r$"%(site+1)
	ax_std_label  ["nph"] = r"$\sigma_{\langle\hat N_{ph,%i}\rangle}$"%(site+1)
	ax_mixed_label["nph"] = r"$\langle\hat N_{ph,%i}\rangle_r \pm %g\sigma_{\langle\hat N_{ph,%i}\rangle}$"%(site+1, num_sigmas, site+1)
	
	ax_mean_label ["x"]       = r"$\langle\hat x_{%i}\rangle_r/\sqrt{\frac{\hbar}{m\omega_0}}$"%(site+1)
	ax_std_label  ["x"]       = r"$\sigma_{\langle\hat x_{%i}\rangle}/\sqrt{\frac{\hbar}{m\omega_0}}$"%(site+1)
	ax_mixed_label["x"]       = r"$\left(\langle\hat x_{%i}\rangle_r \pm %g\sigma_{\langle\hat x_{%i}\rangle}\right)/\sqrt{\frac{\hbar}{m\omega_0}}$"\
	                            %(site+1, num_sigmas, site+1)
	
	ax_mean_label ["p"]       = r"$\langle\hat p_{%i}\rangle_r/\sqrt{\hbar m\omega_0}$"%(site+1)
	ax_std_label  ["p"]       = r"$\sigma_{\langle\hat p_{%i}\rangle}/\sqrt{\hbar m\omega_0}$"%(site+1)
	ax_mixed_label["p"]       = r"$\left(\langle\hat p_{%i}\rangle_r \pm %g\sigma_{\langle\hat p_{%i}\rangle}\right)/\sqrt{\hbar m\omega_0}$"\
	                            %(site+1, num_sigmas, site+1)
	
	# replace by strings without _r means, TODO get means to work and use them
	
	ax_mean_label ["nel"] = r"$\langle \hat{n}_{el,%i}\rangle$"%(site+1)
	ax_std_label  ["nel"] = r"$\sigma_{\langle \hat n_{el,%i}\rangle}$"%(site+1)
	ax_mixed_label["nel"] = r"$\langle\hat n_{el,%i}\rangle \pm %g\sigma_{\langle\hat n_{el,%i}\rangle}$"%(site+1, num_sigmas, site+1)
	
	ax_mean_label ["nph"] = r"$\langle\hat N_{ph,%i}\rangle$"%(site+1)
	ax_std_label  ["nph"] = r"$\sigma_{\langle\hat N_{ph,%i}\rangle}$"%(site+1)
	ax_mixed_label["nph"] = r"$\langle\hat N_{ph,%i}\rangle \pm %g\sigma_{\langle\hat N_{ph,%i}\rangle}$"%(site+1, num_sigmas, site+1)
	
	ax_mean_label ["x"]       = r"$\langle\hat x_{%i}\rangle/\sqrt{\frac{\hbar}{m\omega_0}}$"%(site+1)
	ax_std_label  ["x"]       = r"$\sigma_{\langle\hat x_{%i}\rangle}/\sqrt{\frac{\hbar}{m\omega_0}}$"%(site+1)
	ax_mixed_label["x"]       = r"$\left(\langle\hat x_{%i}\rangle \pm %g\sigma_{\langle\hat x_{%i}\rangle}\right)/\sqrt{\frac{\hbar}{m\omega_0}}$"\
	                            %(site+1, num_sigmas, site+1)
	
	ax_mean_label ["p"]       = r"$\langle\hat p_{%i}\rangle/\sqrt{\hbar m\omega_0}$"%(site+1)
	ax_std_label  ["p"]       = r"$\sigma_{\langle\hat p_{%i}\rangle}/\sqrt{\hbar m\omega_0}$"%(site+1)
	ax_mixed_label["p"]       = r"$\left(\langle\hat p_{%i}\rangle \pm %g\sigma_{\langle\hat p_{%i}\rangle}\right)/\sqrt{\hbar m\omega_0}$"\
	                            %(site+1, num_sigmas, site+1)
	

	fig_nel, ax_nel = plt.subplots(1, 1, figsize=figsize)
	if plot_inset:
		fig_nel_inset, ax_nel_inset = plt.subplots(1, 1, figsize=figsize)
		
		ax_nel_inset2 = fig_nel_inset.add_axes((.42,.79,.33,.165))
		ax_nel_inset2.plot(simulation["EHRENFEST"].time(), simulation["EHRENFEST"].nel(site),
		                   label=sim_labels["EHRENFEST"], linestyle=line_styles["EHRENFEST"],
		                   color=line_colors["EHRENFEST"], linewidth=0.5)
		
		ax_nel_inset.plot([], [], label=sim_labels["EHRENFEST"],
			              linestyle=line_styles["EHRENFEST"], color=line_colors["EHRENFEST"])
		
	fig_nph, ax_nph = plt.subplots(1, 1, figsize=figsize)
	fig_x, ax_x = plt.subplots(1, 1, figsize=figsize)
	fig_p, ax_p = plt.subplots(1, 1, figsize=figsize)
	
	
	for sim_name in used_sim_names:
		if sim_name == "long_exact_tunneling":
			data = np.loadtxt("ref_data/long_exact_dressed_local(0,1)_Lx=2_t0=0.5_omega=1_gamma=2_nel.txt")
			data[:, 0] /= hopping
			data = data[:np.argmax(data[:, 0] > show_t_max)]
			ax_nel.plot(data[:, 0], data[:, site+3], label=sim_labels["EXACT"],
			            linestyle=line_styles["EXACT"], color=line_colors["EXACT"])
			
			plt.figure()
			plt.title("nel diff " + str(site))
			plt.plot(simulation["MCE1"].time()[::2],  data[:, site+3] - simulation["MCE1"].nel(site)[::2], label="MCE1")
			plt.plot(simulation["MCE2localFlip"].time()[::2],  data[:, site+3] - simulation["MCE2localFlip"].nel(site)[::2], label="MCE2")
			plt.legend()
		else:
			sim = simulation[sim_name]
			
			if "long_exact_tunneling" not in used_sim_names or sim_name != "EXACT":
				ax_nel.plot(sim.time(), sim.nel(site), label=sim_labels[sim_name],
				            linestyle=line_styles[sim_name], color=line_colors[sim_name])
			if plot_inset and sim_name != "EHRENFEST":
				ax_nel_inset.plot(sim.time(), sim.nel(site), label=sim_labels[sim_name],
				                  linestyle=line_styles[sim_name], color=line_colors[sim_name])
			ax_nph.plot(sim.time(), sim.nph(site), label=sim_labels[sim_name],
			            linestyle=line_styles[sim_name], color=line_colors[sim_name])
			ax_x  .plot(sim.time(), sim.x  (site), label=sim_labels[sim_name],
			            linestyle=line_styles[sim_name], color=line_colors[sim_name])
			ax_p  .plot(sim.time(), sim.p  (site), label=sim_labels[sim_name],
			            linestyle=line_styles[sim_name], color=line_colors[sim_name])
	
	if plot_inset:
		
		ax_nel_inset.set_xlabel(time_label)
		ax_nel_inset.set_ylabel(ax_mean_label["nel"])
		
		#ax_nel_inset2.set_xlabel(time_label)
		#ax_nel_inset2.set_ylabel(r"$\langle n_{el,%i}\rangle$"%(site+1))
		
		fig_nel_inset.tight_layout(pad=tight_layout_pad)
		
		if save_plots:
			fig_nel_inset.savefig(shared_path + f"nel{site+1}_inset.pdf")
		ax_nel_inset.legend(loc="lower center")
		if save_plots:
			fig_nel_inset.savefig(shared_path + f"nel{site+1}_inset_legend.pdf")
		
		if not show_plots:
			plt.close(fig_nel_inset)
	
	ax_nel.set_xlabel(time_label)
	ax_nel.set_ylabel(ax_mean_label["nel"])
	
	ax_nph.set_xlabel(time_label)
	ax_nph.set_ylabel(ax_mean_label["nph"])
	
	ax_x.set_xlabel(time_label)
	ax_x.set_ylabel(ax_mean_label["x"])
	
	ax_p.set_xlabel(time_label)
	ax_p.set_ylabel(ax_mean_label["p"])
	
	fig_nel.tight_layout(pad=tight_layout_pad)
	fig_nph.tight_layout(pad=tight_layout_pad)
	fig_x  .tight_layout(pad=tight_layout_pad)
	fig_p  .tight_layout(pad=tight_layout_pad)
	
	if save_plots:
		fig_nel.savefig(shared_path + f"nel{site+1}.pdf")
		fig_nph.savefig(shared_path + f"nph{site+1}.pdf")
		fig_x  .savefig(shared_path + f"x{site+1}.pdf")
		fig_p  .savefig(shared_path + f"p{site+1}.pdf")
	
	ax_nel.legend()
	ax_nph.legend()
	ax_x.legend()
	ax_p.legend()
	
	if save_plots:
		fig_nel.savefig(shared_path + f"nel{site+1}_legend.pdf")
		fig_nph.savefig(shared_path + f"nph{site+1}_legend.pdf")
		fig_x  .savefig(shared_path + f"x{site+1}_legend.pdf")
		fig_p  .savefig(shared_path + f"p{site+1}_legend.pdf")
	if not show_plots:
		plt.close(fig_nel)
		plt.close(fig_nph)
		plt.close(fig_x)
		plt.close(fig_p)


# color plots
nel_params = dict(shared_color_params)
nph_params = dict(shared_color_params)

max_nel = max([simulation[sim_name].nel().max() for sim_name in used_sim_names if sim_name not in ["EHRENFEST", "long_exact_tunneling"]])
max_nel_after_wt = max([simulation[sim_name].nel()[int(1./dt):].max() for sim_name in used_sim_names if sim_name not in ["EHRENFEST", "long_exact_tunneling"]])
min_nel = min([simulation[sim_name].nel().min() for sim_name in used_sim_names if sim_name not in ["EHRENFEST", "long_exact_tunneling"]])

max_nph = max([simulation[sim_name].nph().max() for sim_name in used_sim_names if sim_name not in ["long_exact_tunneling"]])
min_nph = min([simulation[sim_name].nph().min() for sim_name in used_sim_names if sim_name not in ["long_exact_tunneling"]])

lattice_sites = np.arange(num_sites_x)+1
tick_spacing = max(1, int(np.floor(num_sites_x/4)))
lattice_ticks = lattice_sites[((num_sites_x-1)//2)%tick_spacing::tick_spacing]

for use_log_color in ([False, True] if num_sites_x > 2 else []):
	if use_log_color:
		min_nel = max(min_nel, min_nel_cut)
		min_nph = max(min_nph, min_nph_cut)
		
		nel_params["norm"] = matplotlib.colors.LogNorm(vmin=min_nel, vmax=max_nel, clip=True)
		nph_params["norm"] = matplotlib.colors.LogNorm(vmin=min_nph, vmax=max_nph, clip=True)
	else:
		if max_nel_after_wt < 0.5 * max_nel:
			max_nel = np.ceil(max_nel_after_wt*10)/10
			
		nel_params["norm"] = matplotlib.colors.Normalize(vmin=min_nel, vmax=max_nel, clip=True)
		nph_params["norm"] = matplotlib.colors.Normalize(vmin=min_nph, vmax=max_nph, clip=True)
	
	
	
	fig_nel_shared, ax_nel_shared = plt.subplots(1, len(used_sim_names),
	                                             figsize=(figwide_width[2], figwide_height))
	fig_nph_shared, ax_nph_shared = plt.subplots(1, len(used_sim_names),
	                                             figsize=(figwide_width[2], figwide_height))
	
	for sim_num in range(len(used_sim_names)):
		sim_name = used_sim_names[sim_num]
		
		if sim_name == "long_exact_tunneling":
			continue
		
		fig_nel, ax_nel = plt.subplots(1, 1, figsize=figsize)
		fig_nph, ax_nph = plt.subplots(1, 1, figsize=figsize)
		fig_x  , ax_x   = plt.subplots(1, 1, figsize=figsize)
		fig_p  , ax_p   = plt.subplots(1, 1, figsize=figsize)
		sim = simulation[sim_name]
		
		if use_log_color:
			nel_cutted = np.where(sim.nel() >= min_nel_cut, sim.nel(), min_nel_cut)
			nph_cutted = np.where(sim.nph() >= min_nph_cut, sim.nph(), min_nph_cut)
		else:
			nel_cutted = sim.nel()
			nph_cutted = sim.nph()
		
		cm_nel = ax_nel.pcolormesh(lattice_sites, sim.time(), nel_cutted, **nel_params)
		ax_nel_shared[sim_num].pcolormesh(lattice_sites, sim.time(), nel_cutted, **nel_params)
		cm_nph = ax_nph.pcolormesh(lattice_sites, sim.time(), nph_cutted, **nph_params)
		ax_nph_shared[sim_num].pcolormesh(lattice_sites, sim.time(), nph_cutted, **nph_params)
		#cm_x   = ax_nph.pcolormesh(lattice_sites, sim.time(), sim.x(), **x_params)
		#cm_p   = ax_nph.pcolormesh(lattice_sites, sim.time(), sim.p(), **p_params)
		
		
		ax_nel.set_xlabel("site i")
		ax_nel.set_xticks(lattice_ticks)
		ax_nel.set_ylim((0, show_t_max))
		ax_nel.set_ylabel(time_label)
		#ax_nel.set_ylabel(r"$\langle n_{el,%i}\rangle$"%(site+1))
		
		ax_nph.set_xlabel("site i")
		ax_nph.set_xticks(lattice_ticks)
		ax_nph.set_ylim((0, show_t_max))
		ax_nph.set_ylabel(time_label)
		#ax_nph.set_ylabel(r"$\langle N_{ph,%i}\rangle$"%(site+1))
		
		fig_nel.tight_layout(pad=tight_layout_pad)
		fig_nph.tight_layout(pad=tight_layout_pad)
		
		if save_plots:
			name_infix = "_color_"
			if use_log_color:
				name_infix = "_log" + name_infix
			fig_nel.savefig(shared_path + "nel" + name_infix + f"{sim_name}.png", dpi=600)
			fig_nph.savefig(shared_path + "nph" + name_infix + f"{sim_name}.png", dpi=600)
		if not show_plots:
			plt.close(fig_nel)
			plt.close(fig_nph)
		
		ax_nel_shared[sim_num].set_title(sim_labels[sim_name])
		ax_nel_shared[sim_num].set_xlabel("site i")
		ax_nel_shared[sim_num].set_xticks(lattice_ticks)
		ax_nel_shared[sim_num].set_ylim((0, show_t_max))
		
		ax_nph_shared[sim_num].set_title(sim_labels[sim_name])
		ax_nph_shared[sim_num].set_xlabel("site i")
		ax_nph_shared[sim_num].set_xticks(lattice_ticks)
		ax_nph_shared[sim_num].set_ylim((0, show_t_max))
		
		if sim_num != 0:
			ax_nel_shared[sim_num].set_yticklabels([])
			ax_nph_shared[sim_num].set_yticklabels([])
	
	ax_nel_shared[0].set_ylabel(time_label)
	ax_nph_shared[0].set_ylabel(time_label)
	
	if save_plots:
		name_suffix = "_color.png"
		if use_log_color:
			name_suffix = "_log" + name_suffix
		
		name_sufffix_colorbar = "_colorbar.png"
		if use_log_color:
			name_sufffix_colorbar = "_log" + name_sufffix_colorbar
		
		nel_sm = matplotlib.cm.ScalarMappable(norm=nel_params["norm"], cmap=nel_params["cmap"])
		nph_sm = matplotlib.cm.ScalarMappable(norm=nph_params["norm"], cmap=nph_params["cmap"])
		
		subplot_placement_shared = (0.1, 0.14, 0.87, 0.92, 0.03)
		cb_box = (0.88, subplot_placement_shared[1], 0.02, subplot_placement_shared[3]-subplot_placement_shared[1])
		fig_nel_shared.subplots_adjust(*subplot_placement_shared)
		fig_nph_shared.subplots_adjust(*subplot_placement_shared)
		
		ax_cb_nel_shared = fig_nel_shared.add_axes(cb_box)
		fig_nel_shared.colorbar(nel_sm, cax=ax_cb_nel_shared, orientation="vertical", label=r"$\langle\hat n_{el, i}\rangle$")
		ax_cb_nph_shared = fig_nph_shared.add_axes(cb_box)
		fig_nph_shared.colorbar(nph_sm, cax=ax_cb_nph_shared, orientation="vertical", label=r"$\langle\hat N_{ph, i}\rangle$")
		
		
		fig_nel_shared.savefig(shared_path + "nel" + name_suffix, dpi=600)
		fig_nph_shared.savefig(shared_path + "nph" + name_suffix, dpi=600)
		
		
		fig_cb_nel, ax_cb_nel = plt.subplots(figsize=(figwide_width[1], figheight_cb))
		fig_cb_nel.colorbar(nel_sm, cax=ax_cb_nel, orientation="horizontal", label=r"$\langle\hat n_{el, i}\rangle$")
		fig_cb_nel.tight_layout(pad=tight_layout_pad)
		fig_cb_nel.savefig(shared_path + "nel" + name_sufffix_colorbar, dpi=600)
		
		fig_cb_nph, ax_cb_nph = plt.subplots(figsize=(figwide_width[1], figheight_cb))
		fig_cb_nph.colorbar(nph_sm, cax=ax_cb_nph, orientation="horizontal", label=r"$\langle\hat N_{ph, i}\rangle$")
		fig_cb_nph.tight_layout(pad=tight_layout_pad)
		fig_cb_nph.savefig(shared_path + "nph" + name_sufffix_colorbar, dpi=600)
	
	if not show_plots:
		plt.close(fig_cb_nel)
		plt.close(fig_cb_nph)
		plt.close(fig_nel_shared)
		plt.close(fig_nph_shared)


fig_ekin, ax_ekin = plt.subplots(1, 1, figsize=figsize)
ax_ekin.set_xlabel(time_label)
ax_ekin.set_ylabel(ax_mean_label["e_kin"])
for sim_num in range(len(used_sim_names)):
	sim_name = used_sim_names[sim_num]
	if sim_name == "long_exact_tunneling":
		continue
	ax_ekin.plot(sim.time(), sim.energy_kinetic(), label=sim_labels[sim_name])

fig_ekin.tight_layout(pad=tight_layout_pad)
if save_plots:
	fig_ekin.savefig(shared_path + "energy_kinetic.pdf")
	
ax_ekin.legend()
if save_plots:
	fig_ekin.savefig(shared_path + "energy_kinetic_legend.pdf")

#TODO use cm_n[el,ph] to generate colorbars

if show_plots:
	plt.show()
