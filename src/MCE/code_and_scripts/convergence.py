from cycler import cycler
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

# in prop_cycle, skip pink and make purple, brown and gray dashed for color blind
default_cycler = (cycler(color=['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:gray']) +
                  cycler(linestyle=['-', '-', '-', '--', '--', ':', ':']))
plt.rc('axes', prop_cycle=default_cycler)

exact_line_style = ":"
exact_color = "k"

version = 1 # 1, 2 (both MCE), "MTE" (so far mine) or "m_MTE", "m_SH", "m_SHDC" (michael)
hopping = 0.5
omega = 1
gamma = 2
num_sites_x = 11
num_sites_y = 1

large_convergence_set = False
compare_show_exact = hopping != 0.5 and True

allow_crop_nel_range = True
crop_nel_threshold = 0.3
crop_nel_start_time = 1.

pancakeSampling = False
pancakeScale = 1.

initial_state = "bare" # "bare" or "dressed"
psi_initial = "local"

init_site = (num_sites_x-1)//2
site = init_site + 0		# shown site

save_plots_compare = True
save_plots_single_runs = not large_convergence_set and True
save_plots_std_nc_dependence = large_convergence_set or False
#cut_times_nph_std_and_energy_vs_nc = [15.8]
#allow_lag_scale_num_configs = True
#log_scale_min_max_ratio_num_configs = 
cut_times_nph_std_and_energy_vs_nc = np.arange(1, 21, 1)
allow_show_mean_diff = True
plot_max_energy_curve = False

local_psi_init_select_eigenstate = False	# None for default values
scale_el = None
allow_config_delete = False
allow_config_clone  = False

allowCompression = True
allowExpansion = True

if version == "MTE":
	if num_sites_x == 2:
		num_config_list = [100, 400, 1000, 2000]#, 1000, 2000]
	else:
		num_config_list = [100, 400, 1000, 2000]
	#num_config_list = [4, 8, 12]
elif version in [1, 2]:
	if num_sites_x == 2:
		num_config_list = [20, 40, 60, 80]
		num_config_list = [1, 2, 5, 10, 20, 40, 80]
		num_config_list = [1, 2, 3, 4, 5, 7, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120]
	else:
		num_config_list = [200, 500, 1000]
		if large_convergence_set:
			num_config_list = np.arange(100, 3001, 100)
			num_config_list = np.concatenate((num_config_list, [1]))
			num_config_list = np.concatenate((num_config_list, [4000, 5000, 6000]))
			num_config_list = np.concatenate((num_config_list, np.arange(10, 100, 10)))
			num_config_list = np.concatenate((num_config_list, np.arange(100, 200, 25)))
		elif hopping != 0.5:
			num_config_list = [1000, 1500, 3000, 4000, 5000, 6000]
		else:
			num_config_list = [1000, 1500, 3000]
		
		num_config_list = list(np.sort(np.unique(num_config_list)))
		
		
	#num_config_list = [4, 8, 12]#, 16, 20]
else: # Michaels data
	if version == "m_MTE":
		num_config_list = [50, 100, 200, 400]
	else:
		num_config_list = [200, 400, 1000, 2000]

if version == 2 and (allow_config_clone or allow_config_delete):
	max_cloning_ratio = 1.5
else:
	max_cloning_ratio = 1

assert(max_cloning_ratio >= 1)

num_config_list = np.sort(np.unique(np.array(num_config_list)))

varying_num_reps = len(num_config_list) >= 8

if type(version) == str and "m_" == version[:2]:
	num_repetitions = np.full(len(num_config_list), 50, dtype=np.int)
	t_max = 20
	dt = 0.01
	show_dist = 1
	load_mode = True
	obs_list = ["nel", "nel_adia", "x"]
else:
	if num_sites_x == 2:
		if varying_num_reps:
			num_repetitions = np.where(num_config_list > 120, 10, 50)
			num_repetitions[num_config_list <= 20] = 100
		else:
			num_repetitions = np.full(len(num_config_list), 50, dtype=np.int)
	else:
		num_repetitions = np.full(len(num_config_list), 10, dtype=np.int)
		if varying_num_reps:
			num_repetitions = np.full(len(num_config_list), 10, dtype=np.int)
			num_repetitions[np.array(num_config_list) < 2000] = 20
			num_repetitions[np.array(num_config_list) < 1500] = 30
			num_repetitions[np.array(num_config_list) <= 500] = 50
	num_repetitions[num_config_list == 1] = 2	# 2 for std
	
	split_simulation_in_half = False
	odd_indices = False # selects which index subset is used
	if split_simulation_in_half:
		num_repetitions //= 2
	
	if num_sites_x == 2 and num_sites_y == 1:
		t_max = 20
		dt = 0.01
	else:
		t_max = 20 # 0.65 for DMRG range
		dt = 0.01
	show_dist = 1
	load_mode = False
	obs_list = ["nph", "x", "p", "e_kin", "e_inter"]  # should contain nel XOR OCDW, otherwise OCDW is dropped in overview plots
	if psi_initial == "CDW":
		obs_list.append("OCDW")
	else:
		obs_list.append("nel")

num_config_list = list(num_config_list)

print(f"{num_config_list=}")
print(f"{num_repetitions=}")

use_phonon_sum = psi_initial == "CDW"
use_e_inter_sum = use_phonon_sum

num_times = int(np.round(t_max/dt)) + 1
num_times = (num_times-1)//show_dist + 1


BO_state = 0 # currently only supported in load_mode
if psi_initial == "local":
	psi_state = (0, init_site) # (y, x)
elif psi_initial == "eigen":
	psi_state = 1
else: # CDW
	psi_state = None

if psi_initial == "CDW":
	if initial_state == "bare":
		short_title = "BCDW"
	else:
		short_title = "DCDW"
else:
	short_title = f"{initial_state} {psi_initial}"

short_title+= f" $L$ = {num_sites_x}, $\\bar\\gamma$ = %g, $\\bar t_0$ = %g"%(gamma/omega, hopping/omega)
short_title_rep = short_title + f", $N_r$ = {num_repetitions[-1]}"

plot_for_paper = False

save_data = True
force_new_simulation = False
show_plots = False


save_path = "convergence_plots/"
if plot_for_paper:
	subplots_in_thesis = False # must stay false to ensure correct data layout (could be changed, but isn't)
	save_path += "paper/"
	annotate_pos = (0.3, 0.85)
	ylabel_coord = (-0.08, 0.5)
else:
	subplots_in_thesis = True
	save_path += "thesis/"
	annotate_pos = (0.35, 0.85)
	ylabel_coord = None
alpha=0.3

if len(num_config_list) >= 8:
	config_list_str = "many" + str(min(num_config_list)) + "-" + str(max(num_config_list))
else:
	if max_cloning_ratio == 1:
		config_list_str = str(num_config_list)
	else:
		config_list_str = "["
		for num_configs in num_config_list:
			config_list_str += str(num_configs) + "-" + str(int(max_cloning_ratio*num_configs)) + ","
		config_list_str = config_list_str[:-1] + "]"


if subplots_in_thesis:
	legends_in_figs = ["nel", "x", "p", "OCDW"]
else:
	legends_in_figs = "all"

# for paper
figslim_width = 5
figslim_height = {1:3, 2:5, 3:7}

# for thesis:
figwide_width = {1:.7*6, 2:6, 3:6, 4:6}
figwide_height = 3

# if the subplots are saved seperatly, this value is used to get figwidth (if not None and onlt for comparison)
fixed_figwidth_count = 2		# None or positive int

num_plots_nph_std_and_energy_vs_nc = 2
figwidth_nph_std_and_energy_vs_nc = figwide_width[num_plots_nph_std_and_energy_vs_nc]/num_plots_nph_std_and_energy_vs_nc


if subplots_in_thesis:
	tight_layout_pad = 0.7
else:
	tight_layout_pad = 1.

tick_params = {"axis":"both", "direction":"in", "left":True, "bottom":True, "right":True, "top":True}


# following one doesn't affect the plot with just energy
plot_rms_etot = False
plot_total_energy = type(version) == int  # only show total energy for MCE
include_energy_to_obs_compare = plot_total_energy
printMoreProgress = True
auto_shift = True


if version == 1 and num_sites_x == 2 and num_sites_y == 1 and initial_state == "dressed"\
   and omega == 1 and hopping == 0.5 and gamma == 2:
	extra_subfolder = ["", "flip_large/", "flip_small/", "noflip_large/", "noflip_small/"][4]
else:
	extra_subfolder = ""

if auto_shift:
	allowed_shift_distances = np.array([1., 2.5, 5., 10.])  # multiplied with an arbitrary power of 10
else:
	shift_distances = {"nel"      : 0.25,
	                   "nel_adia" : 0.25,
	                   "nph"      : 2.,
	                   "x"        : 1.,
	                   "p"        : 1.,
	                   "OCDW"     : 0.25}


def simulate(num_configs, rep=None):
	if version == "MTE":
		num_samples = num_configs
		return data_manager.Simulation_MTE(hopping, omega, gamma,
		                                   num_sites_x, num_sites_y,
		                                   t_max, dt,
		                                   num_samples,
		                                   initial_state, psi_initial, psi_state,
		                                   symmetric_y_init=False,
		                                   printProgress=printMoreProgress,
		                                   calcYSquared=True,
		                                   fullStorage=False,
		                                   extra_subfolder=str(rep))
	else:
		max_num_configs = int(num_configs*max_cloning_ratio)
		return data_manager.Simulation_MCE(version, hopping, omega, gamma,
		                                   num_sites_x, num_sites_y,
		                                   t_max, dt,
		                                   num_configs,
		                                   initial_state, psi_initial, psi_state,
		                                   max_num_configs,
		                                   allow_config_delete=allow_config_delete,
		                                   allow_config_clone=allow_config_clone,
		                                   printProgress=printMoreProgress,
		                                   force_new_simulation=force_new_simulation,
		                                   load_sampled_state=False,
		                                   store_results=save_data,
		                                   extra_subfolder=extra_subfolder+str(rep),
		                                   calc_BO_electron_density=False,
		                                   local_psi_init_select_eigenstate=local_psi_init_select_eigenstate,
		                                   scale_el=scale_el,
		                                   allowCompression=allowCompression, allowExpansion=allowExpansion,
		                                   pancakeSampling=pancakeSampling, pancakeScale=pancakeScale)

ref = None
# only shown in overview, not in paper version
try:
	ref = data_manager.RefData(initial_state, psi_initial, psi_state, num_sites_x, num_sites_y, hopping, omega, gamma)
	ref.set_cut_time(t_max)
except Exception as e:
	print("no reference Data found:")
	print(repr(e))
	ref = None

exact_label = "ED" if num_sites_x <= 3 else "DMRG"

if plot_for_paper:
	time_label = r"$t\ [1/\omega_0]$"
else:
	time_label = r"$\omega_0 t$"
if omega != 1:
	print("WARNING: omega is not 1, so the time labels are wrong")


num_sigmas = 2.

ax_mean_label = {}
ax_std_label = {}
ax_mixed_label = {}

if plot_for_paper:
	ax_mean_label ["nel"] = r"$\langle n_{el,%i}\rangle$"%(site+1)
	ax_std_label  ["nel"] = r"$\sigma_{n_{el,%i}}$"%(site+1)
	ax_mixed_label["nel"] = r"$\langle n_{el,%i}\rangle \pm %g\sigma_{n_{el,%i}}$"%(site+1, num_sigmas, site+1)
	
	ax_mean_label ["nel_adia"] = r"$\langle n_{BO,%i}\rangle$"%(BO_state)
	ax_std_label  ["nel_adia"] = r"$\sigma_{n_{BO,%i}}$"%(BO_state)
	ax_mixed_label["nel_adia"] = r"$\langle n_{BO,%i}\rangle \pm %g\sigma_{n_{BO,%i}}$"%(BO_state, num_sigmas, BO_state)
	
	if use_phonon_sum:
		ax_mean_label ["nph"] = r"$\langle N_{ph}/N\rangle$"
		ax_std_label  ["nph"] = r"$\sigma_{N_{ph}/N}$"
		ax_mixed_label["nph"] = r"$\langle N_{ph}/N\rangle \pm %g\sigma_{N_{ph}/N}$"%(num_sigmas)
	else:
		ax_mean_label ["nph"] = r"$\langle N_{ph,%i}\rangle$"%(site+1)
		ax_std_label  ["nph"] = r"$\sigma_{N_{ph,%i}}$"%(site+1)
		ax_mixed_label["nph"] = r"$\langle N_{ph,%i}\rangle \pm %g\sigma_{N_{ph,%i}}$"%(site+1, num_sigmas, site+1)
	
	ax_mean_label ["x"]   = r"$\langle x_{%i}\rangle\ \left[\sqrt{\frac{\hbar}{m\omega_0}}\right]$"%(site+1)
	ax_std_label  ["x"]   = r"$\sigma_{x_{%i}}\ \left[\sqrt{\frac{\hbar}{m\omega_0}}\right]$"%(site+1)
	ax_mixed_label["x"]   = r"$\langle x_{%i}\rangle \pm %g\sigma_{x_{%i}}\ \left[\sqrt{\frac{\hbar}{m\omega_0}}\right]$"\
	                        %(site+1, num_sigmas, site+1)
	
	ax_mean_label ["p"]   = r"$\langle p_{%i}\rangle\ \left[\sqrt{\hbar m\omega_0}\right]$"%(site+1)
	ax_std_label  ["p"]   = r"$\sigma_{p_{%i}}\ \left[\sqrt{\hbar m\omega_0}\right]$"%(site+1)
	ax_mixed_label["p"]   = r"$\langle p_{%i}\rangle \pm %g\sigma_{p_{%i}}\ \left[\sqrt{\hbar m\omega_0}\right]$"\
	                        %(site+1, num_sigmas, site+1)
	
	etot_label = r"$E_{tot}\ [\hbar\omega_0]$"
	etot_mean_label = r"$\langle E_{tot}\rangle\ [\hbar\omega_0]$"
else:
	ax_mean_label ["nel"] = r"$\langle \hat{n}_{el,%i}\rangle_r$"%(site+1)
	ax_std_label  ["nel"] = r"$\sigma_{\langle \hat n_{el,%i}\rangle}$"%(site+1)
	ax_mixed_label["nel"] = r"$\langle\hat n_{el,%i}\rangle_r \pm %g\sigma_{\langle\hat n_{el,%i}\rangle}$"%(site+1, num_sigmas, site+1)
	
	ax_mean_label ["nel_adia"] = r"$\langle\hat n_{BO,%i}\rangle_r$"%(BO_state)
	ax_std_label  ["nel_adia"] = r"$\sigma_{\langle\hat n_{BO,%i}\rangle}$"%(BO_state)
	ax_mixed_label["nel_adia"] = r"$\langle\hat n_{BO,%i}\rangle_r \pm %g\sigma_{\langle\hat n_{BO,%i}\rangle}$"%(BO_state, num_sigmas, BO_state)
	
	if use_phonon_sum:
		ax_mean_label ["nph"] = r"$\langle\hat N_{ph}/N\rangle_r$"
		ax_std_label  ["nph"] = r"$\sigma_{\langle\hat N_{ph}/N\rangle}$"
		ax_mixed_label["nph"] = r"$\langle\hat N_{ph}/N\rangle_r \pm %g\sigma_{\langle\hat N_{ph}/N\rangle}$"%(num_sigmas)
	else:
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
	
	etot_label = r"$E_{tot}/\hbar\omega_0$"
	etot_mean_label = r"$\langle E_{tot}\rangle_r/\hbar\omega_0$"

ax_mean_label ["OCDW"]   = r"$\langle\mathcal O_{\rm CDW}\rangle_r$"
ax_std_label  ["OCDW"]   = r"$\sigma_{\mathcal O_{\rm CDW}}$"
ax_mixed_label["OCDW"]   = r"$\langle\mathcal O_{\rm CDW}\rangle_r \pm %g\sigma_{\mathcal O_{\rm CDW}}$"%(num_sigmas)

ax_mean_label ["e_kin"]   = r"$\langle E_{kin}\rangle_r/\hbar\omega_0$"
ax_std_label  ["e_kin"]   = r"$\langle E_{kin}\rangle/\hbar\omega_0$"
ax_mixed_label["e_kin"]   = r"$\left(\langle E_{kin}\rangle_r\pm %g\sigma_{\langle E_{kin}\rangle}\right)/\hbar\omega_0$"\
                            %(num_sigmas)

if use_e_inter_sum:
	ax_mean_label ["e_inter"] = r"$\langle E_{el-ph}\rangle_r/\hbar\omega_0$"
	ax_std_label  ["e_inter"] = r"$\langle E_{el-ph}\rangle/\hbar\omega_0$"
	ax_mixed_label["e_inter"] = r"$\left(\langle E_{el-ph}\rangle_r\pm %g\sigma_{\langle E_{el-ph}\rangle}\right)/\hbar\omega_0$"\
	                            %(num_sigmas)
else:
	ax_mean_label ["e_inter"] = r"$\langle E_{el-ph, %i}\rangle_r/\hbar\omega_0$"%(site+1)
	ax_std_label  ["e_inter"] = r"$\langle E_{el-ph, %i}\rangle/\hbar\omega_0$"%(site+1)
	ax_mixed_label["e_inter"] = r"$\left(\langle E_{el-ph, %i}\rangle_r\pm %g\sigma_{\langle E_{el-ph, %i}\rangle}\right)/\hbar\omega_0$"\
	                            %(site+1, num_sigmas, site+1)


if load_mode:
	if omega != 1 or hopping != 2.5 or gamma != 2.5 or psi_initial != "local" or num_sites_x != 2:
		raise ValueError("Michaels data are only for dimer with omega=1, hopping=gamma=2.5")
	
	sim_name = version
	if version == "m_MTE":
		path_version = "Ehrenfest"
	elif version == "m_SH":
		path_version = "Surfhop"
	elif version == "m_SHDC":
		path_version = "Surfhop_DC"
	
	shared_path = save_path + "michael_" + path_version + "/"
	os.system("mkdir -p " + shared_path)
	
	mean_data = {}
	mean_data["nel"] = np.zeros((len(num_config_list), num_times))
	mean_data["nel_adia"] = np.zeros((len(num_config_list), num_times))
	mean_data["x"] = np.zeros((len(num_config_list), num_times))
	
	std_data = {}
	std_data["nel"] = np.zeros((len(num_config_list), num_times))
	std_data["nel_adia"] = np.zeros((len(num_config_list), num_times))
	std_data["x"] = np.zeros((len(num_config_list), num_times))
	
	for c in range(len(num_config_list)):
		num_configs = num_config_list[c]
		load_path = f"ref_data/convergence/{path_version}_convergence_{num_configs}/"
		
		if version == "m_MTE":
			suffix_nel = "_PROB.DAT"
			suffix_nel_adia = "_PROB.DAT"
		else:
			suffix_nel = "_MIXED.DAT"
			suffix_nel_adia = "_DIST.DAT"
			
		#load nel and time
		data = np.loadtxt(load_path + "ENSEMBLE_CALC_DENSITIES_DIA" + suffix_nel)
		time = data[::show_dist, 0]
		mean_data["nel"][c, :] = data[::show_dist, 1+site]
		std_data ["nel"][c, :] = data[::show_dist, 1+site+num_sites_x]
		
		#load nel_adia
		data = np.loadtxt(load_path + "ENSEMBLE_CALC_DENSITIES_ADIA" + suffix_nel_adia)
		time = data[::show_dist, 0]
		mean_data["nel_adia"][c, :] = data[::show_dist, 1+BO_state]
		std_data ["nel_adia"][c, :] = data[::show_dist, 1+BO_state+num_sites_x]
		
		# load displacement x
		data = np.loadtxt(load_path + "ENSEMBLE_CALC_X.DAT")
		mean_data["x"][c, :] = data[::show_dist, 1+site]
		std_data ["x"][c, :] = data[::show_dist, 1+site+num_sites_x]
else:
	mean_data = {}
	std_data = {}
	for obs in obs_list:
		mean_data[obs] = np.zeros((len(num_config_list), num_times))
		std_data[obs]  = np.zeros((len(num_config_list), num_times))
	
	mean_etot = np.zeros((len(num_config_list), num_times))
	mean_eph  = np.zeros((len(num_config_list), num_times))
	mean_ekin = np.zeros((len(num_config_list), num_times))
	mean_eint = np.zeros((len(num_config_list), num_times))
	
	
	rms_etot = np.zeros((len(num_config_list), num_times))
	std_etot = np.zeros((len(num_config_list), num_times))
	std_eph  = np.zeros((len(num_config_list), num_times))
	std_ekin = np.zeros((len(num_config_list), num_times))
	std_eint = np.zeros((len(num_config_list), num_times))
	
	
	shared_path = ""
	
	for c in range(len(num_config_list)):
		num_configs = num_config_list[c]
		print(f"num_traj/{num_configs=}")
		
		if version in [1,2]:
			conf_title = short_title + f", $N_c$ = {num_configs}"
		else:
			conf_title = short_title + ", $N_{traj}$ = " + str(num_configs)
		
		fig_paper = {}
		axs_paper = {}
		line_paper_obs = {}
		line_paper_e   = {}
		if show_plots or save_plots_single_runs:
			for obs in obs_list:
				num_plots = 1 + int(plot_total_energy)
				if plot_for_paper:
					fig, axs = plt.subplots(num_plots, 1, figsize=(figslim_width, figslim_height[num_plots]))
					fig_paper[obs] = np.array([fig])
				else:
					if subplots_in_thesis:
						fig_list = []
						axs = []
						for i in range(num_plots):
							fig, ax = plt.subplots(1, 1, figsize=(figwide_width[num_plots]/num_plots, figwide_height))
							fig_list.append(fig)
							axs.append(ax)
						fig_paper[obs] = fig_list
					else:
						fig, axs = plt.subplots(1, num_plots, figsize=(figwide_width[num_plots], figwide_height))
						fig_paper[obs] = np.array([fig])
						axs = np.array(axs).reshape(-1)
				
				axs_paper[obs] = axs
				if not subplots_in_thesis:
					fig_paper[obs][0].suptitle(conf_title)
		
		fig, axs = plt.subplots(2, 3, figsize=(16, 9))
		ax = {}
		if "nel" in obs_list:
			ax["nel"]  = axs[0, 0]
		else:
			ax["OCDW"] = axs[0, 0]
		
		ax["nph"] = axs[0, 1]
		ax["x"]   = axs[1, 0]
		ax["p"]   = axs[1, 1]
		ax_energy = axs[0, 2]
		ax_etot = axs[1, 2]
		
		
		ax_energy.set_title("energy")
		ax_etot.set_title("etot")
		
		all_data = {}
		for obs in ax.keys():
			ax[obs].set_title(obs)
		for obs in obs_list:
			all_data[obs] = np.zeros((num_repetitions[c], num_times))
		
		all_etot = np.zeros((num_repetitions[c], num_times))
		all_eph  = np.zeros((num_repetitions[c], num_times))
		all_ekin = np.zeros((num_repetitions[c], num_times))
		all_eint = np.zeros((num_repetitions[c], num_times))
		
		config_path = ""
		
		for i in range(num_repetitions[c]):
			sim = simulate(num_configs, (1+split_simulation_in_half)*i+ (odd_indices and split_simulation_in_half))
			if i == 0:
				full_sim_name = sim.get_sim_name()
				if shared_path == "":
					shared_path = save_path
					cropped_sim_name = full_sim_name[:full_sim_name.rfind("=")+1] + config_list_str
					shared_path += cropped_sim_name.replace(" ", "") + f"_num_repetitions={num_repetitions[-1]}/"
					shared_path += extra_subfolder
					shared_path_sh = shared_path.replace("(", "\(").replace(")", "\)")
					print(f"{shared_path=}")
				if config_path == "":
					config_path = shared_path + str(num_configs) + "/"
					config_path_sh= config_path.replace("(", "\(").replace(")", "\)")
					if save_plots_std_nc_dependence or save_plots_compare:
						os.system("mkdir -p " + shared_path_sh)
					if save_plots_single_runs:
						os.system("mkdir -p " + config_path_sh)
				fig.suptitle(full_sim_name.replace("_", " "))
				
				time = sim.time()[::show_dist]
			
			if "nel" in obs_list:
				all_data["nel"] [i] = sim.nel(site)[::show_dist]
			if "nph" in obs_list:
				if use_phonon_sum:
					all_data["nph"] [i] = sim.nph()[::show_dist].sum(axis=1)/sim.num_electrons
				else:
					all_data["nph"] [i] = sim.nph(site)[::show_dist]
			if "x" in obs_list:
				all_data["x"]   [i] = sim.x(site)[::show_dist]
			if "p" in obs_list:
				all_data["p"]   [i] = sim.p(site)[::show_dist]
			if "OCDW" in obs_list:
				all_data["OCDW"][i] = sim.O_CDW()[::show_dist]
			if "e_kin" in obs_list:
				all_data["e_kin"][i] = sim.energy_kinetic()[::show_dist]
			if "e_inter" in obs_list:
				all_data["e_inter"][i] = sim.energy_interaction("all")[::show_dist, 0, site]
			
			all_etot[i] = sim.energy_total()[::show_dist]
			all_eph[i] = sim.omega*sim.nph()[::show_dist].sum(axis=1)
			all_eint[i] = sim.energy_interaction()[::show_dist]
			all_ekin[i] = sim.energy_kinetic()[::show_dist]
			
			if show_plots or save_plots_single_runs:
				ax_energy.set_prop_cycle(None)
				
				ax_energy.plot(time, all_etot[i], alpha=alpha)
				ax_etot.plot(time, all_etot[i], color='tab:blue', alpha=alpha)
				
				ax_energy.plot(time, all_eph[i], alpha=alpha)
				
				ax_energy.plot(time, all_eint[i], alpha=alpha)
				
				ax_energy.plot(time, all_ekin[i], alpha=alpha)
				
				
				for obs in ax.keys():
					ax[obs].plot(time, all_data[obs][i], color='tab:blue', alpha=alpha)
				
				for obs in obs_list:
					line_paper_obs[obs] = axs_paper[obs][0].plot(time, all_data[obs][i], linestyle='-', color='tab:blue', alpha=alpha)
					if plot_total_energy:
						line_paper_e  [obs] = axs_paper[obs][1].plot(time, all_etot[i], linestyle='-', color='tab:blue', alpha=alpha)
			
			print(f"{i+1}/{num_repetitions[c]}")
		
		if ref is not None:
			if "nel" in obs_list:
				ax["nel"].plot(ref.time(), ref.nel(site), 'k', label="exact")
			else:
				ax["OCDW"].plot(ref.time(), ref.O_CDW(), 'k', label="exact")
			ax["nph"].plot(ref.time(), ref.nph(site), 'k')
			ax["x"].plot(ref.time(), ref.x(site), 'k')
			ax["p"].plot(ref.time(), ref.p(site), 'k')
		
		# overview plot of all observables except energy
		for obs in obs_list:
			mean_data[obs][c]  = np.mean(all_data[obs], axis=0)
			std_data[obs][c]  = np.std(all_data[obs], axis=0, ddof=1)
		for obs in ax.keys():
			ax[obs].plot(time, mean_data[obs][c], color='tab:orange', label="mean")
		ax["nph"].legend()
		ax["x"].set_xlabel(time_label)
		
		# energy
		mean_etot[c] = np.mean(all_etot, axis=0)
		mean_eph[c]  = np.mean(all_eph, axis=0)
		mean_ekin[c] = np.mean(all_ekin, axis=0)
		mean_eint[c] = np.mean(all_eint, axis=0)
		
		rms_etot[c] = np.sqrt(np.mean((all_etot-mean_etot[c, 0])**2, axis=0))
		std_etot[c] = np.std(all_etot, axis=0, ddof=1)
		std_eph[c]  = np.std(all_eph, axis=0, ddof=1)
		std_ekin[c] = np.std(all_ekin, axis=0, ddof=1)
		std_eint[c] = np.std(all_eint, axis=0, ddof=1)
		
		
		ax_energy.set_prop_cycle(None)
		ax_energy.plot(time, mean_etot[c], label="$E_{tot}$")
		ax_energy.plot(time, mean_eph[c], label="$E_{ph}$")
		ax_energy.plot(time, mean_eint[c], label="$E_{inter}$")
		ax_energy.plot(time, mean_ekin[c], label="$E_{kin}$")
		ax_energy.plot(time, np.full(time.size, mean_etot[c, 0]), 'k')
		
		ax_energy.legend()
		
		ax_etot.plot(time, mean_etot[c], color="tab:orange", label="mean")
		if plot_rms_etot:
			ax_etot.plot(time, mean_etot[c, 0] + rms_etot[c], color='tab:red', label="rms + mean(t=0)")
		ax_etot.plot(time, np.full(time.size, mean_etot[c, 0]), 'k')
		ax_etot.legend()
		
		fig.tight_layout()
		if save_plots_single_runs:
			fig.savefig(config_path + "overview.pdf")
		
		if save_plots_single_runs or show_plots:
			for obs in obs_list:
				line_paper_obs[obs][0].set_label("single")
				if plot_total_energy:
					line_paper_e[obs][0].set_label("single")
				
				axs_paper[obs][0].plot(time, mean_data[obs][c], color='tab:orange', label="mean")
				axs_paper[obs][0].set_ylabel(ax_mean_label[obs])
				
				axs_paper[obs][0].set_xlabel(time_label)
				
				axs_paper[obs][0].tick_params(**tick_params)
				if ylabel_coord is not None:
					axs_paper[obs][0].get_yaxis().set_label_coords(*ylabel_coord)
				
				if plot_total_energy:
					axs_paper[obs][1].plot(time, mean_etot[c], color='tab:orange', label="mean")
					if plot_rms_etot:
						axs_paper[obs][1].plot(time, mean_etot[c, 0] + rms_etot[c], color='tab:red', label="rms + mean(t=0)")
					
					axs_paper[obs][1].set_xlabel(time_label)
					axs_paper[obs][1].set_ylabel(etot_label)
					axs_paper[obs][1].legend()
					
					axs_paper[obs][1].tick_params(**tick_params)
					if ylabel_coord is not None:
						axs_paper[obs][1].get_yaxis().set_label_coords(*ylabel_coord)
					
					if not subplots_in_thesis:
						axs_paper[obs][0].annotate('a', annotate_pos, xycoords='axes fraction')
						axs_paper[obs][1].annotate('b', annotate_pos, xycoords='axes fraction')
				else:
					axs_paper[obs][0].legend()
				
				for i in range(len(fig_paper[obs])):
					fig = fig_paper[obs][i]
					fig.tight_layout(pad=tight_layout_pad)
					if save_plots_single_runs:
						filename = config_path
						if plot_total_energy and i+1 == len(fig_paper[obs]):
							filename += f"single_runs_energy"
						else:
							filename += f"single_runs_{obs}"
							if len(fig_paper[obs]) != 1:
								filename += f"_ax{i}_{len(fig_paper[obs])}"
						filename += ".pdf"
						if save_plots_single_runs:
							fig.savefig(filename)
						if not show_plots:
							plt.close(fig)
			
			if not show_plots:
				fig_paper[obs] = None
	
	sim_name = sim.get_sim_name()
	sim_name = sim_name[:sim_name.rfind("_")]


# overview
fig_conf, axs_conf = plt.subplots(2, 4, figsize=(16, 9))
fig_conf.suptitle(sim_name.replace("_", " "))
ax = {}
if "nel" in obs_list:
	ax["nel"] = axs_conf[0, 0]
else:
	ax["OCDW"] = axs_conf[0, 0]
ax["nph"] = axs_conf[0, 1]
ax["x"]   = axs_conf[1, 0]
ax["p"]   = axs_conf[1, 1]
ax_rms_etot = axs_conf[0, 2]
ax_rms_etot_log = axs_conf[0, 3]
ax_std_etot = axs_conf[1, 2]
ax_std_etot_log = axs_conf[1, 3]

for obs in obs_list:
	try:
		ax[obs].set_title(f"std {obs}")
		ax[obs].plot(time, std_data[obs].T)
		ax[obs].set_ylabel(ax_std_label[obs])
		ax[obs].set_xlabel(time_label)
	except KeyError:
		print(f"no axes for observable \"{obs}\" reserved, dropping overview plot")

ax_rms_etot.set_title("rms etot")
ax_rms_etot_log.set_title("rms etot")
ax_rms_etot_log.set_yscale("log")

ax_std_etot.set_title("std etot")
ax_std_etot_log.set_title("std etot")
ax_std_etot_log.set_yscale("log")

if not load_mode:
	ax_rms_etot.plot(time, rms_etot.T)
	ax_rms_etot_log.plot(time, rms_etot.T)
	ax_std_etot.plot(time, std_etot.T)
	ax_std_etot_log.plot(time, std_etot.T)


config_label_list = [str(num_configs) for num_configs in num_config_list]
if compare_show_exact and ref is not None:
	config_label_list.append(exact_label)

ax["nph"].legend(config_label_list)
fig_conf.tight_layout(pad=tight_layout_pad)

if version in [1, 2]:
	legend_title = "$N_c$ ="
else:
	legend_title = "$N_{traj}$ ="


def plot_etot(ax):
	ax.plot(time, mean_etot.T)
	ax.set_xlabel(time_label)
	ax.set_ylabel(etot_mean_label)
	if compare_show_exact and ref is not None:
		initial_energy = mean_etot[:, 0].mean()
		ax.plot(ref.time(), np.full(ref.time().size, initial_energy),
		        linestyle=exact_line_style, color=exact_color, label = exact_label)
	ax.legend(config_label_list, title=legend_title, loc='upper left', ncol= 1+ (figwide_height < 3 and len(config_label_list) > 4))

# total energy
if plot_total_energy:
	if version in [1, 2]:
		num_plots = 1 + int(plot_total_energy)
		fig_etot, axs = plt.subplots(num_plots, 1, figsize=(figslim_width, figslim_height[num_plots]))
		fig_etot.suptitle(short_title)
		
		plot_etot(axs[0])
		
		ax_etot_rms  = axs[1]
		ax_etot_rms.plot(time, rms_etot.T)
		ax_etot_rms.set_xlabel(time_label)
		ax_etot_rms.set_ylabel(r"rms $E_{tot}\ [\hbar\omega_0]$")
		ax_etot_rms.legend(config_label_list, title=legend_title)
		
		fig_etot.tight_layout()
		if save_plots_compare:
			fig_etot.savefig(shared_path + "compare_etot.pdf")



# specific observable in comparison
if save_plots_compare or show_plots:
	for obs in obs_list:
		# mean and std seperated
		num_plots = 2 + int(include_energy_to_obs_compare)
		num_plots_width = num_plots
		if plot_for_paper:
			fig, axs = plt.subplots(num_plots, 1, figsize=(figslim_width, figslim_height[num_plots]))
			figs = np.array([fig])
		else:
			if subplots_in_thesis:
				fig_list = []
				axs = []
				num_plots_width = num_plots if fixed_figwidth_count is None else fixed_figwidth_count
				for i in range(num_plots):
					fig, ax = plt.subplots(1, 1, figsize=(figwide_width[num_plots_width]/num_plots_width, figwide_height))
					fig_list.append(fig)
					axs.append(ax)
				figs = fig_list
			else:
				fig, axs = plt.subplots(1, num_plots, figsize=(figwide_width[num_plots], figwide_height))
				figs = np.array([fig])
				axs = np.array(axs).reshape(-1)
		
		if not subplots_in_thesis:
			figs[0].suptitle(short_title_rep)
		
		ax_obs_mean = axs[0]
		ax_obs_std  = axs[1]
		if not subplots_in_thesis:
			ax_obs_mean.annotate('a', annotate_pos, xycoords='axes fraction')
			ax_obs_std.annotate('b', annotate_pos, xycoords='axes fraction')
		
		ax_obs_mean.set_ylabel(ax_mean_label[obs])
		ax_obs_std.set_ylabel(ax_std_label[obs])
		mean_obs = mean_data[obs]
		std_obs  = std_data[obs]
		
		ax_obs_mean.plot(time, mean_obs.T)
		if compare_show_exact and ref is not None:
			if obs == "nel":
				ax_obs_mean.plot(ref.time(), ref.nel(site), linestyle=exact_line_style, color=exact_color, label=exact_label)
			elif obs == "nph":
				ax_obs_mean.plot(ref.time(), ref.nph(site), linestyle=exact_line_style, color=exact_color, label=exact_label)
			elif obs == "x":
				ax_obs_mean.plot(ref.time(), ref.x(site), linestyle=exact_line_style, color=exact_color, label=exact_label)
			elif obs == "p":
				ax_obs_mean.plot(ref.time(), ref.p(site), linestyle=exact_line_style, color=exact_color, label=exact_label)
		ax_obs_mean.set_xlabel(time_label)
		ax_obs_mean.tick_params(**tick_params)
		if ylabel_coord is not None:
			ax_obs_mean.get_yaxis().set_label_coords(*ylabel_coord)
		
		ax_obs_std.plot(time, std_obs.T)
		ax_obs_std.set_xlabel(time_label)
		ax_obs_std.tick_params(**tick_params)
		if ylabel_coord is not None:
			ax_obs_std.get_yaxis().set_label_coords(*ylabel_coord)
		
		if allow_crop_nel_range and obs == "nel":
			crop_nel_time_index = np.argmax(time >= crop_nel_start_time)
			#print(f"{crop_nel_time_index=}")
			crop_nel = not np.any(mean_obs[:, crop_nel_time_index:] > crop_nel_threshold)
			print(crop_nel)
		else:
			crop_nel = False
		
		if crop_nel:
			crop_nel_ymax_unscaled = mean_obs[:, crop_nel_time_index:].max()
			ax_obs_mean.set_ylim(ymax=crop_nel_ymax_unscaled*1.1)
		
		if include_energy_to_obs_compare:
			plot_etot(axs[2])
			if not subplots_in_thesis:
				axs[2].annotate('c', annotate_pos, xycoords='axes fraction')
			axs[2].tick_params(**tick_params)
			if ylabel_coord is not None:
				axs[2].get_yaxis().set_label_coords(*ylabel_coord)
		else:
			if legends_in_figs == "all" or obs in legends_in_figs:
				ax_obs_mean.legend(config_label_list, title=legend_title)
				#ax_obs_std.legend(config_label_list, title=legend_title)
		
		if save_plots_compare or show_plots:
			for i in range(len(figs)):
				fig = figs[i]
				fig.tight_layout(pad=tight_layout_pad)
				if save_plots_compare:
					filename = shared_path
					if include_energy_to_obs_compare and i + 1 == num_plots:
						filename += f"compare_energy"
					else:
						filename += f"compare_{obs}"
						if obs not in ["e_kin"]:
							filename += (f"{site+1}" if obs != "nph" or not use_phonon_sum else "Sum")
						if len(figs) != 1:
							filename += f"_ax{i}"
					filename += f"_numfigs{num_plots_width}.pdf"
					fig.savefig(filename)
				if not show_plots:
					plt.close(fig)
		
		if not show_plots:
			figs = None
		
		# mean and std combined into one stacked plot
		num_plots = 1 + int(include_energy_to_obs_compare)
		num_plots_width = num_plots
		if plot_for_paper:
			fig, axs = plt.subplots(num_plots, 1, figsize=(figslim_width, figslim_height[num_plots]))
			figs = np.array([fig])
		else:
			if subplots_in_thesis:
				fig_list = []
				axs = []
				num_plots_width = num_plots if fixed_figwidth_count is None else fixed_figwidth_count
				for i in range(num_plots):
					fig, ax = plt.subplots(1, 1, figsize=(figwide_width[num_plots_width]/num_plots_width, figwide_height))
					fig_list.append(fig)
					axs.append(ax)
				figs = fig_list
			else:
				fig, axs = plt.subplots(1, num_plots, figsize=(figwide_width[num_plots], figwide_height))
				figs = np.array([fig])
				fig_shaded_obs = np.array(figs).reshape(-1)
		
		if not subplots_in_thesis:
			figs[0].suptitle(short_title_rep)
		
		
		ax_shaded_obs = axs[0]
		ax_shaded_obs.tick_params(**tick_params)
		
		if auto_shift:
			data_range = mean_obs[0].max()-mean_obs[0].min()
			if crop_nel:
				data_range = crop_nel_threshold
			shift_dist = max(0.25*data_range, 2.5*num_sigmas*std_obs.max())
			print(f"{obs}: {shift_dist=}", end=" ")
			mag = np.floor(np.log10(shift_dist))
			normalized_shift = shift_dist*10**(-mag)
			shift_dist = allowed_shift_distances[np.argmin(np.abs(allowed_shift_distances - normalized_shift))]*10**mag
			#shift_dist = np.round(shift_dist*10**(-mag))*10**mag
			print(f"-> {shift_dist}")
		else:
			shift_dist = shift_distances[obs]
		for c in range(len(num_config_list)):
			num_configs = num_config_list[c]
			ax_shaded_obs.fill_between(time,
			                           c*shift_dist + mean_obs[c] - num_sigmas*std_obs[c],
			                           c*shift_dist + mean_obs[c] + num_sigmas*std_obs[c],
			                           alpha=0.5)
			ax_shaded_obs.plot(time, c*shift_dist + mean_obs[c], label=config_label_list[c])
		
		if crop_nel:
			ax_shaded_obs.set_ylim(ymax=(crop_nel_ymax_unscaled+len(num_config_list)*shift_dist)*1.1)
		
		ax_shaded_obs.set_ylabel(ax_mixed_label[obs])
		ax_shaded_obs.set_xlabel(time_label)
		
		if include_energy_to_obs_compare:
			plot_etot(axs[1])
			if not subplots_in_thesis:
				axs[0].annotate('a', annotate_pos, xycoords='axes fraction')
				axs[1].annotate('b', annotate_pos, xycoords='axes fraction')
			axs[1].tick_params(**tick_params)
		else:
			if legends_in_figs == "all" or obs in legends_in_figs:
				ax_shaded_obs.legend(title=legend_title)
		
		for i in range(len(figs)):
			fig = figs[i]
			fig.tight_layout(pad=tight_layout_pad)
			if save_plots_compare:
				filename = shared_path
				if include_energy_to_obs_compare and i + 1 == num_plots:
					filename += f"compare_energy"
				else:
					filename += f"stacked_compare_{obs}"
					if obs not in ["e_kin"]:
						filename += (f"{site+1}" if obs != "nph" or not use_phonon_sum else "Sum")
					if len(figs) != 1:
						filename += f"_ax{i}"
				filename += f"_numfigs{num_plots_width}.pdf"
				fig.savefig(filename)
			if not show_plots:
				plt.close(fig)

if save_plots_compare or save_plots_single_runs or save_plots_std_nc_dependence:
	print("saving to:")
	print(shared_path)

etot_mean_abs_label = r"$|\langle E_{tot}\rangle_r-\langle E_{tot}(\omega_0t=0)\rangle_r|/\hbar\omega_0$"
num_configs_label = r"number of configurations $N_c$"

if save_plots_std_nc_dependence or show_plots:
	# find index at 3rd maximum of phonons, designed for hopping=gamma=2.5, omega=1
	#print(num_config_list)
	#print(std_data["nph"][:, ind])
	#print(mean_etot[:, ind])
	
	color_mean = "tab:blue"
	color_std = "tab:blue"
	color_energy = "tab:orange"
	
	cmap = "Greys"
	
	
	std_nph_label = ax_std_label["nph"]
	diff_mean_nph_label = r"$\Delta \langle N_{ph,%i}\rangle_r$"%(site+1)
	abs_diff_mean_nph_label = "$|" + diff_mean_nph_label[1:-1] + "|$"
	
	init_energy = np.mean(mean_etot[:, 0])
	
	diff_mean_nph = None
	abs_diff_mean_nph = None
	if ref is not None:
		dt_ref = ref.time()[1] - ref.time()[0]
		dt_sim = time[1] - time[0]
		
		from fractions import Fraction
		ratio = Fraction(dt_ref/dt_sim).limit_denominator(ref.time().size)
		stepping_sim = ratio.numerator
		stepping_ref = ratio.denominator
		
		if abs(ref.time()[stepping_ref] - time[stepping_sim]) > 1e-8:
			print(f"\nWARNING: bad time_alignements!\n{stepping_sim=}, {stepping_ref=}, {ref.time()[stepping_ref]=}, {time[stepping_sim]=}\n")
		ref_nph = ref.nph(site)[::stepping_ref]
		sim_nph = mean_data["nph"][:, ::stepping_sim]
		min_len = min(ref_nph.size, sim_nph.shape[1])
		
		diff_mean_nph = sim_nph[:, :min_len] - ref_nph[:min_len].reshape(1, -1)
		abs_diff_mean_nph = np.abs(diff_mean_nph)
	
	for cut_time in cut_times_nph_std_and_energy_vs_nc:
		ind = np.argmax(time >= cut_time)
		#time_arg_label = ("$($" + time_label + "$=%g)$"%(cut_time)).replace("$$", "")
		#etot_mean_time_arg_label = etot_mean_label.replace("/", time_arg_label[1:-1]+"/")
		time_title = (time_label + "$=%g$"%(cut_time)).replace("$$", "")
		
		
		
		plt.figure(figsize=(figwidth_nph_std_and_energy_vs_nc, figwide_height))
		plt.title(time_title)
		ax1 = plt.gca()
		
		show_mean_diff = allow_show_mean_diff and abs_diff_mean_nph is not None and ind//stepping_sim < abs_diff_mean_nph.shape[1]
		
		line_std = ax1.plot(num_config_list,  std_data["nph"][:, ind], '.', label=std_nph_label)
		
		ax1.set_ylabel(std_nph_label, color=color_std)
		
		ax1.tick_params(axis='y', labelcolor=color_std)
		ax1.set_xlabel(num_configs_label)
		ax1.set_ylim(ymin=0)
		
		ax2 = ax1.twinx()
		line_energy = ax2.plot(num_config_list, mean_etot[:, ind], '+',
		                       label=etot_mean_label[:etot_mean_label.rfind("/")] + "$", color=color_energy)
		ax2.tick_params(axis='y', labelcolor=color_energy)
		ax2.set_ylabel(etot_mean_label, color=color_energy)
		if mean_etot[:, ind].min() >= init_energy:
			ax2.set_ylim(ymin=init_energy)
		
		plt.tight_layout(pad=tight_layout_pad)
		
		if save_plots_std_nc_dependence:
			plt.savefig(shared_path + f"max_nph_std_and_energy_vs_nc_wt={cut_time}.pdf")
		
		if show_mean_diff:
			line_mean = ax1.plot(num_config_list, abs_diff_mean_nph[:, ind//stepping_sim], 'x',
			         color=color_mean, label=abs_diff_mean_nph_label)
			ax1.set_ylabel(std_nph_label + ", " + abs_diff_mean_nph_label, color=color_std)
			
			fake_energy = mean_etot[:, ind] - ax2.get_ylim()[0]
			fake_energy /= ax2.get_ylim()[1] - ax2.get_ylim()[0]
			fake_energy *= ax1.get_ylim()[1] - ax1.get_ylim()[0]
			fake_energy += ax1.get_ylim()[0]
			ax1.plot(num_config_list, fake_energy, '.', markersize=0)
			lns = line_mean + line_std + line_energy
			labs = [l.get_label() for l in lns]
			ax1.legend(lns, labs, loc=0)
			
			if save_plots_std_nc_dependence:
				plt.savefig(shared_path + f"max_nph_std_mean_diff_and_energy_vs_nc_wt={cut_time}_legend.pdf")
		
	
	num_plots = 2
	if fixed_figwidth_count is not None:
		num_plots = fixed_figwidth_count
	
	plt.figure(figsize=(figwide_width[num_plots]/num_plots, figwide_height))
	plt.pcolormesh(time, num_config_list, std_data["nph"][:, :], shading="nearest",
	               cmap=cmap, vmin=0)
	plt.xlabel(time_label)
	plt.ylabel(num_configs_label)
	cb = plt.colorbar()
	cb.set_label(ax_std_label["nph"])
	plt.tight_layout(pad=tight_layout_pad)
	
	if save_plots_std_nc_dependence:
		plt.savefig(shared_path + f"color_std_nph{site+1}_vs_nc_time.png", dpi=600)
	
	
	if diff_mean_nph is not None:
		plt.figure(figsize=(figwide_width[num_plots]/num_plots, figwide_height))
		vmax = diff_mean_nph.max()
		vmin = diff_mean_nph.min()
		vrange = max(vmax, -vmin)
		print(f"{vmin=}, {vmax=}")
		
		plt.pcolormesh(ref.time()[::stepping_ref][:min_len], num_config_list, diff_mean_nph, shading="nearest",
		               cmap="bwr", vmin=-vrange, vmax=vrange)
		plt.xlabel(time_label)
		plt.ylabel(num_configs_label)
		cb = plt.colorbar()
		cb.set_label(diff_mean_nph_label)
		plt.tight_layout(pad=tight_layout_pad)
		
		if save_plots_std_nc_dependence:
			plt.savefig(shared_path + f"color_diff_mean_nph{site+1}_vs_nc_time.png", dpi=600)
		
	
	energy_diff = np.abs(mean_etot[:, :]-init_energy)
	
	plt.figure(figsize=(figwide_width[num_plots]/num_plots, figwide_height))
	plt.pcolormesh(time, num_config_list, energy_diff, shading="nearest",
	               cmap=cmap, vmin=0)
	
	if plot_max_energy_curve:
		max_energy_diff_ind = np.argmax(energy_diff, axis=0)
		plt.plot(time, np.array(num_config_list)[max_energy_diff_ind], color="tab:orange")
	
	plt.xlabel(time_label)
	plt.ylabel(num_configs_label)
	cb = plt.colorbar()
	cb.set_label(etot_mean_abs_label)
	plt.tight_layout(pad=tight_layout_pad)
	
	if save_plots_std_nc_dependence:
		plt.savefig(shared_path + "color_energy_vs_nc_time.png", dpi=600)
		
if num_repetitions.min() != num_repetitions.max():
	plt.figure()
	plt.title("num repetitions")
	plt.plot(num_config_list, num_repetitions)
	plt.xlabel(num_configs_label)
	plt.ylabel("num repetitions")
	if save_plots_compare or save_plots_std_nc_dependence:
		plt.savefig(shared_path + "num_repetitions.pdf")

if show_plots:
	plt.show()
