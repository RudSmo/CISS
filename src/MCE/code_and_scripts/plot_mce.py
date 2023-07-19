import numpy as np
import matplotlib.pyplot as plt
import data_manager


use_log_scale = False

def n_el_ph_scaling(y):
	if use_log_scale:
		return np.log10(y)
	return y

def xp_scaling(y):
	sign = np.sign(y)
	return y # sign*(np.abs(y)**(1))

version = 1
allow_config_delete = False
allow_config_clone = True
hopping = .5
omega = 1
gamma = 2.
num_sites_x = 11
num_sites_y = 1

num_configs = 1000
max_num_configs = num_configs
num_reps = 0 # number of simulations which error plots vs energy
nel_cmp_slope = 0. #0.05*4./gamma
nph_cmp_slope = 0. #1.*omega/gamma


dt = 0.01
t_max = 20.
save_step = 1

initial_state = "bare"
psi_initial = "local"

if psi_initial == "local":
	psi_state = (0, (num_sites_x-1)//2) # (y, x)
	phonon_dressed_displaced_sites = None
elif psi_initial == "eigen":
	psi_state = 1
	phonon_dressed_displaced_sites = (0, 0)


calc_BO_electron_density = False
BO_nel_save_dist = 10
load_sampled_state = False
force_new_simulation = False
pancakeScale = 1
pancakeSampling = False
scale_el = None

local_psi_init_select_eigenstate = None
"""
if version == 1:
	local_psi_init_select_eigenstate = False
else:
	local_psi_init_select_eigenstate = True
"""

save_data = False
show_plots = True
save_plots = False

#limits of color plots of nel and nph
if use_log_scale:
	vmin_nel = -6
	vmax_nel = 0
	vmin_nph = -5
	vmax_nph = 1.2
else:
	vmin_nel = 0
	vmax_nel = .5
	vmin_nph = 0
	vmax_nph = None

figsize = (8, 4)
title = initial_state + " " + psi_initial + (str(psi_state) if psi_initial == "eigen" else "")\
        + ": $t_0=%g, \\omega=%g$, $\\gamma=%g,$ conf %i$\\rightarrow$%i"%(hopping, omega, gamma, num_configs, max_num_configs)

ref_path = "wigner_functions/%s_%s_omega=%g_gamma=%g/"\
           %(initial_state, "high" if psi_initial == "eigen" and psi_state == 1 else psi_initial,\
             omega, gamma)


try:
	ref_high = np.loadtxt(ref_path + "PROBABILITIES.DAT")
	ref_time = ref_high[:, 0]/hopping
	ref_nel  = ref_high[:, 4:2:-1]
	ref_q = -np.load(ref_path + "q1.npy")
	ref_pq = -np.load(ref_path + "pq1.npy")
except Exception:# as e:
	#print("no Wigner reference Data found:")
	#print(repr(e))
	ref_high = None

try:
	ref = data_manager.RefData(initial_state, psi_initial, psi_state, num_sites_x, num_sites_y, hopping, omega, gamma)
	ref.set_cut_time(t_max)
except Exception as e:
	print("no reference Data found:")
	print(repr(e))
	ref = None


num_times = int(np.round(t_max/dt)) + 1  # +1 for t=0

sim = data_manager.Simulation_MCE(version, hopping, omega, gamma,
                                  num_sites_x, num_sites_y,
                                  t_max, dt,
                                  num_configs,
                                  initial_state, psi_initial, psi_state,
                                  max_num_configs,
                                  allow_config_delete=allow_config_delete,
                                  allow_config_clone=allow_config_clone,
                                  printProgress=True,
                                  phonon_dressed_displaced_sites=phonon_dressed_displaced_sites,
                                  force_new_simulation=force_new_simulation,
                                  load_sampled_state=load_sampled_state,
                                  store_results=save_data,
                                  extra_subfolder="0",
                                  scale_el=scale_el, pancakeSampling=pancakeSampling, pancakeScale=pancakeScale,
                                  calc_BO_electron_density=calc_BO_electron_density,
                                  BO_nel_save_dist=BO_nel_save_dist,
                                  local_psi_init_select_eigenstate=local_psi_init_select_eigenstate)


if num_sites_x == 3:
	plt.figure()
	plt.title("x sim")
	plt.plot(sim.time(), sim.x(1) - .5*(sim.x(0)+sim.x(2)))

if max_num_configs != num_configs and version == 2:
	plt.figure()
	plt.title("num configs")
	plt.plot(sim.time(), sim.get_num_configs_over_time(), '.-')

plt.figure()
plt.title(r"$nph-(\langle x\rangle^2 + \langle p\rangle^2)/2$")
var_ph = sim.nph(squeezeY=False).sum(axis=(1, 2)) - (np.sum(sim.x(squeezeY=False)**2, axis=(1, 2)) + np.sum(sim.p(squeezeY=False)**2, axis=(1, 2)))/2
plt.plot(sim.time(), var_ph)
if ref is not None:
	var_ph_ref = np.sum(ref.nph(), axis=1) - np.sum(ref.x()**2+ref.p()**2, axis=1)/2
	plt.plot(ref.time(), var_ph_ref, '--')

plt.figure()
plt.title("x, p phasespace " + title)
plt.plot(sim.x(), sim.p())
if ref is not None:
	plt.plot(ref.x(), ref.p(), '--')

plt.figure()
plt.title("x " + title)
plt.plot(sim.time(), sim.x())
if ref is not None:
	plt.plot(ref.time(), ref.x(), '--')

plt.figure()
plt.title("p " + title)
plt.plot(sim.time(), sim.p())
if ref is not None:
	plt.plot(ref.time(), ref.p(), '--')


fig_nph = plt.figure(figsize=figsize)
#plt.title("nph " + title)
plt.plot(sim.time(), sim.nph())
if ref is not None:
	plt.plot(ref.time(), ref.nph(), '--')
plt.ylabel("$n_{ph}$")
plt.xlim((0, t_max))
if use_log_scale:
	plt.yscale("log")

fig, axs = plt.subplots(1 + (ref is not None), 1, squeeze=False, sharex="all", sharey="all")
fig.suptitle("nel " + title)
cm = axs[0,0].pcolormesh(sim.time(), np.arange(num_sites_x), n_el_ph_scaling(sim.nel().T),
                         shading="nearest", vmin=vmin_nel, vmax=vmax_nel)
if ref is not None:
	axs[0,0].set_title("MCE")
	axs[1,0].set_title("ED/DMRG")
	axs[1,0].pcolormesh(ref.time(), np.arange(num_sites_x), n_el_ph_scaling(ref.nel().T),
	                    shading="nearest", vmin=vmin_nel, vmax=vmax_nel)
#plt.colorbar(cm)


fig, axs = plt.subplots(1 + (ref is not None), 1, squeeze=False, sharex="all", sharey="all")
fig.suptitle("nph " + title)
cm = axs[0,0].pcolormesh(sim.time(), np.arange(num_sites_x), n_el_ph_scaling(sim.nph().T),
                         shading="nearest", vmin=vmin_nph, vmax=vmax_nph)
if ref is not None:
	axs[0,0].set_title("MCE")
	axs[1,0].set_title("ED/DMRG")
	axs[1,0].pcolormesh(ref.time(), np.arange(num_sites_x), n_el_ph_scaling(ref.nph().T),
	                    shading="nearest", vmin=vmin_nph, vmax=vmax_nph)
#plt.colorbar(cm)


plt.figure()
scaled_x = xp_scaling(sim.x().T)
vlim = max(scaled_x.max(), -scaled_x.min())
plt.title("x " + title)
plt.pcolormesh(sim.time(), np.arange(num_sites_x), scaled_x,
               shading="nearest", cmap="RdBu", vmin=-vlim, vmax=vlim)
plt.colorbar()

plt.figure()
scaled_p = xp_scaling(sim.p().T)
vlim = max(scaled_p.max(), -scaled_p.min())
plt.title("p " + title)
plt.pcolormesh(sim.time(), np.arange(num_sites_x), scaled_p,
               shading="nearest", cmap="RdBu", vmin=-vlim, vmax=vlim)
plt.colorbar()


fig_nel = plt.figure(figsize=figsize)
#plt.title("nel " + title)
plt.plot(sim.time(), sim.nel())
if ref is not None:
	plt.plot(ref.time(), ref.nel(), '--')
if ref_high is not None:
	plt.plot(ref_time, ref_nel, "--")
plt.xlim((0, t_max))
plt.ylabel("$n_{el}$, conf %i$\\rightarrow$%i"%(num_configs, max_num_configs))
if use_log_scale:
	plt.yscale("log")

if calc_BO_electron_density:
	print("min(sum BO_nel)=", sim.nel_BO().sum(axis=1).min())
	print("max(sum BO_nel)=", sim.nel_BO().sum(axis=1).max())
	plt.figure()
	plt.title("BO nel " + title)
	plt.plot(sim.time(BO_steps=True), sim.nel_BO())
	plt.xlabel("t")
	plt.ylabel("BO nel")
	if num_sites_y == 1 and num_sites_x == 2:
		try:
			if psi_initial == "eigen":
				if psi_state == 0:
					psi_initial_old = "low"
				if psi_state == -1 or psi_state == num_sites_x-1:
					psi_initial_old = "high"
			else:
				psi_initial_old = psi_initial
			wigner_ref = data_manager.Wigner(omega/hopping, gamma/hopping, initial_state, psi_initial_old)
			plt.plot(wigner_ref.time()*hopping, wigner_ref.nel_BO(), '--')
		except Exception as e:
			print("BO occupation number from Wigner failed:", e)
	

plt.figure()
plt.title("kinetic energy " + title)
plt.plot(sim.time(), sim.energy_kinetic())
if ref is not None:
	plt.plot(ref.time(), ref.energy_kinetic(), '--')


plt.figure()
plt.title("electron phonon interaction per site")
plt.plot(sim.time(), sim.energy_interaction(site="all")[:, 0, :])


fig_energy = plt.figure(figsize=(8,4))
e_ph = omega * sim.nph(squeezeY=False).sum(axis=(1,2))
e_kin = sim.energy_kinetic()
e_inter = sim.energy_interaction()
e_tot = e_ph + e_kin + e_inter

plt.plot(sim.time(), np.full(sim.time().size, e_tot[0]), 'k')
plt.plot(sim.time(), e_ph, label="$E_{ph}$")
plt.plot(sim.time(), e_kin, label="$E_{kin}$")
plt.plot(sim.time(), e_inter, label="$E_{inter}$")
plt.plot(sim.time(), e_tot, label="$E_{tot}$")
plt.xlim((0, t_max))
plt.ylabel("energy")
#if ref is not None and dt == ref.time()[1]:
#	plt.plot(ref.time(), omega*(np.sum(obs_nph, axis=(1, 2))-np.sum(ref.nph(), axis=1)), '--', label="omega*(nph -nph_ref)")
#	plt.plot(time, (var_ph-var_ph_ref)*omega, '--', label="omega*(diff var nph)")
plt.legend()


# FIXME old syntax, doesn't work like that
if False:
	fig, axs = plt.subplots(num_sites_y, num_sites_x, figsize=(16,9), squeeze=False)
	fig.suptitle("weights end, " + title)
	fig.tight_layout()
	for j in range(num_sites_y):
		for i in range(num_sites_x):
			cm = axs[j, i].scatter(np.sqrt(2)*hoz[j, i, :num_configs].real, np.sqrt(2)*hoz[j, i, :num_configs].imag, c=np.abs(el_coeff[0, j, i, :num_configs]), vmin=0)
			axs[j, i].set_xlabel("x")
			axs[j, i].set_ylabel("p")
			#fig.colorbar(cm, ax=axs[i])

	plt.figure()
	plt.plot(np.abs(el_coeff[0, 0, :,]).T, ".")

if num_sites_x == 2 and num_sites_y == 1:
	plt.figure()
	plt.title("q")
	plt.plot(sim.time(), sim.q())
	if ref is not None:
		plt.plot(ref.time(), ref.q(), '--')
	if ref_high is not None:
		plt.plot(ref_time, ref_q, "--")
	plt.xlim((0, t_max))
	
	plt.figure()
	plt.title("pq")
	plt.plot(sim.time(), sim.pq())
	if ref is not None:
		plt.plot(ref.time(), ref.pq(), '--')
	if ref_high is not None:
		plt.plot(ref_time, ref_pq, "--")
	plt.xlim((0, t_max))


if save_plots:
	fig_nel.tight_layout()
	fig_nph.tight_layout()
	fig_e_tot.tight_layout()
	
	
	fig_nel.savefig(path+"nel.png", dpi=300)
	fig_nph.savefig(path+"nph.png", dpi=300)
	fig_e_tot.savefig(path+"e_tot.png", dpi=300)


# FIXME still old syntax
if ref is not None and num_reps != 0:
	from scipy.interpolate import interp1d
	
	fig, axs = plt.subplots(3, 2, figsize=(16, 9))
	axs[0, 0].set_title(r"$\Delta n_{el, 0}$")
	axs[0, 0].set_xlabel(r"$\Delta E_{tot}$")
	axs[1, 0].set_xlabel(r"$\max_t \Delta E_{tot}$")
	axs[2, 0].set_xlabel(r"$t\cdot t_0$")
	axs[0, 1].set_title(r"$\Delta n_{ph, 0}$")
	axs[0, 1].set_xlabel(r"$\Delta E_{tot}$")
	axs[1, 1].set_xlabel(r"$\max_t \Delta E_{tot}$")
	axs[2, 1].set_xlabel(r"$t\cdot t_0$")
	
	nel0_mean = obs_nel.copy()
	
	for rep in range(0, num_reps):
		
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
		
		nel_ref_interp = interp1d(ref.time(), ref.nel(), axis=0, kind="cubic")
		nel_diff = nel_ref_interp(time) - obs_nel[:, 0, :]
		axs[0, 0].plot(e_tot, nel_diff)
		axs[1, 0].plot(e_tot_maxed, nel_diff)
		axs[2, 0].plot(time, nel_diff)
		
		nph_ref_interp = interp1d(ref.time(), ref.nph(), axis=0, kind="cubic")
		nph_diff = nph_ref_interp(time) - obs_nph[:, 0, :]
		axs[0, 1].plot(e_tot, nph_diff)
		axs[1, 1].plot(e_tot_maxed, nph_diff)
		axs[2, 1].plot(time, nph_diff)
		
		if rep != num_reps-1:
			hoz, el_coeff = initial_state_MCE(num_sites_x, num_sites_y, num_configs, max_num_configs,
			                                  initial_state, psi_initial,
			                                  hopping, omega, gamma, psi_state=psi_state)
			
			# simulation
			obs_y, obs_nph, obs_nel, obs_el_ekin, obs_el_ph_inter = simulate_system_MCE(hoz, el_coeff, hopping, omega, gamma, num_times, dt, num_configs, printProgress=True)
			nel0_mean += obs_nel
	
	
	axs[0, 0].plot(e_tot, nel_cmp_slope*e_tot, 'k')
	axs[0, 0].plot(e_tot, -nel_cmp_slope*e_tot, 'k')
	axs[1, 0].plot(e_tot_maxed, nel_cmp_slope*e_tot_maxed, 'k')
	axs[1, 0].plot(e_tot_maxed, -nel_cmp_slope*e_tot_maxed, 'k')
	
	axs[0, 1].plot(e_tot, nph_cmp_slope*e_tot, 'k')
	axs[0, 1].plot(e_tot, -nph_cmp_slope*e_tot, 'k')
	axs[1, 1].plot(e_tot_maxed, nph_cmp_slope*e_tot_maxed, 'k')
	axs[1, 1].plot(e_tot_maxed, -nph_cmp_slope*e_tot_maxed, 'k')
	
	fig.tight_layout()
	
	nel0_mean /= num_reps
	plt.figure()
	plt.title("nel0_mean")
	plt.plot(time, nel0_mean[:, 0])
	plt.plot(ref.time(), ref.nel(), '--')

if show_plots or not save_data:
	plt.show()

