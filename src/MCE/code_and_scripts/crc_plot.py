import numpy as np
import matplotlib.pyplot as plt


import data_manager

initial_state = "BCDW"
psi_initial="local"
hopping = 1.
omega_good = .1
gamma_good = .4

omega_bad = 2
gamma_bad = 4

num_sites = 13
num_electrons = num_sites//2
cubic_interpolate = False	# just special case num_sites = 13, num_electrons=6, omega=2, gamma=4

# number of ehrenfest trajectories (randomly chosen)
# use num_samples = 0 for single trajectory starting in rest on all sites
num_samples = 4000

if num_sites == 2:
	t_max = 10.
else:
	t_max = 6.5

mce_good_spacing = 5

if num_sites == 2:
	mce_path_good = "data_export/BCDW_psi_initial=local_omega=%g_gamma=%g_mce_numconf=30/"%(omega_good, gamma_good)
	mce_path_bad = "data_export/BCDW_psi_initial=local_omega=%g_gamma=%g_mce_numconf=40/"%(omega_bad, gamma_bad)
	
	mce_bad_spacing = mce_good_spacing
	mce_style = "x"
	
	sh_good = np.loadtxt("data_export/sh_gamma_%g_omega_%g_electron_densities.dat"%(gamma_good, omega_good))
	sh_good_time = sh_good[:, 0]
	sh_good_cdw = sh_good[:, 3] - sh_good[:, 4]
	
	sh_bad = np.loadtxt("data_export/sh_gamma_%g_omega_%g_electron_densities.dat"%(gamma_bad, omega_bad))
	sh_bad_time = sh_bad[:, 0]
	sh_bad_cdw = sh_bad[:, 3] - sh_bad[:, 4]
	
	sh_good_spacing = 10
	sh_bad_spacing = 10
	
	mce_style = "x"
	sh_style = "+"

exact_style = "-"
mte_style = "o"

mte_good_spacing = 2*mce_good_spacing
mte_bad_spacing = 100 * mce_good_spacing

ref_good = data_manager.RefData(initial_state, psi_initial, num_sites, omega_good, gamma_good)
ref_good.set_cut_time(t_max)

ref_bad = data_manager.RefData(initial_state, psi_initial, num_sites, omega_bad, gamma_bad)
ref_bad.set_cut_time(t_max)

sim_good = data_manager.Simulation(omega_good, gamma_good, t_max,
                                   num_samples=num_samples, num_sites=num_sites,
                                   symmetric_y_init=False,
                                   num_electrons=num_electrons,
                                   printProgress=True,
                                   initial_state=initial_state,
                                   psi_initial=psi_initial,
                                   calcYSquared=True,
                                   fullStorage=False,
                                   hopping=hopping)


sim_bad = data_manager.Simulation(omega_bad, gamma_bad, t_max, time_step=5e-4,
                                   num_samples=num_samples, num_sites=num_sites,
                                   symmetric_y_init=False,
                                   num_electrons=num_electrons,
                                   printProgress=True,
                                   initial_state=initial_state,
                                   psi_initial=psi_initial,
                                   calcYSquared=True,
                                   fullStorage=False,
                                   hopping=hopping)


print("mte_good_dt", sim_good.time()[1], sim_good.time()[1]*mte_good_spacing)
print("mte_bad_dt", sim_bad.time()[1], sim_bad.time()[1]*mte_bad_spacing)

if num_sites == 2:
	figsize = (9,3)
else:
	figsize = (8,3)

fig, axs = plt.subplots(1, 2, figsize=figsize, sharex=False)
ax_good = axs[0]
ax_bad = axs[1]

#text_x = 1.5
#text_y = .98

ax_good.plot(ref_good.time(), ref_good.O_CDW(), exact_style)
ax_good.plot(sim_good.time()[::mte_good_spacing],
             sim_good.O_CDW()[::mte_good_spacing],
             mte_style, mfc="none")
#ax_good.plot([], [])
if num_sites == 2:
	print("sh_good_dt", sh_good_time[1], sh_good_time[1]*sh_good_spacing)
	ax_good.plot(sh_good_time[::sh_good_spacing],
	             sh_good_cdw[::sh_good_spacing],
	             sh_style)
	mce_good_time = np.loadtxt(mce_path_bad + "time")
	print("mce_good_dt", mce_good_time[1], mce_good_time[1]*mce_good_spacing)
	mce_good_nel = np.loadtxt(mce_path_good + "electron_density")
	ax_good.plot(mce_good_time[::mce_good_spacing],
	             mce_good_nel[::mce_good_spacing, 0]-mce_good_nel[::mce_good_spacing, 1],
	             mce_style)
	

ax_good.set_ylim(ymax=1.15)

ax_good.set_xlim(0, t_max)
ax_good.set_xlabel(r"$t\cdot t_0$")

if num_sites == 2:
	ax_good.text(1.5, 1.02, "$L=%i, \\omega_0=%gt_0$, $\\gamma=%gt_0$"%(num_sites, omega_good, gamma_good))
	ax_good.set_ylabel(r"$\langle n_1\rangle - \langle n_2\rangle$")
	
	ax_bad.text(1.5, .97, "$L=%i, \\omega_0=%gt_0$, $\\gamma=%gt_0$"%(num_sites, omega_bad, gamma_bad))
	ax_bad.set_ylabel(r"$\langle n_1\rangle - \langle n_2\rangle$")
else:
	ax_good.text(.1, 1.02, "$L=%i, N=%i, \\omega_0=%gt_0$, $\\gamma=%gt_0$"%(num_sites, num_sites//2, omega_good, gamma_good))
	ax_good.set_ylabel(r"$\mathcal{O}_{\mathsf{CDW}} = \frac{1}{N}\sum_{i=1}^L (-1)^i\langle n_i\rangle$")
	ax_bad.text(1., .97, "$L=%i, N=%i, \\omega_0=%gt_0$, $\\gamma=%gt_0$"%(num_sites, num_sites//2, omega_bad, gamma_bad))
	ax_bad.set_ylabel(r"$\mathcal{O}_{\mathsf{CDW}} = \frac{1}{N}\sum_{i=1}^L (-1)^i\langle n_i\rangle$")


ax_bad.plot(ref_bad.time(), ref_bad.O_CDW(), exact_style, label=("exact" if num_sites == 2 else "DMRG"))
ax_bad.plot(sim_bad.time()[::mte_bad_spacing],
            sim_bad.O_CDW()[::mte_bad_spacing],
            mte_style, label="Ehrenfest" if num_samples <= 1 else "MTE", mfc="none")

if num_sites == 2:
	print("sh_bad_dt", sh_bad_time[1], sh_bad_time[1]*sh_bad_spacing)
	ax_bad.plot(sh_bad_time[::sh_bad_spacing],
	            sh_bad_cdw[::sh_bad_spacing],
	            sh_style, label="surface hopping")
	mce_bad_time = np.loadtxt(mce_path_bad + "time")
	print("mce_bad_dt", mce_bad_time[1], mce_bad_time[1]*mce_bad_spacing)
	mce_bad_nel = np.loadtxt(mce_path_bad + "electron_density")
	ax_bad.plot(mce_bad_time[::mce_bad_spacing],
	            mce_bad_nel[::mce_bad_spacing, 0]-mce_bad_nel[::mce_bad_spacing, 1],
	            mce_style, label="MCE")
	

ax_bad.set_xlim(0, t_max)
ax_bad.set_xlabel(r"$t\cdot t_0$")

if num_sites == 2:
	fig.legend(loc="right", framealpha=1.)
	fig.tight_layout(rect=(0., 0., 0.82, 1.))
else:
	fig.legend(loc="upper center", framealpha=1.)
	fig.tight_layout(rect=(0., 0., 1., 1.))

fig.savefig("../crc_plots/ehrenfest_%s_vs_%s_L=%i_N=%i_(%g,%g)_(%g,%g)_%s_%s%s.pdf"\
            %("single" if num_samples <= 1 else "multi", 
              "ED" if num_sites == 2 else "DMRG",
              num_sites, num_electrons, omega_good, gamma_good, omega_bad, gamma_bad,
              initial_state, psi_initial,
              ("_cube_interpol" if cubic_interpolate else "")))

#plt.show()
