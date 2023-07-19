import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import os

from src.ehrenfest_init import initial_state_MCE
from simulate import simulate_system_MCE

import data_manager

fast = False
plot_nel = False

if fast:
	hopping = 10
	omega = 1
	gamma = 4
	t_max = None
	num_configs = 1200
else:
	hopping = .5
	omega = 1
	gamma = 2
	t_max = 20
	num_configs = 1500


num_traj = 11000

num_sites_x = 11
num_sites_y = 1
num_electrons = 1

max_num_configs = num_configs

dt = 0.01  # 0.0001
initial_state = "BCDW"
psi_initial = "local"

save_plots = True
figsize = (8, 4.5)
cmap="Greys"
x = np.arange(1, num_sites_x+1)
y = np.arange(1, num_sites_y+1)

if plot_nel:
	vmin=-2.5
	vmax=0
else:
	vmin=0
	vmax=None

"""
hoz, el_coeff = initial_state_MCE(num_sites_x, num_sites_y, num_configs, num_configs, num_electrons, 
                                  False, initial_state, psi_initial,
                                  hopping, omega, gamma,
                                  False, 0,
                                  plotOverlapMatrix=False, use_QSS=True, random_hoz=True,
                                  cutoff_eigval_overlap=0,
                                  useDisorderRegularization=False)
"""

ref_path = "ref_data/DMRG2/holtdmrglbo_single_ebcdwL%iM40site_in6Md3000lboMd100t0%.1fgam%.1fomg%.1fdt0.01tot%.1flbocut1e-09cut1e-09DIR/"%(num_sites_x, hopping, gamma, omega, 20 if gamma == 2 else 10)
ref_time = np.loadtxt(ref_path + "time.dat")

ref_nel = np.empty((ref_time.size, num_sites_x))
ref_nph = np.empty((ref_time.size, num_sites_x))
for s in range(num_sites_x):
	ref_nel[:, s] = np.loadtxt(ref_path+"nlocSite%i.dat"%(s+1))
	ref_nph[:, s] = np.loadtxt(ref_path+"NphlocSite%i.dat"%(s+1))


mte_path = "data_export/%s_psi_initial=%s_hopping=%g_omega=%g_gamma=%g_mte_num_traj=%i/"\
           %(initial_state, psi_initial, hopping, omega, gamma, num_traj)

mte_time = np.loadtxt(mte_path+"time")
mte_nel = np.loadtxt(mte_path +"electron_density")
mte_nph = np.loadtxt(mte_path +"Nph")


# simulation
#obs_y, obs_nph, obs_nel, obs_el_ekin, obs_el_ph_inter = simulate_system_MCE(hoz, el_coeff, hopping, omega, gamma, num_times, dt, num_configs, printProgress=True)

path = path = "data_export/" + "%s_psi_initial=%s_hopping=%g_omega=%g_gamma=%g_mce_numconf=%i-%i_num_sites_x=%i_num_sites_y=%i/"\
       %(initial_state, psi_initial, hopping, omega, gamma, num_configs, max_num_configs, num_sites_x, num_sites_y)

obs_nel = np.load(path + "electron_density.npy")
obs_nph = np.load(path + "Nph.npy")
time = np.load(path + "time.npy")
e_tot = np.load(path + "energy_total.npy")

if t_max is None:
	t_max = time[-1]

obs_nel = obs_nel.reshape((time.size, num_sites_y, num_sites_x))
obs_nph = obs_nph.reshape((time.size, num_sites_y, num_sites_x))


fig=plt.figure(figsize=figsize)
ax0 = fig.add_axes([0.07, 0.1, 0.25, 0.84])
ax1 = fig.add_axes([0.39, 0.1, 0.25, 0.84])
ax2 = fig.add_axes([0.64, 0.1, 0.25, 0.84], sharex=ax1)
cbar_ax = fig.add_axes([0.9, 0.1, 0.02, 0.84])

"""
ax0 = plt.subplot(131)
ax1 = plt.subplot(132)
ax2 = plt.subplot(133, sharex = ax1)
"""
if plot_nel:
	ax0.plot(time, obs_nel[:, 0, num_sites_x//2], label="MCE %i"%num_configs)
	ax0.plot(ref_time, ref_nel[:, num_sites_x//2], '--', label="DMRG")
	ax0.plot(mte_time, mte_nel[:, num_sites_x//2], '-.', label="MTE")
	cm = ax1.pcolormesh(x, time, np.log10(obs_nel[:, 0, :]), cmap=cmap, vmin=vmin, vmax=vmax, shading="nearest")
	ax2.pcolormesh(x, ref_time, np.log10(ref_nel), cmap=cmap, vmin=vmin, vmax=vmax, shading="nearest")
else:
	if vmax is None:
		vmax = obs_nph.max()
	ax0.plot(time, obs_nph[:, 0, num_sites_x//2], label="MCE %i"%num_configs)
	ax0.plot(ref_time, ref_nph[:, num_sites_x//2], '--', label="DMRG")
	ax0.plot(mte_time, mte_nph[:, num_sites_x//2], '-.', label="MTE")
	cm = ax1.pcolormesh(x, time, obs_nph[:, 0, :], cmap=cmap, vmin=vmin, vmax=vmax, shading="nearest")
	ax2.pcolormesh(x, ref_time, ref_nph, cmap=cmap, vmin=vmin, vmax=vmax, shading="nearest")

left = num_sites_x//2-.5 + 1
right = num_sites_x//2+.5 + 1
ax1.plot([left, right, right, left, left], [0, 0, t_max, t_max, 0], color="tab:blue")
ax2.plot([left, right, right, left, left], [0, 0, t_max, t_max, 0], '--', color="tab:orange")


ax0.set_xlim((-5*dt, t_max+5*dt))
ax1.set_ylim((t_max+5*dt, -5*dt))
ax2.set_ylim((t_max+5*dt, -5*dt))

#ax0.set_title(r"$L=11, t_0 = \omega_0/2, \gamma=2\omega_0$")
ax1.set_title("MCE %i"%num_configs)
ax2.set_title("DMRG (David)")

ax0.set_xlabel(r"$\omega_0t$")
if plot_nel:
	ax0.set_ylabel(r"$\langle n_{el, %i}\rangle$"%(num_sites_x//2 +1))
	ax0.set_ylim((-0.01, 1.01))
else:
	ax0.set_ylabel(r"$\langle N_{ph, %i}\rangle$"%(num_sites_x//2 +1))
	ax0.set_ylim(ymin=-0.01, ymax = vmax+0.01)
ax0.legend()

ax1.set_ylabel(r"$\omega_0t$")
ax2.set_yticklabels([])

ax1.set_xticks(np.arange(1, num_sites_x+1), minor=True)
ax1.set_xticks(np.arange(2, num_sites_x+1, 4), minor=False)

yticks = np.arange(0, t_max+.5, 4)
ax1.set_yticks(yticks)
ax2.set_yticks(yticks)
ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
ax1.set_xlabel("site")
ax2.set_xlabel("site")

cbar = fig.colorbar(cm, cax=cbar_ax)
if plot_nel:
	cbar.set_label(r"$\log_{10}(\langle n_{el}\rangle)$")
else:
	cbar.set_label(r"$\langle N_{ph}\rangle$")



plt.figure()
plt.title(r"MCE: $E_{tot}(t)-(E_{tot}(0)=%g)$ "%e_tot[0])
e_tot -= e_tot[0]
plt.plot(time, e_tot)

if save_plots:
	path = "../21_05_03_main_research_presentation/images/Multi_Config_Num_%g_%s_psi_initial=%s_hopping=%g_omega=%g_gamma=%g_num_sites_x=%i_num_sites_y=%i"\
	       %(num_configs, initial_state, psi_initial, hopping, omega, gamma, num_sites_x, num_sites_y)
	
	plt.savefig(path + "_etot.png", dpi=300)
	fig.savefig(path + "_%s.png"%("nel" if plot_nel else "nph"), dpi=500)
else:
	plt.show()
