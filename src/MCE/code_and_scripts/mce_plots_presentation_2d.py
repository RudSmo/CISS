import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import os

from src.ehrenfest_init import initial_state_MCE
from simulate import simulate_system_MCE

import data_manager

hopping = 1
omega = 1
gamma = 1
num_sites_x = 11
num_sites_y = 11
num_electrons = 1

num_configs = 1000

dt = 0.1  # 0.0001
t_max = 3.5

initial_state = "BCDW"
psi_initial = "local"

save_plots = True
figsize = (8, 4)
cmap="Greys"
vmin=-2.5
vmax=0
x = np.arange(1, num_sites_x+1)
y = np.arange(1, num_sites_y+1)

num_times = int(np.round(t_max/dt)) + 1  # +1 for t=0
time = np.linspace(0, t_max, num_times)/hopping*omega

time_stepping = num_times//5

hoz, el_coeff = initial_state_MCE(num_sites_x, num_sites_y, num_configs, num_configs, num_electrons, 
                                  False, initial_state, psi_initial,
                                  hopping, omega, gamma,
                                  False, 0,
                                  plotOverlapMatrix=False, use_QSS=True, random_hoz=True,
                                  cutoff_eigval_overlap=0,
                                  useDisorderRegularization=False)

# simulation
obs_y, obs_nph, obs_nel, obs_el_ekin, obs_el_ph_inter = simulate_system_MCE(hoz, el_coeff, hopping, omega, gamma, num_times, dt, num_configs, printProgress=True)



fig, axs = plt.subplots(2, 3, figsize=figsize, sharex=True, sharey=True)
cm = axs[0, 0].pcolormesh(x, y, np.log10(obs_nel[0*time_stepping, :, :]), cmap=cmap, vmin=vmin, vmax=vmax, shading="nearest")
axs[0, 1].pcolormesh(x, y, np.log10(obs_nel[1*time_stepping, :, :]), cmap=cmap, vmin=vmin, vmax=vmax, shading="nearest")
axs[0, 2].pcolormesh(x, y, np.log10(obs_nel[2*time_stepping, :, :]), cmap=cmap, vmin=vmin, vmax=vmax, shading="nearest")

axs[1, 0].pcolormesh(x, y, np.log10(obs_nel[3*time_stepping, :, :]), cmap=cmap, vmin=vmin, vmax=vmax, shading="nearest")
axs[1, 1].pcolormesh(x, y, np.log10(obs_nel[4*time_stepping, :, :]), cmap=cmap, vmin=vmin, vmax=vmax, shading="nearest")
axs[1, 2].pcolormesh(x, y, np.log10(obs_nel[5*time_stepping, :, :]), cmap=cmap, vmin=vmin, vmax=vmax, shading="nearest")

axs[0, 0].set_xlim((.5, num_sites_x+.5))
axs[0, 0].set_ylim((.5, num_sites_y+.5))
axs[0, 0].set_xticks(x)
axs[0, 0].set_yticks(y)
"""
axs[0, 1].set_xlim((.5, num_sites_x+.5))
axs[0, 1].set_ylim((.5, num_sites_y+.5))

axs[0, 2].set_xlim((.5, num_sites_x+.5))
axs[0, 2].set_ylim((.5, num_sites_y+.5))

axs[1, 0].set_xlim((.5, num_sites_x+.5))
axs[1, 0].set_ylim((.5, num_sites_y+.5))

axs[1, 1].set_xlim((.5, num_sites_x+.5))
axs[1, 1].set_ylim((.5, num_sites_y+.5))
"""

axs[0, 0].set_title(r"$\omega_0 t = %g$"%(omega*time[0*time_stepping]))
axs[0, 1].set_title(r"$\omega_0 t = %g$"%(omega*time[1*time_stepping]))
axs[0, 2].set_title(r"$\omega_0 t = %g$"%(omega*time[2*time_stepping]))

axs[1, 0].set_title(r"$\omega_0 t = %g$"%(omega*time[3*time_stepping]))
axs[1, 1].set_title(r"$\omega_0 t = %g$"%(omega*time[4*time_stepping]))
axs[1, 2].set_title(r"$\omega_0 t = %g$"%(omega*time[5*time_stepping]))

fig.subplots_adjust(left=0.05, bottom=.06, right=0.87, top=.94, wspace=0.1)
cbar_ax = fig.add_axes([0.9, 0.06, 0.02, 0.88])
cbar = fig.colorbar(cm, cax=cbar_ax)
cbar.set_label(r"$\log_{10}(\langle n\rangle)$")

plt.figure()
e_tot =   omega * obs_nph.sum(axis=(1, 2))\
        + gamma * obs_el_ph_inter.sum(axis=(1, 2))\
        + hopping * obs_el_ekin.sum(axis=1)
plt.title(r"MCE: $E_{tot}(t)-(E_{tot}(0)=%g)$ "%e_tot[0])
e_tot -= e_tot[0]
plt.plot(time, e_tot)

if save_plots:
	path = "../21_05_03_main_research_presentation/images/Multi_Config_Num_%g_%s_psi_initial=%s_hopping=%g_omega=%g_gamma=%g_num_sites_x=%i_num_sites_y=%i"\
	       %(num_configs, initial_state, psi_initial, hopping, omega, gamma, num_sites_x, num_sites_y)
	
	plt.savefig(path + "_etot.png", dpi=300)
	fig.savefig(path + ".png", dpi=300)


plt.show()
