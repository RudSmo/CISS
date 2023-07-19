import numpy as np
import matplotlib.pyplot as plt
import os

from src.ehrenfest_init import initial_state_MCE
from simulate import simulate_system_MCE

import data_manager

hopping = 1
omega = 2
gamma = 4
num_sites_x = 2
num_sites_y = 1
num_electrons = 1

num_configs_list = [10, 30, 60]
linestyles = ['-', '--', '-.']


dt = 0.1  # 0.0001
t_max = 10

initial_state = "BCDW"
psi_initial = "local"

save_plots = True
figsize = (6, 3)


ref = data_manager.RefData(initial_state, psi_initial, num_sites_x, omega, gamma)
ref.set_cut_time(t_max)

num_times = int(np.round(t_max/dt)) + 1  # +1 for t=0
time = np.linspace(0, t_max, num_times)/hopping*omega

fig, axs = plt.subplots(1, 2, figsize=figsize, sharex=True)
ax_nel = axs[0]
#ax_nph = axs[1]
ax_etot = axs[1]

ax_nel.plot(ref.time()/hopping*omega, ref.nel(0), 'k', label="exact")
#ax_nph.plot(ref.time(), ref.nph(0), 'k', label="exact")
ax_etot.plot(ref.time()/hopping*omega, np.zeros(ref.time().size), 'k', label="exact")

for i in range(len(num_configs_list)):
	hoz, el_coeff = initial_state_MCE(num_sites_x, num_sites_y, num_configs_list[i], num_configs_list[i], num_electrons, 
	                                  False, initial_state, psi_initial,
	                                  hopping, omega, gamma,
	                                  False, 0,
	                                  plotOverlapMatrix=False, use_QSS=True, random_hoz=True,
	                                  cutoff_eigval_overlap=0,
	                                  useDisorderRegularization=False)
	
	# simulation
	obs_y, obs_nph, obs_nel, obs_el_ekin, obs_el_ph_inter = simulate_system_MCE(hoz, el_coeff, hopping, omega, gamma, num_times, dt, num_configs_list[i], printProgress=True)
	
	ax_nel.plot(time, obs_nel[:, 0, 0], linestyles[i], label="MCE %i"%num_configs_list[i])
	#ax_nph.plot(time, obs_nph[:, 0, 0], linestyles[i], label="MCE %i"%num_configs_list[i])
	
	
	e_tot =   omega * obs_nph.sum(axis=(1, 2))\
	        + gamma * obs_el_ph_inter.sum(axis=(1, 2))\
	        + hopping * obs_el_ekin.sum(axis=1)
	#plt.title(r"$\Delta E_{tot}$ " + title)
	e_tot -= e_tot[0]
	ax_etot.plot(time, e_tot/omega, linestyles[i], label="MCE %i"%num_configs_list[i])


ax_nel.set_xlim((-.1, time[-1]+.1))
ax_nel.set_ylim((-.01, 1.01))
ax_nel.set_xlabel(r"$\omega_0 t$")
ax_nel.set_ylabel(r"$\langle n_{el,1}\rangle$")

#ax_nph.set_ylim(ymin=-.05*ref.nph(0).max())

ax_etot.set_xlabel(r"$\omega_0 t$")
ax_etot.set_ylabel("$E_{tot}/\\hbar\\omega_0$")
ax_etot.legend()

fig.tight_layout()


if save_plots:
	path = "../21_05_03_main_research_presentation/images/Multi_Config_Num_%s_psi_initial=%s_hopping=%g_omega=%g_gamma=%g_num_sites_x=%i_num_sites_y=%i.png"\
	       %(initial_state, psi_initial, hopping, omega, gamma, num_sites_x, num_sites_y)
	fig.savefig(path, dpi=300)


plt.show()
