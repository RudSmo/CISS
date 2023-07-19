import numpy as np
import matplotlib.pyplot as plt

from src.ehrenfest_init import initial_state_MCE
from simulate import simulate_system_MCE

import data_manager

initial_state = "BCDW"
psi_initial = "local"
psi_state = None
num_electrons=1
hopping = 0.5
omega = 1
gamma = np.sqrt(2)
num_sites_x = 2
num_sites_y = 1

symmetrize_start = True
num_samples = 1000

num_configs = 80	# must be multiple of 4 if symmetrize_start (x and p are mirrored)
max_num_configs = num_configs
regularize_strength = 1e-9		# increase for larger number of configs

regularize_diag = False


dt_MCE = 0.05
t_max = 20


rubio_exact_time_nel = np.loadtxt("ref_data/exact_prb_rubio.dat")
rubio_MTE_time_nel = np.loadtxt("ref_data/MTE_prb_rubio.dat")

simMTE = data_manager.Simulation(omega, gamma, t_max,
                                 num_samples=num_samples,
                                 num_sites_x=num_sites_x, num_sites_y=num_sites_y,
                                 symmetric_y_init=symmetrize_start,
                                 num_electrons=num_electrons,
                                 printProgress=True,
                                 initial_state=initial_state,
                                 psi_initial=psi_initial,
                                 calcYSquared=False,
                                 fullStorage=False,
                                 hopping=hopping)

simMCE = data_manager.Simulation_MCEv1(hopping, omega, gamma,
                                       num_sites_x, num_sites_y,
                                       t_max, dt_MCE,
                                       num_configs,
                                       initial_state, psi_initial, psi_state,
                                       max_num_configs, printProgress=True)

title = initial_state + " " + psi_initial + ": $t_0=%g, \\omega=%g$, $\\gamma=%g,$ conf %i$\\rightarrow$%i"%(hopping, omega, gamma, num_configs, max_num_configs)

plt.title("rubio compare electron density, " +title)

plt.plot(simMTE.time(), simMTE.O_CDW(), '-', label="own MTE")
plt.plot(simMCE.time(), simMCE.O_CDW(), '-', label="MCE")
plt.plot(rubio_exact_time_nel[:, 0], rubio_exact_time_nel[:, 1], '.', label="rubio exact")
plt.plot(rubio_MTE_time_nel[:, 0], rubio_MTE_time_nel[:, 1], '.', label="rubio MTE")

plt.xlabel("t")
plt.ylabel(r"$\langle \hat{n}_{el,1}\rangle - \langle \hat{n}_{el,2}\rangle$")
plt.legend()
plt.show()



