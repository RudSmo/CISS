import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from src.ehrenfest_init import initial_state_MCE
from simulate import simulate_system_MCE

import data_manager

hopping = 1
omega = 2
gamma = 4
num_sites_x = 2
num_sites_y = 1

num_configs = 40
max_num_configs = num_configs

time_step = 0.05
time_max = 1
time_index = -1

initial_state = "BCDW"
psi_initial = "local"
psi_state = None

markersize = 16
q = np.linspace(-8, 8, 101)
pq = np.linspace(-4, 4, 101)

sim = data_manager.Simulation_MCE(2, hopping, omega, gamma,
                                  num_sites_x, num_sites_y,
                                  time_max, time_step,
                                  num_configs,
                                  initial_state, psi_initial, psi_state,
                                  max_num_configs, fullStorage=True)

time_index = (sim.time().size+time_index)%(sim.time().size)
title = initial_state + " " + psi_initial + ":$t_0=%g$, $\\omega=%g$, $\\gamma=%g$, t=%g "\
        %(hopping, omega, gamma, sim.time()[time_index])

q, pq, wigner = sim.get_Wignerfunction_qp(time_index, q, pq)

dq = q[1]-q[0]
dpq = pq[1]-pq[0]

print("integral wigner=", dq*dpq*np.sum(wigner))

vmax = 1.1*wigner.max()
vmin = -vmax/4.

old_cm = plt.get_cmap('RdGy', 8*256)
newcolors = old_cm(np.linspace(3/8, 1, 5*256))
newcolors[:256, :] = old_cm(np.linspace(0, 1/2, 256))
newcmp = ListedColormap(newcolors)


Q, PQ = np.meshgrid(q, pq)

plt.title("wigner function, " + title)
plt.pcolormesh(Q, PQ, wigner, vmin=vmin, vmax=vmax, cmap=newcmp, shading="nearest") #gouraud
#plt.pcolormesh(Q, PQ, wigner, vmin=-vmax, vmax=vmax, cmap="RdGy", shading="nearest") #gouraud
plt.colorbar()

plt.figure()
plt.title("config locations, " + title)
el_coeff = sim.fullStorage_el_coeff[time_index]
hoz = sim.fullStorage_hoz[time_index]

hoz = hoz[0, 0, :] - hoz[0, 1, :]
el_coeff = np.linalg.norm(el_coeff, axis=(0, 1))
order = np.argsort(el_coeff[:num_configs])

plt.scatter(np.sqrt(2)*hoz[order].real, np.sqrt(2)*hoz[order].imag,
            c=el_coeff[order], vmin=0, vmax=el_coeff[order[-1]], s=markersize)
plt.colorbar()
plt.xlim((q[0], q[-1]))
plt.ylim((pq[0], pq[-1]))

plt.show()

