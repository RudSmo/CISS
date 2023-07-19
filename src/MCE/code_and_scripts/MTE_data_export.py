import numpy as np
import matplotlib.pyplot as plt
import os

import data_manager

initial_state = "DCDW"
psi_initial="local"
hopping = 10
omega = 1
gamma = 4
num_sites_x = 2
num_sites_y = 1
num_electrons = 1

num_samples = 0
symmetrize_start = False

t_max = 100
dt = 0.0025
save_dt = dt #0.025

sim = data_manager.Simulation(omega, gamma, t_max, time_step=dt,
                              num_samples=num_samples,
                              num_sites_x=num_sites_x, num_sites_y=num_sites_y,
                              symmetric_y_init=symmetrize_start,
                              num_electrons=num_electrons,
                              printProgress=True,
                              initial_state=initial_state,
                              psi_initial=psi_initial,
                              calcYSquared=True,
                              fullStorage=False,
                              hopping=hopping)

save_stepping = int(np.round(save_dt/sim.time()[1]))
print("save_stepping", save_stepping)

path = "data_export/" + "%s_psi_initial=%s_hopping=%g_omega=%g_gamma=%g_mte_num_traj=%i/"\
       %(initial_state, psi_initial, hopping, omega, gamma, num_samples)

os.system("mkdir -p "+path)

np.savetxt(path+"time", sim.time()[::save_stepping])
np.savetxt(path+"x", sim.x()[::save_stepping])
np.savetxt(path+"p", sim.p()[::save_stepping])
np.savetxt(path+"Nph", sim.nph()[::save_stepping])
np.savetxt(path+"q", sim.q()[::save_stepping], header="q = x0-x1")
np.savetxt(path+"pq", sim.pq()[::save_stepping], header="pq = p0-p1")
np.savetxt(path+"electron_density", sim.nel()[::save_stepping], header="electronDensity on site ...")
np.savetxt(path+"epot_ph", sim.x2()[::save_stepping]*omega/2, header="potential energy per site")
np.savetxt(path+"ekin_ph", sim.p2()[::save_stepping]*omega/2, header="potential energy per site")
np.savetxt(path+"ekin_el", sim.energy_kinetic()[::save_stepping], header="kinetic energy electron")
np.savetxt(path+"energy_interaction", sim.energy_interaction(site="all")[::save_stepping], header="el-ph interaction per site")

