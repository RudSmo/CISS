import os
import numpy as np
import matplotlib.pyplot as plt
# Correct latex
# latex style
use_tex = True
if use_tex:
	from matplotlib import rc
	rc('text', usetex=use_tex)
	rc('font', family='serif')

import data_manager


initial_state = "dressed"
psi_initial="CDW"
hopping = 10
omega = 1
gamma = 4
num_sites_x = 13
num_sites_y = 1
figsize=(5, 7)
if omega != 1:
	raise ValueError("omega must be 1")
 
# number of ehrenfest trajectories (randomly chosen)
# use num_samples = 0 for single trajectory starting in rest on all sites
num_samples = 4000
symmetrize_start = False


# dt is the save distance, simulation distance may be smaller
dt = 0.002
#t_max = 6.5/hopping #20*np.pi
# error in energy conservation \approx 0.2*omega*gamma^2*dt^2 for BCDW
# DCDW comparable or lower


try:
	ref = data_manager.RefData(initial_state, psi_initial, None, num_sites_x, hopping, omega, gamma)
except Exception as e:
	print("no reference Data found:")
	print(repr(e))
	quit()

t_max = ref.time()[-1]

file_name = "figures/extended_MTE_DMRG/"
os.system("mkdir -p " + file_name)
file_name += f"{initial_state}_hopping=%g_gamma=%g.pdf"%(hopping, gamma)

title = "$L=%i$, $N=%i$, $\\bar\\gamma=%g,$ $\\bar t_0=%g$, "\
        %(num_sites_x, num_sites_x//2, gamma, hopping)
if initial_state == "bare":
	title += "BCDW"
else:
	title += "DCDW"

time_label = "$t\\ [1/\\omega_0]$"

sim = data_manager.Simulation_MTE(hopping, omega, gamma,
                                  num_sites_x, num_sites_y,
                                  t_max, dt,
                                  num_samples,
                                  initial_state, psi_initial, psi_state=None,
                                  phonon_dressed_displaced_sites=None,
                                  calcYSquared=True,
                                  fullStorage=False,
                                  symmetric_y_init=symmetrize_start,
                                  force_new_simulation=False,
                                  load_sampled_state=False,
                                  store_results=True)

MTE_label = "MTE"
DMRG_label = "DMRG"

fig, axs = plt.subplots(3, 1, figsize=figsize, sharex=True)
fig.suptitle(title)
ax_nph  = axs[0]
ax_OCDW = axs[1]
ax_ekin = axs[2]
tick_params = {"axis":"both", "direction":"in", "left":True, "bottom":True, "right":True, "top":True}
ylabel_coord = (-0.08, 0.5)
ax_nph .tick_params(**tick_params)
ax_OCDW.tick_params(**tick_params)
ax_ekin.tick_params(**tick_params)

ax_nph .get_yaxis().set_label_coords(*ylabel_coord)
ax_OCDW.get_yaxis().set_label_coords(*ylabel_coord)
ax_ekin.get_yaxis().set_label_coords(*ylabel_coord)

# phonons
ax_nph.plot(sim.time(), sim.nph().sum(axis=1)/sim.num_electrons, label="MTE")
if ref is not None:
	ax_nph.plot(ref.time(), ref.nph().sum(axis=1)/ref.num_electrons, '--', label=DMRG_label)
ax_nph.set_ylabel(r"$N_{ph}/N$")


# O_CDW
O_CDW, O_CDW_eq = sim.O_CDW(True)
ax_OCDW.plot(sim.time()[[0,-1]], [O_CDW_eq, O_CDW_eq], 'k')
#ax_OCDW.plot(sim.time()[[0,-1]], [0, 0], 'k')
ax_OCDW.plot(sim.time(), O_CDW, label=MTE_label)
if ref is not None:
	ax_OCDW.plot(ref.time(), ref.O_CDW(), '--', label=DMRG_label)

ax_OCDW.set_ylabel(r"$\mathcal{O}_{CDW}$")


# kinetic energy
ax_ekin.plot(sim.time(), sim.energy_kinetic(), label=MTE_label)
if ref is not None:
	ax_ekin.plot(ref.time(), ref.energy_kinetic(), '--', label=DMRG_label)
ax_ekin.set_ylabel("$E_{kin}$")  # r"$\langle H_{el, kin}\rangle$")


annotate_pos = (0.95, 0.5)

axs[0].annotate('a', annotate_pos, xycoords='axes fraction')
axs[1].annotate('b', annotate_pos, xycoords='axes fraction')
axs[2].annotate('c', annotate_pos, xycoords='axes fraction')

ax_OCDW.legend()
ax_ekin.set_xlim((0, sim.time()[-1]))
ax_ekin.set_xlabel(time_label)
fig.subplots_adjust(left=.12, bottom=.06, right=.99, top=.94, wspace=None, hspace=0.)
fig.savefig(file_name)

#plt.show()

