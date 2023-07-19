import numpy as np
import matplotlib.pyplot as plt

from src.ehrenfest_init import initial_state_MCE
from simulate import simulate_system_MCE

import data_manager

initial_state = "BCDW"
psi_initial = "local"
num_electrons=1
hopping = 1.
omega, gamma = .1, .4
#omega, gamma = .4, 1.2
#omega, gamma = 2, 4
num_sites = 2

random_hoz = True
symmetrize_start = True

num_samples = 1000

num_configs = 40	# must be multiple of 4 if symmetrize_start (x and p are mirrored)
regularize_strength = 1e-9		# increase for larger number of configs

regularize_diag = False


dt_MCE = 0.0005
t_max = 10 #np.pi*20


try:
	ref = data_manager.RefData("BCDW", "local", num_sites, omega, gamma)
	ref.set_cut_time(t_max)
	t_max = ref.time()[-1]
except Exception as e:
	print("no reference Data found:")
	print(repr(e))
	ref = None

num_times = int(np.round(t_max/dt_MCE)) + 1  # +1 for t=0

# initialization
electron_sites = []
for i in range(num_electrons):
	site = 2*i+(num_sites+1)//2-num_electrons
	electron_sites.append(site)

psi_init = np.zeros((num_electrons, num_sites), dtype=np.complex)
for i in range(len(electron_sites)):
	psi_init[i, electron_sites[i]] = 1

hoz, el_coeff = initial_state_MCE(num_sites, num_configs, 1, symmetrize_start,
                                  np.zeros((num_sites, 2)), psi_init,
                                  regularize_diag, regularize_strength)

obs_y, obs_nph, obs_nel, obs_el_ekin, obs_el_ph_inter = simulate_system_MCE(hoz, el_coeff, hopping, omega, gamma, num_times, dt_MCE, printProgress=False)

simMTE = data_manager.Simulation(omega, gamma, t_max,
                                 num_samples=num_samples, num_sites=num_sites,
                                 symmetric_y_init=symmetrize_start,
                                 num_electrons=num_electrons,
                                 printProgress=True,
                                 initial_state=initial_state,
                                 psi_initial=psi_initial,
                                 calcYSquared=True,
                                 fullStorage=False,
                                 hopping=hopping)

time = np.linspace(0, t_max, num_times)


fig, axs = plt.subplots(1, 2, figsize=(5, 4))
#fig.suptitle("$\\omega=%g$, $\\gamma=%g$"%(omega, gamma))

print("omega=%g, gamma=%g"%(omega, gamma))
ref_color= "gray"

axs[0].set_title(r"electron density $\langle N_{el, 1}\rangle$")
axs[0].plot(ref.time(), ref.nel(0))
#axs[0].plot(ref.time(), ref.nel(1))

axs[0].plot(time, obs_nel[:, 0], "--")
#axs[0].plot(time, obs_nel[:, 1], "--")
axs[0].plot(simMTE.time(), simMTE.nel(0), "-.")
#axs[0].plot(simMTE.time(), simMTE.nel(1), "-.")

axs[0].set_xlabel("time $t*t_0$")
#axs[0].set_ylabel(r"$\langle N_{el, 1}\rangle$")

axs[0].set_xlim((0, t_max))
axs[0].set_ylim((0, 1))


axs[1].set_title(r"phonon number $\langle N_{ph}\rangle$")

axs[1].plot(ref.time(), ref.nph(0), label="exact")
axs[1].plot(time, obs_nph[:, 0], "--", label="MCE %i configurations"%num_configs)
axs[1].plot(simMTE.time(), simMTE.nph(0), "-.", label="MTE %i trajectories"%num_samples)

axs[1].set_prop_cycle(None)
axs[1].plot(ref.time(), ref.nph(1))
axs[1].plot(time, obs_nph[:, 1], "--")
axs[1].plot(simMTE.time(), simMTE.nph(1), "-.")

axs[1].set_xlabel("time $t*t_0$")
#axs[1].set_ylabel(r"$\langle N_{ph, 1}\rangle$")

axs[1].set_xlim((0, t_max))
axs[1].set_ylim(ymin=0)



"""
axs[2].set_title(r"phonon number site 2 $\langle N_{ph, 2}\rangle$")
axs[2].set_xlabel("time $t*t_0$")
axs[2].set_ylabel(r"$\langle N_{ph, 2}\rangle$")

axs[2].set_xlim((0, t_max))
axs[2].set_ylim(ymin=0)
"""
fig.legend(loc="lower left")

fig.tight_layout(rect=(0, 0.15, 1, 1))
#fig.savefig("../21_01_29_mini_presentation/MCE_vs_MTE_vs_exakt_omega=%g, gamma=%g.png"%(omega, gamma), dpi=300)
plt.show()
