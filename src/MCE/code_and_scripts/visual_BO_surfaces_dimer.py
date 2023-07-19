import numpy as np
import matplotlib.pyplot as plt
import data_manager

# latex style
use_tex = True
if use_tex:
	from matplotlib import rc
	rc('text', usetex=use_tex)
	rc('font', family='serif')

hopping = 2.5
gamma   = 2.5

def BO_surfaces(q):
	half_q2 = .5*q**2
	e_diff = np.sqrt(q**2*gamma**2 + hopping**2)
	return np.stack((half_q2 - e_diff, half_q2+e_diff))

q_range  = 4

label_q  = r"$q/\sqrt{\frac{\hbar}{m\omega_0}}$"
label_pq = r"$p_q/\sqrt{\hbar m\omega_0}$"

q_grid = np.linspace(-q_range, q_range, num=100)
energy_q = BO_surfaces(q_grid)

plt.figure(figsize=(.7*6, 3))
plt.plot(q_grid, energy_q[0], label="$E_{-}^{BO}(q)$")
plt.plot(q_grid, energy_q[1], '--', label="$E_{+}^{BO}(q)$")

plt.legend()
plt.xlabel(label_q)
plt.ylabel(r"$E^{BO}/\hbar\omega_0$")

plt.tight_layout(pad=0.7)
plt.savefig("../thesis/figures/BO_surfaces_dimer_t0=%g_omega=1_gamma=%g.pdf"%(hopping, gamma))
plt.show()
