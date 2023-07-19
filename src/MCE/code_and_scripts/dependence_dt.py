import re
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from uncertainties import ufloat

import data_manager

def calc_errors(dts, omega, gamma, initial_state, psi_initial="local",
                num_sites=2, num_electrons=1, num_samples=0):
	
	show_fit_plots = False
	print_omega_tmax = False
	err_wg = np.empty((omega.size, gamma.size))

	for w in range(omega.size):
		#lowest_freq = hopping if omega[w] == 0 else min(hopping, omega[w])
		t_max = 20.*np.pi  # 4.*np.pi/lowest_freq
		if print_omega_tmax:
			print(omega[w], t_max)
		for g in range(gamma.size):
			err = np.empty(dts.size)
			
			# get multiple samples of the time
			for i in range(dts.size):
				sim = data_manager.Simulation(omega[w], gamma[g], t_max, dts[i],
								              num_samples=num_samples, num_sites=num_sites,
								              symmetric_y_init=True,
								              num_electrons=num_electrons,
								              printProgress=False,
								              initial_state=initial_state,
								              psi_initial=psi_initial,
								              calcYSquared=True,
								              fullStorage=False,
								              hopping=1.)

				energy_total = sim.energy_total()
				energy_total -= energy_total[0]
				err[i] = np.max(np.abs(energy_total))
				
			curv_opt, curv_var = curve_fit(lambda x, a: a*x**2, dts, err, p0=err[0]/dts[0]**2)
			if show_fit_plots:
				curv = ufloat(curv_opt, np.sqrt(curv_var))
				dts_fine = np.linspace(0, dts.max(), 100)
				plt.figure(figsize=(4, 3))
				#plt.title(r"$\omega=%g, \gamma=%g$"%(omega[w], gamma[g]))
				plt.plot(dts_fine, curv_opt*dts_fine**2, '-', color="tab:orange",
				         label="${:.1uSL}\\cdot (\\Delta t \\cdot t_0)^2$".format(curv))
				plt.plot(dts, err, '.', label="simulation")  #, color="tab:blue")
				plt.xlabel("$\\Delta t\\cdot t_0$")
				plt.ylabel("$\\sigma(\\omega=%gt_0, \\gamma=%gt_0, \\Delta t)$"%(omega[w], gamma[g]))
				plt.xlim(xmin=0)
				plt.ylim(ymin=0)
				plt.legend()
				plt.tight_layout()
				plt.savefig("../report_research_and_dev/figures/error_dt_parabola.pdf")
				plt.show()
			
			err_wg[w, g] = curv_opt
		
		if not print_omega_tmax:
			print("%i/%i"%(w+1, omega.size))
	
	return err_wg



# BCDW is high around omega = 2, DCDW is weird

ref_scale = 0.2

omega = np.linspace(0.1, 3., 20)
gamma = np.linspace(0.1, 4., 20)

dts = np.linspace(0.05, 0., 10, endpoint=False)

err_wg_BCDW = calc_errors(dts, omega, gamma, "BCDW")
err_wg_DCDW = calc_errors(dts, omega, gamma, "DCDW")

cb_label = "$\\sigma(\\omega, \\gamma, \\Delta t)/(\\Delta t\\cdot t_0)^2$"
cb_label = "$\\sigma/(\\Delta t\\cdot t_0)^2$"

fig, axs = plt.subplots(1, 2, figsize=(7., 3))
axs[0].set_title("BCDW")
c = axs[0].pcolormesh(omega, gamma, err_wg_BCDW.T, shading="nearest", vmin=0)
cb = fig.colorbar(c, ax=axs[0])
cb.set_label(cb_label)
axs[0].set_xlabel(r"$\omega/t_0$")
axs[0].set_ylabel(r"$\gamma/t_0$")

axs[1].set_title("DCDW")
c = axs[1].pcolormesh(omega, gamma, err_wg_DCDW.T, shading="nearest", vmin=0)
cb = fig.colorbar(c, ax=axs[1])
cb.set_label(cb_label)
axs[1].set_xlabel(r"$\omega/t_0$")
axs[1].set_ylabel(r"$\gamma/t_0$")
fig.tight_layout(pad=0)
fig.savefig("../report_research_and_dev/figures/error_dt_omega_gamma.png", dpi=400)

quit()

normalized_err = err_wg_BCDW/gamma.reshape(1, -1)**2/omega.reshape(-1, 1)
print("(%g+-%g)*omega*gamma**2*dt**2"%(normalized_err.mean(), normalized_err.std()))

fig, axs = plt.subplots(1, 2)
axs[0].set_title("raw")
c = axs[0].pcolormesh(omega, gamma, err_wg_BCDW.T, shading="nearest", vmin=0)
fig.colorbar(c, ax=axs[0])
axs[0].set_xlabel(r"$\omega/t_0$")
axs[0].set_ylabel(r"$\gamma/t_0$")

axs[1].set_title("divided by gamma^2 and omega")
c = axs[1].pcolormesh(omega, gamma, normalized_err.T, shading="nearest", vmin=0)
cb = fig.colorbar(c, ax=axs[1])
cb.set_label(cb_label)
axs[1].set_xlabel(r"$\omega/t_0$")
axs[1].set_ylabel(r"$\gamma/t_0$")


fig, axs = plt.subplots(1, 2)
axs[0].set_title("raw")
axs[0].plot(omega, err_wg_BCDW, '.-')
axs[0].set_xlabel(r"$\omega/t_0$")

axs[1].set_title("div by gamma**2")
axs[1].plot(omega, err_wg_BCDW/gamma.reshape(1, -1)**2, '.-')
w = np.linspace(0, omega.max())
axs[1].plot(w, ref_scale*w, 'k')
axs[1].set_xlabel(r"$\omega/t_0$")


fig, axs = plt.subplots(1, 2)
axs[0].set_title("raw")
axs[0].plot(gamma, err_wg_BCDW.T, '.-')
axs[0].set_xlabel(r"$\gamma$")

axs[1].set_title("div by omega")
axs[1].plot(gamma, err_wg_BCDW.T/omega.reshape(1, -1), '.-')
g = np.linspace(0, gamma.max())
axs[1].plot(g, ref_scale*g**2, 'k')
axs[1].set_xlabel(r"$\gamma$")




#plt.legend()
plt.show()

