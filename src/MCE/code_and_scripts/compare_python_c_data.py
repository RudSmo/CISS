import numpy as np
import matplotlib.pyplot as plt

from ehrenfest_holstein_1d_c_sim import simulate_system
from ehrenfest_holstein_1d_py_sim import simulate_system_py, generate_initial_classical_state

hopping = 1.
omega, gamma = .5, 1.
mass_ho = 1/omega 	# mass of particles in harmonic oscillator
l0 = np.sqrt(1/(mass_ho*omega))
# number of ehrenfest trajectories (randomly chosen)
# use num_samples = 0 for single trajectory starting in rest on all sites
num_samples = 0
# generates 3*num_samples to enforce symmetry
symmetrize_start = True

# dt is the size of an entire step, NOT the half step as in the notes
dt = 0.025
t_max = 4 #10/gamma
t, dt = np.linspace(0, t_max, int(t_max/dt)+1, retstep=True)


# initial electron state
num_sites = 2
num_electrons = 1
psi_initial = np.zeros((num_electrons, num_sites), dtype=np.complex)
for i in range(num_electrons):
	psi_initial[i, 2*i+(num_sites+1)//2-num_electrons] = 1



y = generate_initial_classical_state(num_sites, num_samples, symmetrize_start, mass_ho, omega)

# comparison between C and Python version
y_all_times_py, density_py, Nph_py = simulate_system_py(y.copy(), psi_initial.copy(),
		                                                hopping, omega, mass_ho, gamma, t)	
y_all_times, density, Nph, y_squared = simulate_system(y.copy(), psi_initial.copy(),
		                                    hopping, omega, mass_ho, gamma, t, True, True)
print(Nph[-1])
print(Nph_py[-1])

plt.figure()
plt.plot(y_all_times**2)

plt.figure()
plt.plot(y_squared)

plt.figure()
plt.title("difference y_c-y_py")
diff_y = y_all_times-y_all_times_py
amp = max(diff_y.max(), abs(diff_y.min()))
plt.pcolormesh(diff_y, vmin=-amp, vmax=amp, cmap="bwr")
plt.colorbar()

plt.figure()
plt.title("diff density d_c-d_py")
diff_d = density-density_py
amp = max(diff_d.max(), abs(diff_d.min()))
plt.pcolormesh(diff_d, vmin=-amp, vmax=amp, cmap="bwr")
plt.colorbar()

plt.figure()
plt.title("diff density d_c-d_py")
plt.plot(np.abs(diff_d[:, 0]), '.-')
plt.plot(np.abs(diff_d[:, 1]), '.-')
plt.yscale("log")

plt.figure()
plt.title("densities d_c and d_py")
plt.plot(density[:, 0], label="d_c")
plt.plot(density_py[:, 0], label="d_py")

plt.figure()
plt.title("diff phonons Nph_c-Nph_py")
diff_Nph = Nph-Nph_py
amp = max(diff_Nph.max(), abs(diff_Nph.min()))
plt.pcolormesh(diff_Nph, vmin=-amp, vmax=amp, cmap="bwr")
plt.colorbar()

v_min = 0
v_max = 1

plt.figure()
plt.title("density d_c")
plt.pcolormesh(density, vmin=v_min, vmax=v_max)
plt.colorbar()

plt.figure()
plt.title("density d_py")
plt.pcolormesh(density_py, vmin=v_min, vmax=v_max)
plt.colorbar()


plt.figure()
plt.title("density Nph_c")
plt.pcolormesh(Nph, vmin=0)
plt.colorbar()

plt.figure()
plt.title("density Nph_py")
plt.pcolormesh(Nph_py, vmin=0)
plt.colorbar()

plt.show()
