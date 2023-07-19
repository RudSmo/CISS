import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import data_manager


num_sites = 2 # only 2 is supported
omega, gamma = .1, .4

dt = 0.025
t_max = 20.

num_samples = 100
num_bins = 50

save_video = False

initial_state="DCDW"
psi_initial = "local"

ref = holstein.RefData(initial_state, psi_initial, num_sites, omega, gamma)
ref.set_cut_time(t_max)

sim = data_manager.Simulation(omega, gamma, t_max, dt,
	                          num_samples, num_sites=num_sites,
	                          symmetric_y_init=True, y_initial=None,
	                          num_electrons=1, initial_state=initial_state,
	                          psi_initial=psi_initial,
	                          fullStorage=True, printProgress=True)
"""
q_centered = sim.get_q_full() - sim.get_q_mean()
q_centered_sq = q_centered**2
q_var = q_centered_sq.sum(axis=0)/(q_centered.shape[0]-1)
print(q_centered.shape, q_centered_sq/(2*  q_var.reshape(1, -1)).mean(axis=0))
print(q_var.shape)

p_q_centered = sim.get_p_q_full() - sim.get_p_q_mean()
p_q_centered_sq = p_q_centered**2
p_q_var = (p_q_centered_sq).sum(axis=0)/(p_q_centered.shape[0]-1)

kl_div =   (  q_centered_sq/(2*  q_var.reshape(1, -1))).mean(axis=0) \
         #+ (p_q_centered_sq/(2*p_q_var.reshape(1, -1))).mean(axis=0) \
         #+ .5*np.log(q_var*(2*np.pi)) + .5*np.log(p_q_var*(2*np.pi)) \
         #- .5*np.log(q_var*(2*np.pi*np.e)) - .5*np.log(p_q_var*(2*np.pi*np.e))
#kl_div -= kl_div[0]
plt.plot(sim.time(), kl_div)

plt.figure()
plt.plot(sim.time(), .5*np.log(q_var*(2*np.pi)) + .5*np.log(p_q_var*(2*np.pi)))


plt.figure()
plt.pcolormesh(q_centered)
plt.colorbar()

plt.figure()
plt.hist(q_centered[:, -10], bins=num_bins)

plt.show()
quit()
"""

hist_data, hist_edges = sim.get_histogram_xp(num_bins)
q_hist_data, q_hist_edges = sim.get_histogram_qpq(num_bins)



fig, axs = plt.subplots(2,3, figsize=(15,10))
axs = axs.flatten()
fig.suptitle("$\\omega$=%1.1f, $\\gamma$=%1.1f, t=%.3f"%(omega, gamma, sim.time()[0]))

# time evolution of mean
axs[0].plot(sim.x(0), sim.p(0), "r")
axs[1].plot(sim.x(1), sim.p(1), "r")
axs[2].plot(sim.q(), sim.pq(), "r")
axs[0].plot(ref.x(0), ref.p(0), "g--")
axs[1].plot(ref.x(1), ref.p(1), "g--")
axs[2].plot(ref.x(0)-ref.x(1), ref.p(0)-ref.p(1), "g--")

# indicator of current mean
mid0, = axs[0].plot(sim.x(0)[0], sim.p(0)[0], ".r")
mid1, = axs[1].plot(sim.x(1)[0], sim.p(1)[0], ".r")
midq, = axs[2].plot(sim.q()[0], sim.pq()[0], ".r")

# full colorplot of distribution
c0 = axs[0].pcolormesh(hist_edges[0], hist_edges[1], hist_data[0, 0, :, :].T, cmap="Greys")
c1 = axs[1].pcolormesh(hist_edges[0], hist_edges[1], hist_data[0, 1, :, :].T, cmap="Greys")
cq = axs[2].pcolormesh(q_hist_edges[0], q_hist_edges[1], q_hist_data[0, :, :].T, cmap="Greys")


hist_x_line, = axs[3].plot(hist_edges[0, :-1], hist_data[0, 0, :, :].sum(axis=1), label="x0")
hist_p_line, = axs[3].plot(hist_edges[1, :-1], hist_data[0, 0, :, :].sum(axis=0), label="p0")

hist_q_line,  = axs[4].plot(q_hist_edges[0, :-1], q_hist_data[0, :, :].sum(axis=1), label="q")
hist_pq_line, = axs[4].plot(q_hist_edges[1, :-1], q_hist_data[0, :, :].sum(axis=0), label="pq")

axs[5].plot(sim.time(), sim.get_q_kurtosis(), color="tab:blue", label="q")
midq_kurt, = axs[5].plot(sim.time()[0], sim.get_q_kurtosis()[0], ".", color="tab:blue")
axs[5].plot(sim.time(), sim.get_p_q_kurtosis(), color="tab:orange", label="pq")
midpq_kurt, = axs[5].plot(sim.time()[0], sim.get_p_q_kurtosis()[0], ".", color="tab:orange")

# some labels
axs[0].set_title("site 0")
axs[0].set_xlabel("x0/l0")
axs[0].set_ylabel("p0*l0")

axs[1].set_title("site 1")
axs[1].set_xlabel("x1/l0")
axs[1].set_ylabel("p1*l0")

axs[2].set_title("q=x0-x1 - basis")
axs[2].set_xlabel("q/l0")
axs[2].set_ylabel("p_q*l0")


axs[3].set_xlabel("x0 or p0")
axs[3].set_ylabel("P(x0 or p0)")
axs[3].axes.legend()

axs[4].set_xlabel("q or p_q")
axs[4].set_ylabel("P(q or p_q)")
axs[4].axes.legend()

axs[5].set_xlabel("t")
axs[5].set_ylabel("excess kurtosis of q or p_q")
axs[5].axes.legend()

fig.tight_layout()

def update(frame):
	fig.suptitle("$\\omega$=%1.1f, $\\gamma$=%1.1f, t=%.3f"%(omega, gamma, sim.time()[frame]))
	c0.set_array(hist_data[frame, 0, :, :].T.ravel())
	c1.set_array(hist_data[frame, 1, :, :].T.ravel())
	cq.set_array(q_hist_data[frame, :, :].T.ravel())
	
	mid0.set_data(sim.x(0)[frame], sim.p(0)[frame])
	mid1.set_data(sim.x(1)[frame], sim.p(1)[frame])
	midq.set_data(sim.q()[frame], sim.pq()[frame])
	
	hist_x_line.set_data(hist_edges[0, :-1], hist_data[frame, 0, :, :].sum(axis=1))
	hist_p_line.set_data(hist_edges[1, :-1], hist_data[frame, 0, :, :].sum(axis=0))
	
	hist_q_line.set_data(q_hist_edges[0, :-1], q_hist_data[frame, :, :].sum(axis=1))
	hist_pq_line.set_data(q_hist_edges[1, :-1], q_hist_data[frame, :, :].sum(axis=0))
	
	midq_kurt.set_data(sim.time()[frame], sim.get_q_kurtosis()[frame])
	midpq_kurt.set_data(sim.time()[frame], sim.get_p_q_kurtosis()[frame])

shown_times = np.arange(sim.time().size)
time_stepping = max(1, sim.time().size//160)
shown_times = shown_times[::time_stepping]

ani = FuncAnimation(fig, update, frames=np.arange(shown_times.size), blit=False,
                    repeat=True, repeat_delay=1000, interval=10000//shown_times.size)

if save_video:
	print("saving to file...")
	ani.save(file_prefix+"holstein_dimer_%s_phase_space_hist_omega=%1.1f_gamma=%1.2f_nosym4e4.mp4"%(initial_state, omega, gamma))

plt.show()

