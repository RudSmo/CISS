import numpy as np
import matplotlib.pyplot as plt

import data_manager

from src.ehrenfest_init import initial_state_MCE
from simulate import simulate_system_MCE

import os
from multiprocessing import Process

initial_state="BCDW"
psi_initial  ="local"

omega = 1
gamma = 1
hopping= 1
t_max = 1.
video_duration = t_max
last_frame_hold_duration = 0.5
save_plots = False
save_data = False

save_init_and_quit = True

num_sites_x = 5
num_sites_y = 5 # num_sites_x
figsize = (16, 4)
num_electrons = 1 # either 1 or negative (cdw) if num_sites_y > 1, negativ will fail in MCE

num_samples = -1000  # <0 to not calc MTE

num_configs = 2000
max_num_configs = num_configs
dt = 0.01

dt_mte = min(0.05, np.inf if omega*gamma == 0. else np.sqrt(10*1e-6/omega/gamma**2))
mte_time_stepping = int(np.ceil(dt/dt_mte))
dt_mte = dt/mte_time_stepping

if num_samples >= 0:
	sim = data_manager.Simulation(omega, gamma, t_max, time_step=dt_mte,
	                              num_samples=num_samples,
	                              num_sites_x=num_sites_x,
	                              num_sites_y=num_sites_y,
	                              symmetric_y_init=True,
	                              num_electrons=num_electrons,
	                              printProgress=True,
	                              initial_state=initial_state,
	                              psi_initial=psi_initial,
	                              calcYSquared=True, hopping=hopping)
else:
	sim = None


if num_sites_x == 6 and num_sites_y == 6:
	ref_path = "ref_data/DMRG2/2dholtdvplbo_single_ebcdwLx6Ly6M30site_in16Md3000lboMd100t01.0gam-2.0omg1.0dt0.01tot20.0lbocut1e-09cut1e-09DIR/"
elif num_sites_x == 5 and num_sites_y == 5:
	ref_path = "ref_data/DMRG2/2dholtdvplbo_single_ebcdwLx5Ly5M20site_in8Md3000lboMd100t01.0gam-2.0omg1.0dt0.01tot20.0lbocut1e-09cut1e-09DIR/"
elif num_sites_x == 11 and num_sites_y == 1:
	ref_path = "ref_data/DMRG2/holtdmrglbo_single_ebcdwL11M30site_in6Md3000lboMd100t01.0gam-2.0omg1.0dt0.01tot20.0lbocut1e-09cut1e-09DIR/"
else:
	ref_path = None

if ref_path is not None:
	ref_time = np.loadtxt(ref_path + "time.dat")
	ref_nel = np.empty((ref_time.size, num_sites_x*num_sites_y))
	ref_nph = np.empty((ref_time.size, num_sites_x*num_sites_y))
	for i in range(num_sites_x*num_sites_y):
		ref_nel[:, i] = np.loadtxt(ref_path + "nlocSite%i.dat"%(i+1))
		ref_nph[:, i] = np.loadtxt(ref_path + "NphlocSite%i.dat"%(i+1))
else:
	ref_time = None
	ref_nel = None
	ref_nph = None


num_times = int(np.round(t_max/dt)) + 1  # +1 for t=0

hoz, el_coeff = initial_state_MCE(num_sites_x, num_sites_y, num_configs, max_num_configs, num_electrons, 
                                  False, initial_state, psi_initial,
                                  1, omega, gamma,
                                  False, 0.,
                                  plotOverlapMatrix=False, use_QSS=True, random_hoz=True,
                                  useDisorderRegularization=True)

if save_init_and_quit:
	print("saving data of size %i Bytes"%(2*num_sites_x*num_sites_y*max_num_configs*16))
	file_raw_start = open("data_export/raw_start.dat", "wb")
	hoz.tofile(file_raw_start)
	el_coeff.tofile(file_raw_start)
	file_raw_start.close()
	quit()


obs_y, obs_nph, obs_nel, obs_el_ekin, obs_el_ph_inter = simulate_system_MCE(hoz, el_coeff, hopping, omega, gamma, num_times, dt, num_configs, printProgress=True)


time = np.linspace(0, t_max, num_times)

plt.figure()
e_tot =   omega * obs_nph.sum(axis=(1, 2))\
        + gamma * obs_el_ph_inter.sum(axis=(1, 2))\
        + hopping * obs_el_ekin.sum(axis=1)
plt.title(r"MCE: $E_{tot}(t)-(E_{tot}(0)=%g)$ "%e_tot[0])
e_tot -= e_tot[0]
plt.plot(time, e_tot)

title = r"%s, $\omega=%g, \gamma=%g$: "%("bare" if initial_state == "BCDW" else "dressed", omega, gamma)


fig_nel, axs = plt.subplots(num_sites_y, num_sites_x, figsize=figsize, squeeze=False)
fig_nel.suptitle(title + "nel")
for x in range(num_sites_x):
	for y in range(num_sites_y):
		axs[y, x].plot(time, obs_nel[:, y, x])
		if sim is not None:
			axs[y, x].plot(sim.time(), sim.nel(x, y), '--')
		if ref_time is not None:
			axs[y, x].plot(ref_time, ref_nel[:, y*num_sites_x+x], "--")
		axs[y, x].set_xlim((0, t_max))
fig_nel.tight_layout()

fig_nph, axs = plt.subplots(num_sites_y, num_sites_x, figsize=figsize, squeeze=False)
fig_nph.suptitle(title + "nph")
for x in range(num_sites_x):
	for y in range(num_sites_y):
		axs[y, x].plot(time, obs_nph[:, y, x],)
		if sim is not None:
			axs[y, x].plot(sim.time(), sim.nph(x, y), '--')
		if ref_time is not None:
			axs[y, x].plot(ref_time, ref_nph[:, y*num_sites_x+x], "--")
		axs[y, x].set_xlim((0, t_max))
fig_nph.tight_layout()


fig_x, axs = plt.subplots(num_sites_y, num_sites_x, figsize=figsize, squeeze=False)
fig_x.suptitle(title + "x")
for x in range(num_sites_x):
	for y in range(num_sites_y):
		axs[y, x].plot(time, obs_y[:, y, x, 0])
		if sim is not None:
			axs[y, x].plot(sim.time(), sim.x(x, y), '--')
fig_x.tight_layout()

fig_p, axs = plt.subplots(num_sites_y, num_sites_x, figsize=figsize, squeeze=False)
fig_p.suptitle(title + "p")
for x in range(num_sites_x):
	for y in range(num_sites_y):
		axs[y, x].plot(time, obs_y[:, y, x, 1])
		if sim is not None:
			axs[y, x].plot(sim.time(), sim.p(x, y), '--')
		axs[y, x].set_xlim((0, t_max))
fig_p.tight_layout()


plt.figure()
plt.title("ekin")
plt.plot(time, obs_el_ekin.sum(axis=1)*1)
if sim is not None:
	plt.plot(sim.time(), sim.energy_kinetic()*1, '--')
plt.xlim((0, t_max))


def update(frame, fig, sim_mte, mte_nel, obs_nel, mce_nel, mte_time_stepping):
	if sim_mte is not None:
		mte_nel.set_array(sim_mte.nel()[mte_time_stepping*frame].ravel())
	mce_nel.set_array(obs_nel[frame].ravel())

def draw_frames(shown_times, time_indices, sim_mte, obs_nel, mte_time_stepping, time):
	#global sim, wigner
	fig, axs = plt.subplots(1,2, figsize=(16, 9))
	axs = axs.flatten()
	fig.suptitle("%s %s: $\\omega$=%1.1f, $\\gamma$=%1.1f, t=%.3f"%(initial_state, psi_initial, omega, gamma, time[0]))
	
	axs[0].set_title("MTE")
	
	axs[1].set_title("MCE")
	
	if sim_mte is not None:
		mean_max = sim_mte.nel(squeezeY=False).max(axis=(1, 2)).mean()
		mte_nel = axs[0].pcolormesh(sim_mte.nel()[0], # shading="nearest",
		                            cmap="viridis", vmin=0, vmax=mean_max)
	else:
		mte_nel = None
	
	mean_max = obs_nel.max(axis=(1, 2)).mean()
	mce_nel  = axs[1].pcolormesh(obs_nel[0], # shading="nearest",
	                            cmap="viridis", vmin=0, vmax=mean_max)
	
	fig.subplots_adjust(left=0.05, bottom=.1, right=0.93, top=.9, wspace=0.1)
	cbar_ax = fig.add_axes([0.94, 0.1, 0.02, 0.8])
	fig.colorbar(mce_nel, cax=cbar_ax)
	
	for i in time_indices:
		frame = shown_times[i]
		fig.suptitle("%s %s: $\\omega$=%g, $\\gamma$=%g, t=%.3f"\
		             %(initial_state, psi_initial, omega, gamma, time[frame]))
		update(frame, fig, sim_mte, mte_nel, obs_nel, mce_nel, mte_time_stepping)
		fig.savefig("images_for_video/pic%i.png"%(i), dpi=120)
		print("%i/%i"%(i+1, shown_times.size))

if save_data or save_plots:
	path = "data_export/" + "%s_psi_initial=%s_omega=%g_gamma=%g_mce_numconf=%i-%i_num_sites_x=%i_num_sites_y=%i/"\
	       %(initial_state, psi_initial, omega, gamma, num_configs, max_num_configs, num_sites_x, num_sites_y)
	os.system("mkdir -p "+path)

if save_data:
	np.save(path+"time", time[:])
	np.save(path+"x", obs_y[:, :, :, 0])
	np.save(path+"p", obs_y[:, :, :, 1])
	np.save(path+"electron_density", obs_nel[:, :, :])  # , header="site0\tsite1")
	np.save(path+"energy_total", e_tot[:])  #, header="total energy difference to start (should be 0)")
	np.save(path+"Nph", obs_nph[:, :, :])
	
	if sim is not None:
		np.save(path+"time_MTE", sim.time())
		np.save(path+"nel_MTE", sim.nel())

if save_plots:
	shown_times = np.arange(time.size)
	path = "video_2d_spreading/omega=%g_gamma=%g_%ix%i_samples=%i_configs=%i/"\
	       %(omega, gamma, num_sites_x, num_sites_y, num_samples, num_configs)
	os.system("mkdir -p " + path)
	video_name = path+"nel.mp4"
	
	fig_nel.savefig(path + "nel.png", dpi=120)
	fig_nph.savefig(path + "nph.png", dpi=120)
	fig_x.savefig(path + "x.png", dpi=120)
	fig_p.savefig(path + "p.png", dpi=120)
	
	os.system("rm -f images_for_video/*")
	os.system("mkdir -p images_for_video")
	print("creating images for video:")

	# create all images
	time_indices = np.arange(shown_times.size)
	process = []
	num_processes = 4
	for i in range(num_processes):
		p = Process(target=draw_frames, args=(shown_times, time_indices[i::num_processes], sim, obs_nel, mte_time_stepping, time))
		p.start()
		process.append(p)

	for p in process:
		p.join()
	
	# create the video
	os.system("ffmpeg -y -r %g -s 1920x1080 "
	          "-i images_for_video/pic%%d.png -vcodec libx265 -crf 20 "
	          "-vf tpad=stop_mode=clone:stop_duration=%g %s"\
	          %(time.size/video_duration, last_frame_hold_duration, video_name))  #sim.time().size


plt.show()
