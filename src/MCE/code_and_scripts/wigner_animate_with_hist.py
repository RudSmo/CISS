import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import ListedColormap
from multiprocessing import Process

import data_manager

# TODO this only works if sim.dt*4 = wigner.dt
def update(frame, fig, sim, wigner, cqpqE, cqpqW, midqE_0, midqW_0, midqE_1, midqW_1):
	time_scale = 4
	fig.suptitle("%s %s: $\\omega$=%1.1f, $\\gamma$=%1.1f, t=%.3f"\
	             %(initial_state, psi_initial, omega, gamma, sim.time()[time_scale*frame]))
	wigner.load(frame)	
	cqpqE.set_array(q_hist_data[time_scale*frame, :, :].T.ravel())
	cqpqW.set_array(wigner.wigner.ravel())
	
	midqE_0.set_data(sim.q()[time_scale*frame], sim.pq()[time_scale*frame])
	midqE_1.set_data(sim.q()[time_scale*frame], sim.pq()[time_scale*frame])
	midqW_0.set_data(wigner.q1[frame], wigner.pq1[frame])
	midqW_1.set_data(wigner.q1[frame], wigner.pq1[frame])

def draw_frames(shown_times, time_indices):
	#global sim, wigner
	fig, axs = plt.subplots(1,2, figsize=(16, 9))
	axs = axs.flatten()
	fig.suptitle("%s %s: $\\omega$=%1.1f, $\\gamma$=%1.1f, t=%.3f"%(initial_state, psi_initial, omega, gamma, sim.time()[0]))

	axs[0].set_title("MT Ehrenfest")
	axs[0].set_xlabel("q")
	axs[0].set_ylabel("p_q")

	axs[1].set_title("Wigner (Michael)")
	axs[1].set_xlabel("q")
	axs[1].set_ylabel("p_q")


	if custom_ax_range is not None:
		axs[0].set_xlim(custom_ax_range[0])
		axs[0].set_ylim(custom_ax_range[1])
		axs[1].set_xlim(custom_ax_range[0])
		axs[1].set_ylim(custom_ax_range[1])
	elif use_wigner_ax_lim:
		axs[0].set_xlim((wigner.q_sampled[0], wigner.q_sampled[-1]))
		axs[0].set_ylim((wigner.p_q_sampled[0], wigner.p_q_sampled[-1]))
		axs[1].set_xlim((wigner.q_sampled[0], wigner.q_sampled[-1]))
		axs[1].set_ylim((wigner.p_q_sampled[0], wigner.p_q_sampled[-1]))
	else:
		axs[0].set_xlim((q_hist_edges[0,0], q_hist_edges[0,-1]))
		axs[0].set_ylim((q_hist_edges[1,0], q_hist_edges[1,-1]))
		axs[1].set_xlim((q_hist_edges[0,0], q_hist_edges[0,-1]))
		axs[1].set_ylim((q_hist_edges[1,0], q_hist_edges[1,-1]))
	
	cqpqW = axs[1].pcolormesh(wigner.q_sampled, wigner.p_q_sampled, wigner.wigner, shading="gouraud",
		                      cmap=newcmp, vmin=min_density, vmax=max_density)
	cqpqE = axs[0].pcolormesh(q_hist_edges[0], q_hist_edges[1], q_hist_data[0, :, :].T,
		                      cmap=newcmp, vmin=min_density, vmax=max_density)

	#fig.colorbar(cqpqW, ax=axs[2])
	fig.subplots_adjust(left=0.05, bottom=.1, right=0.93, top=.9, wspace=0.1)
	cbar_ax = fig.add_axes([0.94, 0.1, 0.02, 0.8])
	fig.colorbar(cqpqW, cax=cbar_ax)

	# time evolution of mean
	axs[0].plot(sim.q(), sim.pq(), "r-")
	axs[0].plot(wigner.q1, wigner.pq1, "c--")
	axs[1].plot(sim.q(), sim.pq(), "r-")
	axs[1].plot(wigner.q1, wigner.pq1, "c--")

	# indicator of current mean
	midqW_0, = axs[0].plot(wigner.q1[0], wigner.pq1[0], ".c", markersize=8)
	midqW_1, = axs[1].plot(wigner.q1[0], wigner.pq1[0], ".c", markersize=8)
	midqE_0, = axs[0].plot(sim.q()[0], sim.pq()[0], ".r", markersize=8)
	midqE_1, = axs[1].plot(sim.q()[0], sim.pq()[0], ".r", markersize=8)
	
	for i in time_indices:
		frame = shown_times[i]
		update(frame, fig, sim, wigner, cqpqE, cqpqW, midqE_0, midqW_0, midqE_1, midqW_1)
		fig.savefig("images_for_video/pic%i.png"%(i), dpi=120)
		print("%i/%i"%(i+1, shown_times.size))


if __name__ == "__main__":
	initial_state="BCDW"
	psi_initial = "local"
	omega, gamma = 2., 4.

	dt = 0.025
	t_max = 10.

	use_wigner_ax_lim = True
	custom_ax_range = ((-8, 8), (-4, 4))
	#custom_ax_range = None

	num_samples = 40000
	num_bins = 100

	save_video = True
	time_stepping = 1
	num_processes = 4
	file_name = "Wigner_vs_MTEhrenfest_holstein_dimer_%s_%s_omega=%g_gamma=%g_sym%i"\
		        %(initial_state, psi_initial, omega, gamma, num_samples)

	wigner = data_manager.Wigner(omega, gamma, initial_state, psi_initial)


	sim = data_manager.Simulation(omega, gamma, t_max,
			                      num_samples=num_samples, num_sites=2,
			                      symmetric_y_init=False, y_initial=None,
			                      num_electrons=1,
			                      initial_state=initial_state,
		                          psi_initial=psi_initial,
			                      fullStorage=True, printProgress=True)
	
	shown_times = np.arange(wigner.time.size)
	shown_times = shown_times[::time_stepping]
	
	q_hist_data, q_hist_edges = sim.get_histogram_qpq(num_bins)

	sim.y_full = None
	sim.q_full = None

	plt.figure()
	plt.title("%s %s: $\\omega=%g$, $\\gamma=%g$"%(initial_state, psi_initial, omega, gamma))
	plt.xlabel("time")
	plt.ylabel(r"$\langle n \rangle$")
	plt.tight_layout()
	plt.plot(sim.time(), sim.nel())
	try:
		plt.gca().set_prop_cycle(None)
		plt.plot(wigner.time, wigner.nel(), '--')
	except FileNotFoundError:
		pass
	if save_video:
		plt.savefig(wigner.path + file_name + "_electron_density.png")


	# set up the custom color map
	max_density = 1.1*wigner.wigner.max()
	min_density = -max_density/4.

	#cm = plt.get_cmap("viridis")
	old_cm = plt.get_cmap('RdGy', 8*256)
	newcolors = old_cm(np.linspace(3/8, 1, 5*256))
	newcolors[:256, :] = old_cm(np.linspace(0, 1/2, 256))
	newcmp = ListedColormap(newcolors)


	if save_video:
		video_name = wigner.path+file_name+".mp4"
		os.system("rm -f images_for_video/*")
		os.system("mkdir images_for_video")
		print("creating images for video:")
		
		# create all images
		time_indices = np.arange(shown_times.size)
		process = []
		for i in range(num_processes):
			p = Process(target=draw_frames, args=(shown_times, time_indices[i::num_processes],))
			p.start()
			process.append(p)
		
		for p in process:
			p.join()
		
		# create the video
		os.system("ffmpeg -y -r %g -s 1920x1080 -i images_for_video/pic%%d.png -vcodec libx264 -crf 20 %s"%(10/time_stepping, video_name))
	else:
		ani = FuncAnimation(fig, update, frames=shown_times, blit=False,
				            repeat=True, repeat_delay=0, interval=int(t_max*1000/(shown_times.size-1)))

		plt.show()
