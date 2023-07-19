import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os

from src.ehrenfest_init import initial_state_MCE
from simulate import simulate_system_MCE

import data_manager
from multiprocessing import Process

original_backend = matplotlib.rcParams['backend']

version = 2
hopping = 10
omega = 1
gamma = 4
num_sites_x = 7
num_sites_y = 1
plot_q = num_sites_y == 1 and num_sites_x == 2 and True

num_configs = 200
max_num_configs = 600

time_step = 0.01
time_max = 5

initial_state = "BCDW"
psi_initial = "local"
psi_state = None
"""
vmax = None #0.2
xlim = (-3, 7)
plim = (-5, 5)
"""
vmax = None #min(1., 20/num_configs)
markersize = 9
color = "BO_state" # index, BO_state or coeff

if color not in ["index", "BO_state", "coeff"]:
	raise ValueError(f"invalid coloring scheme \"{color}\"")

print("running simulation...")
sim = data_manager.Simulation_MCE(version, hopping, omega, gamma,
                                  num_sites_x, num_sites_y,
                                  time_max, time_step,
                                  num_configs,
                                  initial_state, psi_initial, psi_state,
                                  max_num_configs, fullStorage=True)

title = sim.get_sim_name()

output_path = sim.get_storage_path(override_export_path="video_config_spreading")
output_path_sh = output_path.replace("(", "\(").replace(")", "\)")

print("saving to", output_path)

video_name = output_path_sh+"spreading.mp4"
video_duration = 10
last_frame_hold_duration = 0.5

os.system("rm -f images_for_video/*")
os.system("mkdir -p images_for_video")
os.system("mkdir -p " + output_path_sh)



def draw_frames(time_indices, sim, output_path):
	matplotlib.use('Agg')
	fig = plt.figure()
	fig, axs = plt.subplots(num_sites_y, num_sites_x+int(plot_q), figsize=(16,9), squeeze=False)
	fig.suptitle(title)
	fig.tight_layout()
	xmargin = 0.1
	pmargin = 0.1
	xlim = (np.sqrt(2)*np.nanmin(sim.fullStorage_hoz.real) - xmargin, np.sqrt(2)*np.nanmax(sim.fullStorage_hoz.real) + xmargin)
	plim = np.nanmax(np.abs(np.sqrt(2)*sim.fullStorage_hoz.imag))
	plim = (-plim - pmargin, plim + pmargin)
	
	if plot_q:
		q  = sim.fullStorage_hoz[:, 0, 0, :].real - sim.fullStorage_hoz[:, 0, 1, :].real
		pq = sim.fullStorage_hoz[:, 0, 0, :].imag - sim.fullStorage_hoz[:, 0, 1, :].imag
		qlim = (np.sqrt(2)*np.nanmin(q) - xmargin, np.sqrt(2)*np.nanmax(q))
		pqlim = np.sqrt(2)*np.nanmax(np.abs(pq))
		pqlim = (-pqlim - pmargin, pqlim + pmargin)
		del q, pq
	
	
	for t in time_indices:
		# hoz and el_coeff have been changed in place
		fig.suptitle(title + " weights t=%5.2f"%(sim.time()[t]))
		
		el_coeff = sim.fullStorage_el_coeff[t, :, :, :sim.fullStorage_num_configs[t]]
		hoz = sim.fullStorage_hoz[t, :, :, :sim.fullStorage_num_configs[t]]
		
		
		if color == "BO_state":
			BO_state = sim.get_config_highest_occupied_eigenstate(t)
		
		for j in range(num_sites_y):
			for i in range(num_sites_x):
				axs[j, i].clear()
				if color == "index":
					axs[j, i].plot(np.sqrt(2)*hoz[j, i, :].real.reshape(1, -1),
					               np.sqrt(2)*hoz[j, i, :].imag.reshape(1, -1),
					               '.', ms=markersize)
				elif color == "coeff":
					order = np.argsort(np.abs(el_coeff[j, i, :]))
					cm = axs[j, i].scatter(np.sqrt(2)*hoz[j, i, order].real, np.sqrt(2)*hoz[j, i, order].imag,
					                       c=np.abs(el_coeff[j, i, order]), vmin=0, vmax=vmax, s=markersize)
				elif color == "BO_state":
					for n in range(num_sites_x*num_sites_y):
						axs[j, i].plot(np.sqrt(2)*hoz[j, i, BO_state == n].real.reshape(-1),
						               np.sqrt(2)*hoz[j, i, BO_state == n].imag.reshape(-1),
						               '.', ms=markersize)
				#axs[j, i].set_xlabel("x")
				#axs[j, i].set_ylabel("p")
				axs[j, i].set_xlim(xlim)
				axs[j, i].set_ylim(plim)
		
		if plot_q:
			axs[0, -1].clear()
			axs[0, -1].set_xlim(qlim)
			axs[0, -1].set_ylim(pqlim)
			if color == "index":
				axs[0, -1].plot(np.sqrt(2)*(hoz[0, 0, :].real-hoz[0, 1, :].real).reshape(1, -1),
				                np.sqrt(2)*(hoz[0, 0, :].imag-hoz[0, 1, :].imag).reshape(1, -1), '.', ms=markersize)
			elif color == "coeff":
				axs[0, -1].scatter(np.sqrt(2)*(hoz[0, 0, :].real-hoz[0, 1, :].real).reshape(1, -1),
				                   np.sqrt(2)*(hoz[0, 0, :].imag-hoz[0, 1, :].imag).reshape(1, -1),
				                   c=np.abs(el_coeff[0, 0, :]), vmin=0, vmax=None, s=markersize)
			elif color == "BO_state":
				for n in range(num_sites_x*num_sites_y):
					axs[0, -1].plot(np.sqrt(2)*(hoz[0, 0, BO_state == n].real-hoz[0, 1, BO_state == n].real).reshape(-1),
					                np.sqrt(2)*(hoz[0, 0, BO_state == n].imag-hoz[0, 1, BO_state == n].imag).reshape(-1), '.', ms=markersize)
		
		fig.savefig("images_for_video/pic%i.png"%(t), dpi=120)
		if time_indices[0] == 0:
			print("%i/%i"%(t+1, time_indices[-1]))
	plt.close(fig)


print("plotting...")
if np.max(sim.nel()[-1]) == 0:
	# sim aborted
	first_bad_index = np.argmax(np.max(sim.nel(squeezeY=False), axis=(1,2)) == 0)
	time_indices = np.arange(first_bad_index)
	print(f"sim aborted, {first_bad_index=}")
else:
	# normal case
	time_indices = np.arange(sim.time().size)
if False:
	draw_frames(time_indices, sim, output_path)
else:
	process = []
	num_processes = 4
	for i in range(num_processes):
		p = Process(target=draw_frames, args=(time_indices[i::num_processes], sim, output_path))
		p.start()
		process.append(p)

	for p in process:
		p.join()

matplotlib.use(original_backend)

vmin = 0
if vmax is None:
	vmax = np.nanmax(np.abs(sim.fullStorage_el_coeff))


fig, axs = plt.subplots(num_sites_y, num_sites_x, figsize=(16,9), squeeze=False)
fig.suptitle("weights $|c_{ij}|$, vmax=%g"%vmax)
for j in range(num_sites_y):
		for i in range(num_sites_x):
			axs[j,i].pcolormesh(sim.time(), np.arange(max_num_configs),
			                    np.abs(sim.fullStorage_el_coeff[:, j, i, :]).T,
			                    vmin=vmin, vmax=vmax,
			                    shading="nearest")
fig.tight_layout()
fig.savefig(output_path + "conf_weights.png", dpi=300)

e_ph = omega * sim.nph(squeezeY=False).sum(axis=(1,2))
e_kin = sim.energy_kinetic()
e_inter = sim.energy_interaction()
e_tot = e_ph + e_kin + e_inter

fig_energy = plt.figure(figsize=(8,4))

plt.title(title + " energy")
plt.plot(sim.time(), e_ph, label="$E_{ph}$")
plt.plot(sim.time(), e_kin, label="$E_{kin}$")
plt.plot(sim.time(), e_inter, label="$E_{inter}$")
plt.plot(sim.time(), e_tot, label="$E_{tot}$")
plt.xlim((0, time_max))
plt.ylabel("energy")
plt.legend()

plt.savefig(output_path + "energy.png", dpi=300)

plt.figure()
plt.title("num configs")
plt.plot(sim.time(), sim.fullStorage_num_configs)


# create the video
ffmpeg_call = "ffmpeg -y -framerate %g -s 1920x1080 "\
              "-i images_for_video/pic%%d.png -crf 20 "\
              "-vf tpad=stop_mode=clone:stop_duration=%g -vcodec libx264 %s"\
              %(sim.time().size/video_duration, last_frame_hold_duration, video_name)
              #
print(ffmpeg_call)
os.system(ffmpeg_call)  #sim.time().size

plt.show()

