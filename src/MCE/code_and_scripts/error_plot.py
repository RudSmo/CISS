import numpy as np
import matplotlib.pyplot as plt

initial_state = "DCDW"
omega, gamma = .1, 0.4
plot_mean_std = True	# requires constant number of repetitions on full rectangular grid of num_sites, num_samples

data = np.loadtxt("error_ehrenfest/%s_omega=%g_gamma=%g.dat"%(initial_state, omega, gamma))
"""
data stored as [:, i] with i beeing
0 sites
1 num_samples
2 dt
3 error_nel
4 error_nph
5 error_x
6 error_p
7 error_ekin (is nan iff num_sites==2)
"""
set_num_sites = np.unique(data[:, 0])
set_num_samples = np.unique(data[:, 1])
set_time_steps = np.unique(data[:, 2])

title = initial_state + " " + "local" + ": $\\omega=%g$, $\\gamma=%g$"%(omega, gamma)

for dt in set_time_steps:
	# creating and labeling figures and axis:
	fig_const_L, axs_const_L = plt.subplots(2, 2)
	fig_const_L.suptitle(title + ": dt=%g, constant site number"%dt)
	
	axs_const_L[0, 0].set_xlabel("number of samples")
	axs_const_L[0, 0].set_ylabel(r"error $\langle n_{el}\rangle$")
	
	axs_const_L[0, 1].set_xlabel("number of samples")
	axs_const_L[0, 1].set_ylabel(r"error $\langle n_{ph}\rangle$")
	
	axs_const_L[1, 0].set_xlabel("number of samples")
	axs_const_L[1, 0].set_ylabel(r"error $\langle x_l\rangle$")
	
	axs_const_L[1, 1].set_xlabel("number of samples")
	axs_const_L[1, 1].set_ylabel(r"error $\langle p_l\rangle$")
	
	
	
	
	fig_const_S, axs_const_S = plt.subplots(2, 2)
	fig_const_S.suptitle(title + ": dt=%g, constant sample number"%dt)
	
	axs_const_S[0, 0].set_xlabel("number of sites")
	axs_const_S[0, 0].set_ylabel(r"error $\langle n_{el}\rangle$")
	
	axs_const_S[0, 1].set_xlabel("number of sites")
	axs_const_S[0, 1].set_ylabel(r"error $\langle n_{ph}\rangle$")
	
	axs_const_S[1, 0].set_xlabel("number of sites")
	axs_const_S[1, 0].set_ylabel(r"error $\langle x_l\rangle$")
	
	axs_const_S[1, 1].set_xlabel("number of sites")
	axs_const_S[1, 1].set_ylabel(r"error $\langle p_l\rangle$")
	
	
	# actual plotting:
	
	if plot_mean_std:
		data_dt = data[data[:, 2] == dt]
		
		all_inds = np.arange(data_dt.shape[0])
		
		if data_dt.shape[0]%(set_num_sites.size*set_num_samples.size) != 0:
			raise ValueError("dimension mismatch of data_dt for dt=%g"%dt)
		
		indices = np.empty((set_num_sites.size, set_num_samples.size, data_dt.shape[0]//set_num_sites.size//set_num_samples.size), dtype=np.int)
		
		for i in range(set_num_sites.size):
			for j in range(set_num_samples.size):
				indices[i, j] = all_inds[np.logical_and(data_dt[:, 0] == set_num_sites[i], data_dt[:, 1] == set_num_samples[j])]
		
		
		for i in range(set_num_sites.size):
			which = indices[i, :]
					
			axs_const_L[0, 0].errorbar(data_dt[which, 1][:, 0], data_dt[which, 3].mean(axis=1), data_dt[which, 3].std(ddof=1, axis=1), fmt='.', label="%i sites"%set_num_sites[i])
			axs_const_L[0, 1].errorbar(data_dt[which, 1][:, 0], data_dt[which, 4].mean(axis=1), data_dt[which, 4].std(ddof=1, axis=1), fmt='.', label="%i sites"%set_num_sites[i])
			axs_const_L[1, 0].errorbar(data_dt[which, 1][:, 0], data_dt[which, 5].mean(axis=1), data_dt[which, 5].std(ddof=1, axis=1), fmt='.', label="%i sites"%set_num_sites[i])
			axs_const_L[1, 1].errorbar(data_dt[which, 1][:, 0], data_dt[which, 6].mean(axis=1), data_dt[which, 6].std(ddof=1, axis=1), fmt='.', label="%i sites"%set_num_sites[i])
		
		
		for j in range(set_num_samples.size):
			which = indices[:, j]
					
			axs_const_S[0, 0].errorbar(data_dt[which, 0][:, 0], data_dt[which, 3].mean(axis=1), data_dt[which, 3].std(ddof=1, axis=1), fmt='.', label="%i samples"%set_num_samples[j])
			axs_const_S[0, 1].errorbar(data_dt[which, 0][:, 0], data_dt[which, 4].mean(axis=1), data_dt[which, 4].std(ddof=1, axis=1), fmt='.', label="%i samples"%set_num_samples[j])
			axs_const_S[1, 0].errorbar(data_dt[which, 0][:, 0], data_dt[which, 5].mean(axis=1), data_dt[which, 5].std(ddof=1, axis=1), fmt='.', label="%i samples"%set_num_samples[j])
			axs_const_S[1, 1].errorbar(data_dt[which, 0][:, 0], data_dt[which, 6].mean(axis=1), data_dt[which, 6].std(ddof=1, axis=1), fmt='.', label="%i samples"%set_num_samples[j])
	
	else:
		for num_sites in set_num_sites:
			which = np.logical_and(data[:, 0] == num_sites, data[:, 2] == dt)
			
			axs_const_L[0, 0].plot(data[which, 1], data[which, 3], '.', label="%i sites"%num_sites)
			axs_const_L[0, 1].plot(data[which, 1], data[which, 4], '.', label="%i sites"%num_sites)
			axs_const_L[1, 0].plot(data[which, 1], data[which, 5], '.', label="%i sites"%num_sites)
			axs_const_L[1, 1].plot(data[which, 1], data[which, 6], '.', label="%i sites"%num_sites)

		for num_samples in set_num_samples:
			which = np.logical_and(data[:, 1] == num_samples, data[:, 2] == dt)
			
			axs_const_S[0, 0].plot(data[which, 0], data[which, 3], '.', label="%i samples"%num_samples)
			axs_const_S[0, 1].plot(data[which, 0], data[which, 4], '.', label="%i samples"%num_samples)
			axs_const_S[1, 0].plot(data[which, 0], data[which, 5], '.', label="%i samples"%num_samples)
			axs_const_S[1, 1].plot(data[which, 0], data[which, 6], '.', label="%i samples"%num_samples)
	
	
	for ax in axs_const_L.flatten():
		ax.legend()
		ax.set_ylim(ymin=0)
	fig_const_L.tight_layout()
	
	for ax in axs_const_S.flatten():
		ax.legend()
		ax.set_ylim(ymin=0)
	fig_const_S.tight_layout()

plt.show()
