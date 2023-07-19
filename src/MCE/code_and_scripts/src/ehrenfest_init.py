import numpy as np
from scipy import stats#, optimize
from scipy.special import erfinv
import multiprocessing as mp
import os # for random number generation

#x_shift must be broadcastable to (num_sites_y, num_sites_x)
def initial_state_MTE(num_sites_x, num_sites_y, num_samples, symmetrize_start, omega, x_shift):
	scale_x = 1./np.sqrt(2)
	scale_p = 1./np.sqrt(2)
	
	if num_samples == 0:
		y = np.zeros((1, num_sites_y, num_sites_x, 2), dtype=np.float)
	else:
		y = np.empty(((4 if symmetrize_start else 1)*num_samples, num_sites_y, num_sites_x, 2), dtype=np.float)
		
		y[:num_samples, :, :, 0] = np.random.normal(loc=0,
		                                            scale=scale_x,
		                                            size=(num_samples, num_sites_y, num_sites_x))
		y[:num_samples, :, :, 1] = np.random.normal(loc=0,
		                                            scale=scale_p,
		                                            size=(num_samples, num_sites_y, num_sites_x))
		
		if symmetrize_start:
			# enforce symmetry by creating 4 counter parts
			y[1*num_samples:2*num_samples, :, :, 0] = -y[:num_samples, :, :, 0]
			y[1*num_samples:2*num_samples, :, :, 1] =  y[:num_samples, :, :, 1]
			
			y[2*num_samples:3*num_samples, :, :, 0] =  y[:num_samples, :, :, 0]
			y[2*num_samples:3*num_samples, :, :, 1] = -y[:num_samples, :, :, 1]
			
			y[3*num_samples:4*num_samples, :, :, 0] = -y[:num_samples, :, :, 0]
			y[3*num_samples:4*num_samples, :, :, 1] = -y[:num_samples, :, :, 1]
		else:
			# eliminate random varaition from mean 0
			y -= y.mean(axis=0).reshape(1, num_sites_y, num_sites_x, 2)
	
	# scale to variance
	if y.shape[0] >= 2:
		y[:, :, :, 0] *= scale_x/np.sqrt((y[:, :, :, 0]**2).sum(axis=0)/y.shape[0]).reshape(1, num_sites_y, num_sites_x)
		y[:, :, :, 1] *= scale_p/np.sqrt((y[:, :, :, 1]**2).sum(axis=0)/y.shape[0]).reshape(1, num_sites_y, num_sites_x)
	# shift to desired mean given by x_shift
	y[:, :, :, 0] += x_shift.reshape(1, *x_shift.shape)
	return y

# helper function to align a numpy array to any alignment (currently unused)
def aligned(a, alignment):
	if (a.ctypes.data % alignment) == 0:
		return a
	assert alignment % a.itemsize == 0
	extra = alignment // a.itemsize
	buf = np.empty(a.size + extra, dtype = a.dtype)
	ofs = (-buf.ctypes.data % alignment) // a.itemsize
	aa = buf[ofs:ofs + a.size].reshape(a.shape)
	np.copyto(aa, a)
	assert aa.ctypes.data % alignment == 0
	return aa


def overlap_matrix(hoz):
	num_sites = hoz.shape[0]
	num_configs = hoz.shape[1]
	hoz_half_norm = np.sum(hoz.real**2 + hoz.imag**2, axis=0)/2
	
	R_overlap_zz = np.zeros((num_configs, num_configs), dtype=np.complex)
	for n in range(num_sites):
		R_overlap_zz += np.outer(np.conjugate(hoz[n]), hoz[n])
	
	R_overlap_zz -= hoz_half_norm.reshape(1, -1)
	R_overlap_zz -= hoz_half_norm.reshape(-1, 1)
	
	R_overlap_zz = np.exp(R_overlap_zz)
	return R_overlap_zz


# eval_overlap_integral only has an effect if psi_init != "local".
# in that case, it determines whether the electroniv overlap integral of the electron eigenstate
# is approximated at the center or (monte carlo) integrated over the full electronic hilbertspace
def _trial_initial_state_MCE(scale_arr, num_sites_x, num_sites_y, num_configs, max_num_configs,
                             hopping, omega, gamma, hoz_initial, hoz,
                             psi_init, psi_state, eval_overlap_integral, returnNormOnly,
                             mc_error_tolerance = 0.0001, flip_initial_state_for_half_configs=False,
                             pancakeSampling=False, pancakeScale=1.,
                             local_psi_init_select_eigenstate=False, printProgress=True,
                             returnCondition=False):
	# TODO allow to disable teh bias in psi?
	scale_ph = scale_arr[0]
	scale_el = scale_arr[1]
	num_sites = num_sites_x*num_sites_y
	
	# setup for overlap integration
	use_own_monte_carlo = num_sites > 3
	overlap_integral_num_processes = 4
	
	mc_max_num_reps = 1000  # max number, might quit early if tolerance is fullfilled
	mc_min_num_reps = 1
	mc_batch_size = 2**14
	
	use_quad_integral = not use_own_monte_carlo
	quad_int_points_per_site = 100
	
	quad_integral_weighted_x = True
	quad_integral_weighted_range = .98
	quad_integral_weighted_scale = 3.
	
	quad_integral_linear_range = 4.
	
	
	show_integral_result = False
	show_integrand = num_sites == 2 and False
	
	psi_eigen_use_optimal = False # slight preference to False
	use_new_integrand_shift = False # prefer False
	
	if not use_quad_integral and not use_own_monte_carlo and psi_init != "local":
		raise ValueError("no integration scheme specified")
	
	
	R_overlap_zz = overlap_matrix(hoz[:, :num_configs])
	
	if psi_init == "local":
		psi_initial = np.zeros((num_sites_y, num_sites_x), dtype=np.complex)
		psi_initial[psi_state] = 1
		psi_initial = psi_initial.flatten()
	
	# initialization of electron coordinates
	if psi_init == "local" and not local_psi_init_select_eigenstate:
		psi_bias = np.zeros((num_sites_y, num_sites_x), dtype=np.complex)
		psi_bias[psi_state] = 1
		
		"""
		X, Y = np.meshgrid(np.arange(num_sites_x), np.arange(num_sites_y))
		x = X[el_grid].flatten()
		y = Y[el_grid].flatten()
		
		#for i in range(len(y)):
		psi_bias[y[0], x[0]] = 1
		"""
		
		psi_bias = psi_bias.reshape(num_sites, 1)
		
	else:
		coupling_strength = gamma*np.sqrt(2)
		psi_bias = np.empty((num_sites, max_num_configs), dtype=np.complex)
		
		H = np.zeros((num_sites, num_sites))
		diag_ind = np.arange(num_sites)
		if num_sites_x == 1 or num_sites_y == 1:
			H[diag_ind[:-1], diag_ind[1:]] = -hopping
			H[diag_ind[1:], diag_ind[:-1]] = -hopping
		else:
			# 2d Hamiltonian, just free hopping part
			# TODO rewrite in the pythonic way?
			for y in range(num_sites_y):
				for x in range(num_sites_x):
					site = y*num_sites_x + x
					if x != num_sites_x-1:
						H[site+1, site] = -hopping
						H[site, site+1] = -hopping
					if y != num_sites_y-1:
						H[site+num_sites_x, site] = -hopping
						H[site, site+num_sites_x] = -hopping
		
		if local_psi_init_select_eigenstate and printProgress:
			popcount = np.zeros(num_sites)
		
		for c in range(num_configs):
			H[diag_ind, diag_ind] = (-2*gamma)*hoz[:, c].real
			
			# find eigenvalues and eigenstates
			eigval, eigvecs = np.linalg.eigh(H)
			order = np.argsort(eigval)
			
			# use either lowest or highest eigenstates as psi_bias
			if psi_init == "eigen" and psi_state < num_sites and 0 <= psi_state:
				psi_bias[:, c] = eigvecs[:,psi_state]
			elif psi_init == "low":
				psi_bias[:, c] = eigvecs[:, 0]
			elif psi_init == "high":
				psi_bias[:, c] = eigvecs[:, -1]
			elif psi_init == "local" and local_psi_init_select_eigenstate:
				eigvec_prob = np.abs(eigvecs.T @ psi_initial)**2
				if scale_ph == 0. and num_configs <= num_sites:
					order = np.argsort(eigvec_prob)[::-1]
					psi_bias[:, :num_configs] = eigvecs[:, order[:num_configs]]
					break
				else:
					# select randomly one of the eigenstates according to projection
					eigvec_prob_cumsum = np.cumsum(eigvec_prob)
					state = np.argmax(eigvec_prob_cumsum > np.random.random())
					#print(eigvec_prob_cumsum, f"state={state}")
					psi_bias[:, c] = eigvecs[:, state]
					if printProgress:
						popcount[state] += 1
					#print(psi_bias[:, c], state)
			else:
				raise ValueError("invalid psi_bias")
		
		if local_psi_init_select_eigenstate and printProgress:
			print("selected eigenstates", popcount)
		
		
		H[diag_ind, diag_ind] = (-2*gamma)*hoz_initial.flatten().real
		psi_initial_mid_point_conj = np.linalg.eigh(H)[1][:, psi_state]
	
	# random electronic basis states
	el_coeff_displacement = np.random.normal(loc=0, scale=1, size=(num_sites, num_configs, 2))\
	                        .view(dtype=np.complex).reshape(num_sites, num_configs)
	
	el_coeff_displacement -= np.mean(el_coeff_displacement, axis=1, keepdims=True)
	if num_configs > 1:
		el_coeff_displacement *= scale_el/np.std(el_coeff_displacement, axis=1, keepdims=True, ddof=1)
	
	el_coeff = np.empty((num_sites, max_num_configs), dtype=np.complex)
	el_coeff[:, :num_configs] = el_coeff_displacement
	if flip_initial_state_for_half_configs and psi_bias.shape[1] == 1 and not local_psi_init_select_eigenstate:
			el_coeff[:, :num_configs//2] += psi_bias[:, :]
			el_coeff[:, num_configs//2:] += psi_bias[::-1, :]
	else:
		el_coeff[:, :num_configs] += psi_bias[:, :num_configs]
	
	# normalize
	el_coeff[:, :num_configs] /= np.linalg.norm(el_coeff[:, :num_configs], axis=0).reshape(1, num_configs)
	
	# total_initial_overlap
	if psi_init == "local" or not eval_overlap_integral:
		# product state -> calculate seperatly
		
		# electronic contribution 
		if psi_init == "local":
			el_initial_overlap = np.sum(np.conjugate(el_coeff[:, :num_configs]) * psi_initial.reshape(-1, 1), axis=0)
		else:
			el_initial_overlap = np.sum(np.conjugate(el_coeff[:, :num_configs]) * psi_bias[:, :num_configs], axis=0)
		
		# phononic contribution
		hoz_initial_overlap = np.zeros(num_configs, dtype=np.complex)

		for l in range(num_sites):
			hoz_initial_overlap += np.conjugate(hoz[l, :num_configs])*hoz_initial[l]

		hoz_initial_overlap -= np.sum(np.abs(hoz[:, :num_configs])**2, axis=0)/2
		hoz_initial_overlap -= np.sum(np.abs(hoz_initial)**2)/2

		hoz_initial_overlap = np.exp(hoz_initial_overlap)
		
		total_initial_overlap = el_initial_overlap * hoz_initial_overlap
	
	else: # psi_init == "eigen"
		if num_sites_y != 1:
			raise ValueError("overlap integral with eigenstates are only supported for 1D systems")
		
		def integrand(x_shifted, hoz, hoz_initial, el_coeff, configs=None):
			if configs is None:
				configs = np.arange(num_configs)
			# x_shifted contains sqrt(2)*x, scaled relative to phonon coordinates
			# x_shifted.shape == (num_evals, num_sites)
			num_evals = x_shifted.shape[0]
			
			if use_new_integrand_shift:
				# new version assumes shift to 2*Re(z_init,n):
				x_restored = x_shifted + 2*hoz_initial.real.reshape(1, num_sites)
				
				H = np.zeros((num_evals, num_sites, num_sites))
				diag_ind = np.arange(num_sites)
				H[:, diag_ind[:-1], diag_ind[1:]] = -hopping
				H[:, diag_ind[1:], diag_ind[:-1]] = -hopping
				H[:, diag_ind, diag_ind] = (-gamma)*x_restored
				
				eigvals, eigvecs = np.linalg.eigh(H)
				
				x_restored = x_restored.reshape(num_evals, 1, num_sites)
				
				state_init = eigvecs[:, :, psi_state].reshape(num_evals, 1, num_sites)
				state_init[np.sum(psi_initial_mid_point_conj.reshape(1, 1, num_sites)*state_init, axis=-1) < 0, :] *= -1
				#state_init = np.ones((1, 1, num_sites))
				
				x_shifted = x_shifted.reshape(num_evals, 1, num_sites)
				
				
				phononic_weights = np.exp(np.sum(-.5*((x_shifted + hoz_initial.real.reshape(1, num_sites))-hoz[:, configs].real.T.reshape(1, num_configs, num_sites))**2\
				                                 +1j*x_restored*(hoz_initial.imag - hoz[:, configs].imag).T.reshape(1, configs.size, num_sites)\
				                                 -.5*(hoz_initial.real - hoz[:, configs].real).T.reshape(1, configs.size, num_sites)**2 \
				                                 +1j*(hoz[:, configs].real*hoz[:, configs].imag).T.reshape(1, configs.size, num_sites)\
				                                 -1j*(hoz_initial.real*hoz_initial.imag).T.reshape(1, 1, num_sites), axis=2))
				
				# return shape = (num_evals, configs.size, num_sites)
				return phononic_weights.reshape(num_evals, configs.size, 1)*state_init/(2*np.pi)**(num_sites/2)
				
			else:
				# OLD version assumes shift to Re(z_jn) + Re(z_init,n):
				
				x_shifted = x_shifted.reshape(num_evals, 1, num_sites)
				
				# shift x back to the correct value (shifted center of gaussian to 0 in the integration scheme)
				x_restored = x_shifted\
					         + (hoz[:, configs].real+hoz_initial.real).T.reshape(1, configs.size, num_sites)
				#x_restored.shape == (num_evals, configs.size, num_sites)
				
				# last term in integral, electron overlap
				# TODO only works in 1D
				# calculate Born-Oppenheimer surface at each x
				H = np.zeros((num_evals, configs.size, num_sites, num_sites))
				diag_ind = np.arange(num_sites)
				H[:, :, diag_ind[:-1], diag_ind[1:]] = -hopping
				H[:, :, diag_ind[1:], diag_ind[:-1]] = -hopping
				
				# hoz[:, configs].shape == (num_sites, configs.size)
				H[:, :, diag_ind, diag_ind] = (-gamma)*x_restored
				eigvals, eigvecs = np.linalg.eigh(H)
				
				if psi_init == "eigen" and psi_state < num_sites and 0 <= psi_state:
					state_init = eigvecs[:, :, :, psi_state]
				elif psi_init == "low":
					state_init = eigvecs[:, :, :, 0]
				elif psi_init == "high":
					state_init = eigvecs[:, :, :, -1]
				else:
					raise Exception("invalid value of psi_init, but should have thrown before")
				
				# fix sign of eigenvector
				for c in range(configs.size):
					state_init[np.sum(psi_initial_mid_point_conj.reshape(1, num_sites)*(state_init[:, c, :]), axis=-1) < 0, c, :] *= -1
				
				# hoz_initial.shape = (num_sites, 1)
				# phononic prefactor
				phononic_weights = np.exp(np.sum(-.5*x_shifted**2\
				                                 +1j*x_restored*(hoz_initial.imag - hoz[:, configs].imag).T.reshape(1, configs.size, num_sites)\
				                                 -.5*(hoz_initial.real - hoz[:, configs].real).T.reshape(1, configs.size, num_sites)**2 \
				                                 +1j*(hoz[:, configs].real*hoz[:, configs].imag).T.reshape(1, configs.size, num_sites)\
				                                 -1j*(hoz_initial.real*hoz_initial.imag).T.reshape(1, 1, num_sites), axis=2))
				#phononic_weights.shape = (num_evals, configs.size)
				
				# return shape = (num_evals, configs.size, num_sites)
				return phononic_weights.reshape(num_evals, configs.size, 1)*state_init/(2*np.pi)**(num_sites/2)
		
		
		if use_own_monte_carlo:
			def gauss(x):
				norm = stats.norm(0, 1)
				return np.prod(norm.pdf(x), axis=1)
			
			optimal_el_coeff = np.zeros((num_configs, num_sites, 2))
			optimal_el_coeff_var = np.full((num_configs, num_sites, 2), np.inf)
			
			all_configs = np.arange(num_configs)
			
			def mc_integrate_batch(queue_data_transfer, mc_batch_size, configs, num_sites, hoz, hoz_initial, el_coeff):
				np.random.seed(int.from_bytes(os.urandom(4), byteorder='little'))
				
				x = np.random.normal(loc=0, scale=1, size=(mc_batch_size, num_sites))
				result = (integrand(x, hoz, hoz_initial, el_coeff, configs)/gauss(x).reshape(-1, 1, 1))\
				         .view(dtype=np.float).reshape(mc_batch_size, configs.size, num_sites, 2)
				
				
				optimal_el_coeff_run_mean = np.mean(result, axis=0)
				optimal_el_coeff_run_var = np.var(result, axis=0, ddof=1)/mc_batch_size
				#print(optimal_el_coeff_run_mean, optimal_el_coeff_run_var)
				if queue_data_transfer is None:
					return optimal_el_coeff_run_mean, optimal_el_coeff_run_var
				else:
					queue_data_transfer.put((optimal_el_coeff_run_mean, optimal_el_coeff_run_var))
			
			print("running Monte Carlo integral")
			q = mp.Queue()
			for i in range(mc_max_num_reps):
				# only sample those configs with high error estimate
				if i < mc_min_num_reps:
					configs = all_configs
				else:
					configs = all_configs[np.any(optimal_el_coeff_var > mc_error_tolerance**2, axis=(1,2))]
				if configs.size == 0:
					print("MC integral %i steps"%i)
					break
				print(i, configs.size)
				
				processes = []
				for i in range(overlap_integral_num_processes):
					p = mp.Process(target=mc_integrate_batch, daemon=True,
					               args=(q, mc_batch_size, configs, num_sites, hoz, hoz_initial, el_coeff))
					p.start()
					processes.append(p)
				
				for i in range(overlap_integral_num_processes):
					optimal_el_coeff_run_mean, optimal_el_coeff_run_var = q.get()
					
					optimal_el_coeff[configs] = (optimal_el_coeff[configs]/optimal_el_coeff_var[configs]\
					                             + optimal_el_coeff_run_mean/optimal_el_coeff_run_var)
					optimal_el_coeff_var[configs] = 1/(1/optimal_el_coeff_var[configs]+1/optimal_el_coeff_run_var)
					optimal_el_coeff[configs] *= optimal_el_coeff_var[configs]
				
				for p in processes:
					p.join()
				
			optimal_el_coeff = optimal_el_coeff.view(dtype=np.complex).reshape(num_configs, num_sites)
			
			mc_optimal_el_coeff = optimal_el_coeff
		
		if use_quad_integral:
			# setup grid for integration
			if quad_integral_weighted_x:
				eval_points = quad_integral_weighted_scale*erfinv(np.linspace(-quad_integral_weighted_range, quad_integral_weighted_range, quad_int_points_per_site))
				dx = np.empty(eval_points.size)
				dx[1:] = (eval_points[1:]-eval_points[:-1])
				dx[0] = dx[-1]
				
				dx[1:-1] = .5*(eval_points[2:]-eval_points[:-2])
				dx[0] = dx[1]
				dx[-1] = dx[-2]
			else:
				quad_int_range = (-quad_integral_linear_range, quad_integral_linear_range)
				eval_points, dx_scalar = np.linspace(*quad_int_range, quad_int_points_per_site, retstep=True)
				dx = np.full(num_sites, dx_scalar)
			
			# recursive integrate over all sites, the last 2 are calculated in one go
			def quad_integrate(dim, prev_x, eval_points, dx, configs, printProgress=True, printNumber=20):
				if printProgress and dim == num_sites:
					print("running quadrature integral")
				
				direct_eval_dim = 2
				if dim <= direct_eval_dim:
					x = np.empty((eval_points.size**direct_eval_dim, num_sites))
					x[:, num_sites-direct_eval_dim:] = np.array(np.meshgrid(*([eval_points]*direct_eval_dim))).T.reshape(-1, direct_eval_dim)
					x[:, :num_sites-direct_eval_dim] = prev_x.reshape(1, -1)
					
					result = integrand(x, hoz, hoz_initial, el_coeff, configs)
					for s in range(direct_eval_dim):
						result = np.sum( dx.reshape(quad_int_points_per_site, 1, 1)
						                *result.reshape(eval_points.size, -1, num_sites), axis=0)
				else:
					result = np.zeros((configs.size, num_sites), dtype=np.complex)
					curr_x = np.empty(prev_x.size+1)
					curr_x[:-1] = prev_x
					
					print_dist = int(np.ceil(eval_points.size/printNumber))
					for i in range(eval_points.size):
						curr_x[-1] = eval_points[i]
						result += dx[i]*quad_integrate(dim-1, curr_x, eval_points, dx, configs)
						if printProgress and dim == num_sites and (i+1)%print_dist == 0:
							print("%i/%i"%(i+1, eval_points.size))
				
				return result
			
			def mp_quad_integrate(queue_data_transfer, num_sites, hoz, hoz_initial, el_coeff,
			                      eval_points, dx, configs, printProgress=True, printNumber=20):
				optimal_el_coeff = quad_integrate(num_sites, np.zeros(0), eval_points, dx, configs, printProgress, printNumber)
				queue_data_transfer.put((configs, optimal_el_coeff))
			
			
			optimal_el_coeff = np.zeros((num_configs, num_sites), dtype=np.complex)
			
			q = mp.Queue()
			processes = []
			for i in range(overlap_integral_num_processes):
				configs = np.arange(i, num_configs, overlap_integral_num_processes)
				p = mp.Process(target=mp_quad_integrate, daemon=True,
				               args=(q, num_sites, hoz, hoz_initial, el_coeff, eval_points, dx, configs, i == 0))
				p.start()
				processes.append(p)
			
			for i in range(overlap_integral_num_processes):
				configs, optimal_el_coeff_p = q.get()
				optimal_el_coeff[configs] = optimal_el_coeff_p
			
			for p in processes:
				p.join()
			
			
			if show_integrand:
				import matplotlib.pyplot as plt
				#config = 0
				# recalc, but needed to not break the program
				result = integrand(x, hoz, hoz_initial, el_coeff)
				for config in range(num_configs):
					for s in range(num_sites):
						fig, axs = plt.subplots(2, 1)
						fig.suptitle("overlap integrand, config %i, site %i"%(config, s))
						real_values = result.reshape(quad_int_points_per_site, quad_int_points_per_site, num_configs, num_sites)[:, :, config, s].real.T
						vmax = max(abs(real_values.min()), abs(real_values.max()))
						cm = axs[0].pcolormesh(eval_points, eval_points, real_values, shading="nearest", vmax=vmax, vmin=-vmax, cmap="RdBu")
						cb = fig.colorbar(cm, ax=axs[0])
						cb.set_label("Real")
						
						imag_values = result.reshape(quad_int_points_per_site, quad_int_points_per_site, num_configs, num_sites)[:, :, config, s].imag.T
						vmax = max(abs(imag_values.min()), abs(imag_values.max()))
						cm = axs[1].pcolormesh(eval_points, eval_points, imag_values, shading="nearest", vmax=vmax, vmin=-vmax, cmap="RdBu")
						cb = fig.colorbar(cm, ax=axs[1])
						cb.set_label("Imag")
						
						# expected center points
						axs[0].plot(0, 0, 'k+')
						axs[1].plot(0, 0, 'k+')
						
						axs[0].plot(hoz[0, config].real - hoz_initial.real.ravel()[0], hoz[1, config].real - hoz_initial.real.ravel()[1], 'kx')
						axs[1].plot(hoz[0, config].real - hoz_initial.real.ravel()[0], hoz[1, config].real - hoz_initial.real.ravel()[1], 'kx')
					
					plt.show()
			
			
			quad_optimal_el_coeff = optimal_el_coeff
		
		if psi_eigen_use_optimal:
			el_coeff[:, :num_configs] = optimal_el_coeff.T
			el_coeff[:, :num_configs] /= np.linalg.norm(el_coeff[:, :num_configs], axis=0).reshape(1, num_configs)
			
			el_coeff[:, :num_configs] += el_coeff_displacement
			el_coeff[:, :num_configs] /= np.linalg.norm(el_coeff[:, :num_configs], axis=0).reshape(1, num_configs)
			
		total_initial_overlap = np.sum(np.conjugate(el_coeff[:, :num_configs].T).reshape(num_configs, num_sites) * optimal_el_coeff, axis=1)
		
		if show_integral_result:
			import matplotlib.pyplot as plt
			# shapes of *_optimal_el_coeff = (num_configs, num_sites)
			ms = 9
			fig, axs = plt.subplots(1, num_sites)
			#plt.plot([0, num_configs-1], [0, 0], 'k')
			#plt.plot(results.real.T, '.', alpha=0.4)
			fig.suptitle("real-config")
			if use_own_monte_carlo:
				for s in range(num_sites):
					axs[s].errorbar(np.arange(num_configs), mc_optimal_el_coeff[:, s].real,\
					             yerr=np.sqrt(optimal_el_coeff_var[:, s, 0]),
					             marker='x', ms=ms, ls='none')
			if use_quad_integral:
				for s in range(num_sites):
					axs[s].plot(quad_optimal_el_coeff.real[:, s], '+', markersize=ms)
			
			
			fig, axs = plt.subplots(1, num_sites)
			fig.suptitle("imag-config")
			#plt.plot([0, num_configs-1], [0, 0], 'k')
			#plt.plot(results.imag.T, '.', alpha=0.4)
			if use_own_monte_carlo:
				for s in range(num_sites):
					axs[s].errorbar(np.arange(num_configs), mc_optimal_el_coeff[:, s].imag,\
					             yerr=np.sqrt(optimal_el_coeff_var[:, s, 1]),\
					             marker='x', ms=ms, ls='none')
			if use_quad_integral:
				for s in range(num_sites):
					axs[s].plot(quad_optimal_el_coeff.imag[:, s], '+', markersize=ms)
			
			fig, axs = plt.subplots(1, num_sites)
			fig.suptitle("imag-real")
			if use_own_monte_carlo:
				for s in range(num_sites):
					axs[s].errorbar(mc_optimal_el_coeff[:, s].real, mc_optimal_el_coeff[:, s].imag,\
					             xerr=np.sqrt(optimal_el_coeff_var[:, s, 0]),\
					             yerr=np.sqrt(optimal_el_coeff_var[:, s, 1]),\
					             marker='x', ms=ms, ls='none')
			if use_quad_integral:
				for s in range(num_sites):
					axs[s].plot(quad_optimal_el_coeff.real[:, s], quad_optimal_el_coeff.imag[:, s], '+', markersize=ms)
			
			if use_own_monte_carlo and use_quad_integral:
				fig, axs = plt.subplots(1, num_sites)
				fig.suptitle("diff real-config")
				for s in range(num_sites):
					axs[s].plot([0, num_configs], [0,0], 'k')
					axs[s].errorbar(np.arange(num_configs), mc_optimal_el_coeff[:, s].real-quad_optimal_el_coeff.real[:, s],\
					                yerr=np.sqrt(optimal_el_coeff_var[:, s, 0]),\
					                marker='x', ms=ms, ls='none')
				
				fig, axs = plt.subplots(1, num_sites)
				fig.suptitle("diff imag-config")
				for s in range(num_sites):
					axs[s].plot([0, num_configs], [0,0], 'k')
					axs[s].errorbar(np.arange(num_configs), mc_optimal_el_coeff[:, s].imag-quad_optimal_el_coeff.imag[:, s],\
					                yerr=np.sqrt(optimal_el_coeff_var[:, s, 1]),\
					                marker='x', ms=ms, ls='none')
				
			plt.show()
		
		
		#print("total_initial_overlap", total_initial_overlap)
		#print("abs(total_initial_overlap)", np.abs(total_initial_overlap))
		#print("arg(total_initial_overlap)", np.angle(total_initial_overlap))
	
	
	el_overlap_matrix = np.zeros((num_configs, num_configs), dtype=np.complex)
	for n in range(num_sites):
		el_overlap_matrix += np.outer(np.conjugate(el_coeff[n, :num_configs]), el_coeff[n, :num_configs])
	
	
	total_overlap_matrix = el_overlap_matrix * R_overlap_zz# yes, elementwise multiply
	
	coupling_coeff = np.linalg.solve(total_overlap_matrix, total_initial_overlap)
	
	el_coeff_scaled = el_coeff[:, :num_configs] * coupling_coeff.reshape(1, num_configs)
	
	# calculate total norm
	total_norm = np.sum(np.conjugate(el_coeff_scaled).reshape(num_sites, 1, num_configs)
	                    @ (R_overlap_zz.reshape(1, num_configs, num_configs)
	                       @ el_coeff_scaled.reshape(num_sites, num_configs, 1))).real
	
	
	if returnNormOnly:
		return total_norm
	else:
		#print("min weight", np.abs(el_coeff[:, :num_configs]).min())
		#R_overlap_zz = np.abs(R_overlap_zz)
		#print("min distance:", np.sort(np.sqrt(-2*np.log(R_overlap_zz)).flatten())[num_configs:num_configs+20:2])
		
		hoz[:, num_configs:] = 0
		el_coeff[:, num_configs:] = 0
		
		retVal = total_norm, hoz.reshape(   num_sites_y, num_sites_x, max_num_configs), \
		         el_coeff.reshape(1, num_sites_y, num_sites_x, max_num_configs), coupling_coeff
		
		if returnCondition:
			condition_ph = np.linalg.cond(R_overlap_zz)
			condition_tot = np.linalg.cond(total_overlap_matrix)
			
			retVal = (*retVal, condition_tot, condition_ph)
		
		return retVal

#@profile
def initial_state_MCE(version, num_sites_x, num_sites_y, num_configs, max_num_configs,
                      hopping, omega, gamma, psi_init, psi_state,
                      phonon_dressed_displaced_sites,
                      max_norm_diff=2e-5,
                      min_norm_diff=5e-6,
                      local_psi_init_select_eigenstate=False, # True recommended for MCEv2
                      printProgress=True,
                      returnFinalParams=False,
                      seed=None, allowCompression=True, allowExpansion=None,
                      scale_ph=None, scale_el=None,
                      flip_initial_state_for_half_configs=False,
                      scale_fact_ph=0.9, scale_fact_el=0.9,
                      pancakeSampling=False, pancakeScale=1.,
                      calcCondition=False):# ONLY used to optimize the tunneling dimer
	
	if pancakeSampling and psi_init != "local":
		print("pancake sampling is only supported for local init")
		pancakeSampling = False
	
	if allowExpansion is None:
		allowExpansion = allowCompression
	
	plot_sampled = False
	if version != 1 and version != 2:
		raise ValueError("Version of MCE simulation must be either 1 or 2")
	
	num_sites = num_sites_x*num_sites_y
	repulsive_configs = num_sites <= 3
	#seed = 1792412846
	if seed is not None:
		np.random.seed(seed)
	
	if local_psi_init_select_eigenstate and num_configs <= num_sites:
		if scale_el is not None and scale_el != 0. or scale_ph is not None and scale_ph != 0.:
			print("WARNING: scale_el and scale_ph must be 0 (or None) if ocal_psi_init_select_eigenstate and num_configs <= num_sites."\
			      " Ignoring given value and using 0 instead.")
		scale_ph = 0.
		scale_el = 0.
		
	if scale_ph is None:
		scale_ph = 1.
		# TODO better auto select of initial compression and/or add dynamic expansion
		#if scale_ph < 2. and (local_psi_init_select_eigenstate and (scale_el is None or scale_el < 0.1) and num_configs >= 30*num_sites_x*num_sites_y
		#                      or num_configs >= 60*num_sites_x*num_sites_y):
		#	scale_ph = 2.
	
	if scale_el is None:
		if psi_init == "eigen":
			scale_el = .05 # if converged, this should increase, but small for small num_configs/num_sites
		elif local_psi_init_select_eigenstate:
			scale_el = 0.
		#elif flip_initial_state_for_half_configs:
		#	scale_el = 0.2
		else:
			scale_el = scale_ph/2 #num_sites
	
	
	
	total_norm = 0.
	
	y = np.zeros((num_sites, max_num_configs, 2))
	y[:, :num_configs] = np.random.normal(loc=0, scale=1, size=(num_sites, num_configs, 2))
	# TODO try symmetry
	
	#y[:, :num_configs] -= np.mean(y[:, :num_configs], axis=1, keepdims=True)
	#if num_configs > 1:
	#	y[:, :num_configs] /= np.std(y[:, :num_configs], axis=1, ddof=0, keepdims=True)
	
	#fig, axs = plt.subplots(1, num_sites, squeeze=False)
	#for i in range(num_sites):
	#	axs[0, i].plot(y[i, :num_configs, 0], y[i, :num_configs, 1], '+')
	
	used_params = {"max_norm_diff":max_norm_diff,
	               "local_psi_init_select_eigenstate":local_psi_init_select_eigenstate,
	               "flip_initial_state_for_half_configs":flip_initial_state_for_half_configs,
	               "repulsive_configs":repulsive_configs}
	
	if repulsive_configs and num_configs > 1:
		repulsion_scale = 2
		repulsion_factor = 0.05
		repulsion_num_iter = 2
		
		used_params["repulsion"] = {"num_iter":repulsion_num_iter,
		                            "factor":repulsion_factor,
		                            "scale":repulsion_scale}
		
		for i in range(repulsion_num_iter):
			distances = y[:, :num_configs].transpose(1, 0, 2).reshape(num_configs, 2*num_sites)
			distances = distances.reshape(num_configs, 1, 2*num_sites) - distances.reshape(1, num_configs, 2*num_sites)
			dist_norm = np.linalg.norm(distances, axis=2)**2
			dist_norm[np.diag_indices(num_configs)] = np.inf
			dist_norm = dist_norm.reshape(num_configs, num_configs, 1)
			min_dist = dist_norm.min()
			if printProgress:
				print("min dist", min_dist)
			
			force = distances/dist_norm*np.exp(-repulsion_scale*np.abs(dist_norm))
			force_on_config = np.sum(force, axis=1).reshape(num_configs, num_sites, 2).transpose(1, 0, 2)
			
			y[:, :num_configs] += repulsion_factor*force_on_config
			y[:, :num_configs] /= np.std(y[:, :num_configs], axis=1, ddof=0, keepdims=True)
		
		used_params["repulsion"]["min_dist"] = min_dist
		
	#for i in range(num_sites):
	#	axs[0, i].plot(y[i, :num_configs, 0], y[i, :num_configs, 1], 'x')
	#plt.show()
	
	y[:, :num_configs] -= np.mean(y[:, :num_configs], axis=1, keepdims=True)
	if num_configs > 1:
		y[:, :num_configs] /= np.std(y[:, :num_configs], axis=1, ddof=0, keepdims=True)
	
	y[:, :num_configs] /= 2 # bring to correct scale in hoz language
	
	hoz_raw = y.view(dtype=np.complex).reshape(num_sites, max_num_configs)
	
	
	# initialization of phonon coordinates
	y_initial = np.zeros((num_sites_y, num_sites_x, 2), dtype=np.float)
	
	x_displacement_dressed = gamma/omega*np.sqrt(2)
	
	phonon_dressed_displaced_sites_set = set()
	if phonon_dressed_displaced_sites is not None:
		for site in phonon_dressed_displaced_sites:
			phonon_dressed_displaced_sites_set.add(tuple(site))
	
	for site in phonon_dressed_displaced_sites_set:
		if len(site) == 2 and 0 <= site[0] < num_sites_y and 0 <= site[1] < num_sites_x:
			y_initial[site[0], site[1], 0] = x_displacement_dressed
	
	#print("y_initial", y_initial)
	hoz_initial = (np.array(y_initial)/np.sqrt(2)).view(dtype=np.complex).reshape((num_sites, 1))
	
	# TODO optimize compression finding routine, focus on "eigen"
	# also allow spreading? copression was rescaled by 1/2 to correct wrong wigner scale
	if printProgress:
		print("total_norm\tscale_ph\tscale_el")
	
	allowCompression &= psi_init != "eigen" and num_configs != 1\
	                    and not(local_psi_init_select_eigenstate and num_configs <= num_sites)
	
	seen_too_small_norm = False
	
	if pancakeSampling:
		X, Y = np.meshgrid(np.arange(num_sites_x), np.arange(num_sites_y))
		dist_to_init = np.sqrt((X-psi_state[1])**2 + (Y-psi_state[0])**2)
		pancake_weight = np.exp(-np.abs(dist_to_init)/pancakeScale).reshape(-1, 1)
		#pancake_weight = np.where(np.abs(dist_to_init) <= pancakeScale, 1, 0).reshape(-1, 1)
		
		#print(dist_to_init.shape, dist_to_init)
	else:
		pancake_weight = 1.
	
	while 1-total_norm > max_norm_diff or total_norm == 0.:
		scale_arr = [scale_ph, scale_el*pancake_weight]
		
		scale_ph_sites = scale_ph*pancake_weight
		
		if flip_initial_state_for_half_configs:
			hoz = np.empty(hoz_raw.shape, dtype=np.complex)
			hoz[:, :num_configs//2] = hoz_initial + scale_ph_sites*hoz_raw[:, :num_configs//2]
			hoz[:, num_configs//2:] = hoz_initial[::-1, :] + scale_ph_sites*hoz_raw[:, num_configs//2:]
			#hoz[:, :num_configs//2] -= hoz_initial
			#hoz[:, num_configs//2:] += hoz_initial[::-1, :]
		else:
			hoz = hoz_initial + scale_ph_sites*hoz_raw
		
		
		retVal =\
		_trial_initial_state_MCE(scale_arr, num_sites_x, num_sites_y, num_configs, max_num_configs,
		                         hopping, omega, gamma, hoz_initial, hoz,
		                         psi_init, psi_state, eval_overlap_integral=True, returnNormOnly=False,
		                         flip_initial_state_for_half_configs=flip_initial_state_for_half_configs,
		                         local_psi_init_select_eigenstate=local_psi_init_select_eigenstate,
		                         pancakeSampling=pancakeSampling, pancakeScale=pancakeScale,
		                         printProgress=printProgress, returnCondition=calcCondition)
		
		if calcCondition:
			total_norm, hoz_sim, el_coeff_sim, coupling_coeff, condition_tot, condition_ph = retVal
		else:
			total_norm, hoz_sim, el_coeff_sim, coupling_coeff = retVal
		
		if seed is None:
			used_params["seed"] = None
		else:
			used_params["seed"] = int(seed)
		used_params["scale_ph"] = scale_ph
		used_params["scale_el"] = scale_el
		used_params["pancakeSampling"] = pancakeSampling
		used_params["pancakeScale"] = pancakeScale
		
		used_params["total_norm"] = total_norm
		if calcCondition:
			used_params["condition_tot"] = condition_tot
			used_params["condition_ph"]  = condition_ph
		
		if printProgress:
			print("%10f\t%10f\t%10f"%(total_norm, scale_ph, scale_el), flush=True)
		
		if not allowCompression:
			break
		
		if 1-total_norm > min_norm_diff:
			seen_too_small_norm = True
		if not seen_too_small_norm and allowExpansion:
			total_norm = 0
			scale_ph /= scale_fact_ph
			scale_el /= scale_fact_el
		else:
			scale_ph *= scale_fact_ph
			scale_el *= scale_fact_el
	
	#print(el_coeff_sim.shape)
	#print(np.sort(np.sum(np.abs(el_coeff_sim[0, :, :, :min(num_configs, 20)])**2, axis=(0,1)))*num_configs)
	
	if plot_sampled:
		import matplotlib.pyplot as plt
		fig, axs = plt.subplots(num_sites_y, num_sites_x, squeeze=False)
		fig.suptitle("hoz sampled")
		fig.tight_layout()
		for j in range(num_sites_y):
			for i in range(num_sites_x):
				cm = axs[j, i].scatter(np.sqrt(2)*hoz_sim[j, i, :num_configs].real,
				                       np.sqrt(2)*hoz_sim[j, i, :num_configs].imag,
				                       c=np.abs(el_coeff_sim[0, j, i, :num_configs]), vmin=0)
				axs[j, i].set_xlabel("x")
				axs[j, i].set_ylabel("p")
				#fig.colorbar(cm, ax=axs[i])
		
		fig, axs = plt.subplots(num_sites_y, num_sites_x, squeeze=False)
		fig.suptitle("hoz sampled")
		fig.tight_layout()
		for j in range(num_sites_y):
			for i in range(num_sites_x):
				cm = axs[j, i].plot(np.sqrt(2)*hoz_sim[j, i, :num_configs].real.reshape(1+flip_initial_state_for_half_configs, -1).T,
					                np.sqrt(2)*hoz_sim[j, i, :num_configs].imag.reshape(1+flip_initial_state_for_half_configs, -1).T, '.')
				axs[j, i].set_xlabel("x")
				axs[j, i].set_ylabel("p")
				#fig.colorbar(cm, ax=axs[i])
		
		fig, axs = plt.subplots(num_sites_y, num_sites_x, squeeze=False)
		fig.suptitle("el_coeff sampled")
		fig.tight_layout()
		for j in range(num_sites_y):
			for i in range(num_sites_x):
				cm = axs[j, i].plot(el_coeff_sim[0, j, i, :num_configs].real.reshape(1+flip_initial_state_for_half_configs, -1).T,
					                el_coeff_sim[0, j, i, :num_configs].imag.reshape(1+flip_initial_state_for_half_configs, -1).T, '.')
				axs[j, i].set_xlabel("x")
				axs[j, i].set_ylabel("p")
		
		plt.figure()
		plt.title("coupling coeff")
		plt.plot(coupling_coeff.real.reshape(1, -1), coupling_coeff.imag.reshape(1, -1), '.')
		
		if version == 2:
			plt.show()
	
	#print("hoz_sim=\n", hoz_sim[:, :, :num_configs])
	#print("el_coeff_sim=\n", el_coeff_sim[0, :, :, :num_configs])
	#print("coupling_coeff=", coupling_coeff[:num_configs])
	
	if version == 1:
		el_coeff_sim[:, :, :, :num_configs] *= coupling_coeff.reshape(1, 1, 1, num_configs)
		
		if plot_sampled:
			fig, axs = plt.subplots(num_sites_y, num_sites_x, squeeze=False)
			fig.suptitle("el_coeff scaled")
			fig.tight_layout()
			for j in range(num_sites_y):
				for i in range(num_sites_x):
					cm = axs[j, i].plot(el_coeff_sim[0, j, i, :num_configs].real.reshape(1+flip_initial_state_for_half_configs, -1).T,
					                    el_coeff_sim[0, j, i, :num_configs].imag.reshape(1+flip_initial_state_for_half_configs, -1).T, '.')
			plt.show()
		
		retVal = hoz_sim, el_coeff_sim, None
	elif version == 2:
		R_overlap_zz = overlap_matrix(hoz[:, :num_configs])
		el_coeff = np.empty((max_num_configs, 1, num_sites_y, num_sites_x), dtype=np.complex)
		el_coeff[:num_configs, :, :, :] = el_coeff_sim[:, :, :, :num_configs].transpose(3, 0, 1, 2)
		
		hoz = np.empty((max_num_configs, num_sites_y, num_sites_x), dtype=np.complex)
		hoz[:num_configs, :, :] = hoz_sim[:, :, :num_configs].transpose(2, 0, 1)
		
		#coupling_coeff[:num_configs] = 1/np.sqrt(num_configs)
		
		coupling_coeff_large_array = np.empty(max_num_configs, dtype=np.complex)
		coupling_coeff_large_array[:num_configs] = coupling_coeff
		
		retVal = hoz, el_coeff, coupling_coeff_large_array
	if returnFinalParams:
		retVal = (*retVal, used_params)
	return retVal
