import numpy as np
import matplotlib.pyplot as plt

import data_manager

# latex style
use_tex = False
if use_tex:
	from matplotlib import rc
	rc('text', usetex=use_tex)
	rc('font', family='serif')

initial_state = "bare"
psi_initial="CDW"
hopping = 10
omega = 1
gamma = 4
num_sites_x = 13
num_sites_y = 1
num_electrons = num_sites_x//2
psi_state = None
 
# number of ehrenfest trajectories (randomly chosen)
# use num_samples = 0 for single trajectory starting in rest on all sites
num_samples = 800
symmetrize_start = False

show_compare_single_traj = False	# only availbale for 2 sites gamma=0.4 and dt=0.05
show_inv_gamma_times = False

show_energy_violation = False

# dt is the save distance, simulation distance may be smaller
dt = 0.01
t_max = 6.5/hopping #20*np.pi
# error in energy conservation \approx 0.2*omega*gamma^2*dt^2 for BCDW
# DCDW comparable or lower

title = initial_state + " " + psi_initial + ": $t_0=%g, \\omega_0=%g$, $\\gamma=%g$"%(hopping, omega, gamma)

try:
	ref = data_manager.RefData(initial_state, psi_initial, psi_state, num_sites_x, num_sites_y, hopping, omega, gamma)
	ref.set_cut_time(t_max)
except Exception as e:
	print("no reference Data found:")
	print(repr(e))
	ref = None

if ref is not None:
	t_max = ref.time()[-1]

"""
hopping, omega, gamma, time_max,
             num_samples=-1, num_sites_x=0, num_sites_y=1,
             psi_initial, initial_state, psi_state=None, num_electrons=-1, 
             time_step=None, tol_energy_error_per_site=1e-6,
             symmetric_y_init=True, y_initial=None,
             fullStorage=False, calcYSquared=True, hopping=1., printProgress=True):
"""

sim = data_manager.Simulation_MTE(hopping, omega, gamma,
                                  num_sites_x, num_sites_y,
                                  t_max, dt,
                                  num_samples,
                                  initial_state, psi_initial, psi_state,
                                  phonon_dressed_displaced_sites=None,
                                  calcYSquared=True,
                                  fullStorage=False,
                                  symmetric_y_init=symmetrize_start,
                                  force_new_simulation=False,
                                  load_sampled_state=False,
                                  store_results=True)


num_sites = num_sites_x*num_sites_y
if num_sites > 2:
	t_edges = np.linspace(sim.time()[0]-dt/2, sim.time()[-1]+dt/2, sim.time().size+1)
	site_edges = np.linspace(-.5, num_sites-.5, num_sites+1)
	if ref is not None:
		t_edges_ref = np.linspace(ref.time()[0]-dt/2, ref.time()[-1]+dt/2, ref.time().size+1)

# plot electronic density
plt.figure(figsize=(5, 5))
plt.title(title + ", electronic density")
for i in range(num_sites):
	plt.plot(sim.time(), sim.nel()[:, i], label="site %i"%i)

if ref is not None:
	plt.gca().set_prop_cycle(None)
	for i in range(num_sites):
		plt.plot(ref.time(), ref.nel(i), "--", label="ref site %i"%i)

if show_inv_gamma_times and gamma != 0:
	for i in range(1, int(np.round(t_max*gamma))):
		plt.plot([i/gamma, i/gamma], [0, 1], 'k')

plt.ylabel(r"$\langle n_{el}\rangle$")
plt.xlim((0,t_max))
plt.ylim((0,1))
if num_sites <= 2:
	plt.legend()

if num_sites > 2:
	fig_nel, axs_nel = plt.subplots(2, 1, squeeze=False, figsize=(5, 8))
	axs_nel = axs_nel.flatten()
	fig_nel.suptitle(title + ", electronic density")
	axs_nel[0].set_title("MTE Stefan")
	axs_nel[1].set_title("ref")
	cm = axs_nel[0].pcolormesh(t_edges, site_edges, sim.nel().T, vmin=0., vmax=1.)
	axs_nel[0].set_ylabel("site")
	if ref is not None:
		axs_nel[1].pcolormesh(t_edges_ref, site_edges, ref.nel().T, vmin=0., vmax=1.)
	axs_nel[1].set_xlabel("time")
	axs_nel[1].set_ylabel("site")
	cb = fig_nel.colorbar(cm, ax=axs_nel[:], shrink=0.8)
	cb.set_label("electronic density")


# plot harmonic oscillator movement
one_row = num_sites <= 2 or ref is None
fig_osc, axs_osc = plt.subplots(1 if one_row else 2, 2,
                                squeeze=False, figsize=(10, 5 if one_row else 8))
fig_osc.suptitle(title + ", oscillator movements")

axs_osc[0, 0].set_xlabel("time")
axs_osc[0, 0].set_ylabel("position x/l0")
axs_osc[0, 1].set_xlabel("time")
axs_osc[0, 1].set_ylabel("momentum p*l0")

if num_sites <= 2 :
	for i in range(num_sites):
		axs_osc[0, 0].plot(sim.time(), sim.x()[:, i], label="site %i"%i)
		axs_osc[0, 1].plot(sim.time(), sim.p()[:, i])
		
	if ref is not None:
		axs_osc[0, 0].set_prop_cycle(None)
		axs_osc[0, 1].set_prop_cycle(None)
		for i in range(num_sites):
			axs_osc[0, 0].plot(ref.time(), ref.x(i), '--')
			axs_osc[0, 1].plot(ref.time(), ref.p(i), '--')
	
	if show_inv_gamma_times and gamma != 0:
		for i in range(1, int(np.round(t_max*gamma))):
			axs_osc[0, 0].plot([i/gamma, i/gamma], [0, 1], 'k')
			axs_osc[0, 1].plot([i/gamma, i/gamma], [0, 1], 'k')
	
	fig_osc.legend()
	fig_osc.tight_layout(rect=(0,0,1,.95))
else:
	# displacement
	c_min = sim.x().min()
	c_max = sim.x().max()
	
	if ref is not None:
		c_min = min(c_min, ref.x().min())
		c_max = min(c_max, ref.x().max())
	
	cm = axs_osc[0, 0].pcolormesh(t_edges, site_edges, sim.x().T, vmin=c_min, vmax=c_max)
	if ref is not None:
		axs_osc[1, 0].pcolormesh(t_edges_ref, site_edges, ref.x().T, vmin=c_min, vmax=c_max)
		axs_osc[1, 0].set_xlabel("time")
		axs_osc[1, 0].set_ylabel("ref position x/l0")
		
		fig_osc.colorbar(cm, ax=axs_osc[:, 0], shrink=0.8)
	else:
		fig_osc.colorbar(cm, ax=axs_osc[:, 0])
	
	# momentum
	c_min = sim.p().min()
	c_max = sim.p().max()
	
	if ref is not None:
		c_min = min(c_min, ref.p().min())
		c_max = min(c_max, ref.p().max())
	
	cm = axs_osc[0, 1].pcolormesh(t_edges, site_edges, sim.p().T, vmin=c_min, vmax=c_max)
	if ref is not None:
		axs_osc[1, 1].pcolormesh(t_edges_ref, site_edges, ref.p().T, vmin=c_min, vmax=c_max)
		axs_osc[1, 1].set_xlabel("time")
		axs_osc[1, 1].set_ylabel("ref momentum p*l0")
		
		fig_osc.colorbar(cm, ax=axs_osc[:, 1], shrink=0.8)
	else:
		fig_osc.colorbar(cm)


# p-x phasespace
plt.figure()
plt.title(title + ", p-x phasespace")
for i in range(num_sites):
	plt.plot(sim.x()[:, i], sim.p()[:, i], label="site %i"%i)

if ref is not None:
	plt.gca().set_prop_cycle(None)
	plt.plot(ref.x(), ref.p(), '--')

plt.legend()
plt.xlabel("x/l0")
plt.ylabel("p*l0")


run_label = "ehrenfest %s"%("single" if num_samples==0 else
	                        "%s%i samples"%("4*" if symmetrize_start else "", num_samples))
plt.figure()
plt.title(title + ", $\\mathcal{O}_{\\rm{CDW}}$")
O_CDW, O_CDW_eq = sim.O_CDW(True)
plt.plot(sim.time(), O_CDW, label=run_label)
if ref is not None:
	plt.plot(ref.time(), ref.O_CDW(), '--', label="ref")
plt.plot(sim.time()[[0,-1]], [O_CDW_eq, O_CDW_eq], 'k')
plt.plot(sim.time()[[0,-1]], [0, 0], 'k')

if show_inv_gamma_times and gamma != 0:
	for i in range(1, int(np.round(t_max*gamma))):
		plt.plot([i/gamma, i/gamma], [-1.1, 1.1], 'k')

plt.ylabel(r"$\mathcal{O}_{\rm{CDW}}$")
plt.xlabel("$t t_0$")
plt.xlim((0, sim.time()[-1]))
#plt.ylim((1.05*O_CDW.min(), 1.05*O_CDW.max()))
plt.legend()


# phonons
plt.figure()
plt.title(title + ", N$_{\\mathrm{ph}}$")
# TODO 
plt.plot(sim.time(), sim.nph())#, label=run_label)
plt.plot(sim.time(), sim.nph().sum(axis=1)/sim.num_electrons, 'k', label="N$_{ph}$/N")
if ref is not None:
	plt.gca().set_prop_cycle(None)
	plt.plot(ref.time(), ref.nph(), "--")#, label="ref")
	plt.plot(ref.time(), ref.nph().sum(axis=1)/ref.num_electrons, '--k', label="ref N$_{ph}$/N")
plt.ylabel(r"N$_{\mathrm{ph}}$")
plt.xlabel("$t t_0$")
plt.xlim((0, sim.time()[-1]))
#plt.legend()

# kinetic energy
plt.figure()
plt.title(title + ", kinetic energy")
plt.plot(sim.time(), sim.energy_kinetic(), label=run_label)
if ref is not None:
	try:
		plt.plot(ref.time(), ref.energy_kinetic(), '--', label="ref")
	except RuntimeError:
		pass
plt.xlabel("$t t_0$")
plt.xlim((0, sim.time()[-1]))
plt.ylabel(r"$\langle H_{el, kin}\rangle$")
plt.legend()

if show_energy_violation:
	# total energy
	plt.figure()
	plt.title("total energy - total energy(0) per site")
	energy_total = sim.energy_total()
	plt.plot(sim.time(), (energy_total - energy_total[0])/num_sites)


if num_sites == 2:
	# plots in relative q-space
	fig_osc_q, axs_osc_q = plt.subplots(1, 2, figsize=(10, 5))
	fig_osc_q.suptitle(title + ", q and pq")
	axs_osc_q[0].plot(sim.time(), sim.q())
	axs_osc_q[1].plot(sim.time(), sim.pq())
		
	if ref is not None:
		axs_osc_q[0].plot(ref.time(), ref.q(), '--')
		axs_osc_q[1].plot(ref.time(), ref.pq(), '--')
	
	if show_inv_gamma_times and gamma != 0:
		for i in range(1, int(np.round(t_max*gamma))):
			axs_osc_q[0].plot([i/gamma, i/gamma], [0, 1], 'k')
			axs_osc_q[1].plot([i/gamma, i/gamma], [0, 1], 'k')
	
	axs_osc_q[0].set_xlabel("time")
	axs_osc_q[0].set_ylabel("(q=x1-x0)/l0")
	axs_osc_q[1].set_xlabel("time")
	axs_osc_q[1].set_ylabel("(pq=p1-p0)*l0")
	#fig_osc_q.legend()
	fig_osc_q.tight_layout(rect=(0,0,1,.95))

	plt.figure()
	plt.title(title + ", pq-q phase space")
	plt.plot(sim.q(), sim.pq(), '-')
	if ref is not None:
		plt.plot(ref.q(), ref.pq(), "--", label=run_label)
	plt.xlabel("(q=x1-x0)/l0")
	plt.ylabel("(pq=p1-p0)*l0")
	
	
		
if show_compare_single_traj:
	# load miroslav_single_traj data
	# _me1 = miroslav ehrenfest single trajectory
	t_me1  = np.load("miroslav_single_traj/ehrenfest_single_time.npy")
	tn_me1 = np.load("miroslav_single_traj/ehrenfest_single_time_density.npy")
	q_me1  = np.load("miroslav_single_traj/ehrenfest_single_q.npy")
	pq_me1 = np.load("miroslav_single_traj/ehrenfest_single_pq.npy")
	n0_me1 = np.load("miroslav_single_traj/ehrenfest_single_n1.npy")
	n1_me1 = np.load("miroslav_single_traj/ehrenfest_single_n2.npy")

	# rescaling to my reference units
	q_me1 /= np.sqrt(2)
	pq_me1 *= 2*np.sqrt(2)
	
	cut_me1 = np.argmax(t_me1 > t_max)
	if cut_me1 == 0:
		cut_me1 = t_me1.size
		
		
	cut_n_me1 = np.argmax(tn_me1 > t_max)
	if cut_n_me1 == 0:
		cut_n_me1 = tn_me1.size
		
	plt.figure(2)
	plt.plot(tn_me1[:cut_n_me1], n0_me1[:cut_n_me1], '--', color="tab:green", label="me1 site 0")
	plt.plot(tn_me1[:cut_n_me1], n1_me1[:cut_n_me1], '--', color="tab:red",label="me1 site 1")
	plt.legend()
	
	plt.figure(7)
	plt.plot(tn_me1[:cut_n_me1], n1_me1[:cut_n_me1] - n0_me1[:cut_n_me1], '--', label="me1 ref")
	plt.legend()
	
	# q(t)
	plt.figure(8)
	plt.plot(t_me1[:cut_me1], q_me1[:cut_me1], '--', label="miro single ehrenfest")
	plt.legend()
	
	# pq-q phasespace
	plt.figure(9)
	plt.plot(q_me1[:cut_me1], pq_me1[:cut_me1], '--', label="miro single ehrenfest")
	plt.legend()
	
	if np.allclose(t[:cut_me1], t_me1[:cut_me1]):
		plt.figure()
		plt.title("q - q_me1")
		plt.plot(t_me1[:cut_me1], sim.q()[:cut_me1]-q_me1[:cut_me1])
		plt.xlabel("$ t t_0$")
		plt.ylabel("(q=x1-x0)/l0")
		
		plt.figure()
		plt.title("pq - pq_me1")
		plt.plot(t_me1[:cut_me1], sim.pq()[:cut_me1]-pq_me1[:cut_me1])
		plt.xlabel("$ t t_0$")
		plt.ylabel("(pq=p1-p0)*l0")
		
		plt.figure()
		plt.title("density difference ni - ni_me1")
		plt.plot(tn_me1[:cut_n_me1:2], sim.nel()[:, 0] - n1_me1[:cut_n_me1:2], label="site 0")
		plt.plot(tn_me1[:cut_n_me1:2], sim.nel()[:, 1] - n0_me1[:cut_n_me1:2], label="site 1")
		plt.xlabel("$ t t_0$")
		plt.ylabel("ni - ni_me1")
		plt.legend()

#plt.figure(7)
#plt.tight_layout()
#name = "figures/O_CDW_omega=%g_gamma=%g_sym%s_samples=%s.png"%(omega, gamma, symmetrize_start, num_samples_set)
#plt.savefig(name, dpi=500)

plt.show()
	
