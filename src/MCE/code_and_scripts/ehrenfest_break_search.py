from os import listdir, mkdir
from os.path import isfile, join
import re # regular expression parser
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

import data_manager


num_sites = 5
num_electrons = num_sites//2
initial_state = "BCDW"
psi_initial="local"

if num_sites == 2:
	path_to_ref = "ref_data/dimer_" + initial_state
	if initial_state != "BCDW" and initial_state != "DCDW":
		raise ValueError("initial state must be \"BCDW\" or \"DCDW\"")
	ref_file_or_folder = [f for f in listdir(path_to_ref) \
		                     if not isfile(join(path_to_ref, f))
		                 ]
	#holst_[freeFerm,pol]_L2M40tot50.0dt0.1t01.0gam-%.3fomg%.3fpb0/OBSL2M40t01.0omega%.1fgam-%.1fdt0.10tot50.PB0.bin
else:
	path_to_ref = "ref_data/DMRG/" + initial_state + "_omega=0.1_gamma=0.4/"
	prefix_dmrg13 = "Holstein_model_CDW_" + ("" if initial_state == "BCDW" else "dressed_") +"decay_"
	ref_file_or_folder = [f for f in listdir(path_to_ref)
		                     if isfile(join(path_to_ref, f)) and f[len(prefix_dmrg13):len(prefix_dmrg13)+4] == "time"]
	
	prefix_dmrg13 = path_to_ref + prefix_dmrg13
	
	phonon_suffix = "_truncOM1e-07_TEBDOM"+ ("" if initial_state == "BCDW" else "_li") +".txt"
	#Holstein_model_CDW_decay_time_L13_N6_t1_omega1_gamma1_tau0.025_time_steps260_bdim2000_trunc1e-07_cutph10_truncOM1e-07_TEBDOM.txt
	unique_params = np.unique([f[:-len(phonon_suffix)-2] for f in ref_file_or_folder])
	remaining_files = []
	for p in unique_params:
		fileset = [f for f in ref_file_or_folder if p in f]
		
		max_num_ph = 0
		for f in fileset:
			num_ph = int(f[-len(phonon_suffix)-2:-len(phonon_suffix)])
			if num_ph > max_num_ph:
				max_num_ph = num_ph
		remaining_files.append(fileset[0][:-len(phonon_suffix)-2] + str(max_num_ph) + phonon_suffix)
	
	ref_file_or_folder = remaining_files
	
	#suffix_dmrg13 = "_L13_N6_t1_omega%g_gamma%g_tau0.025_time_steps260_bdim2000_trunc1e-07_cutph%i_truncOM1e-07_TEBDOM.txt"%(omega, gamma, num_phonon_dmrg13)
#else:
#	raise ValueError("no reference data for given size and electron number existing")

ref_file_or_folder.sort()

print(len(ref_file_or_folder), "files found")

num_samples = 4000
symmetrize_start = False

eps_abs_O_CDW = 0.02
eps_abs_p = 0.01		# TODO p needs also needs a relative error

plot_data = True
plot_power_instead_var_skew_kurt = True
savePlots = True		# saved plots are not shown as window
printProgress = False
saveResultToFile = False
folder_prefix = ("moments_" if plot_power_instead_var_skew_kurt else "vsk_")\
                +initial_state + "_figs_compareL%i_num_samples=%i/"%(num_sites, num_samples)

selectParameterSet = False
selectedOmega = [2]
selectedGamma = [4]


plotQCoordinates = False and num_sites == 2
# whether higher moments of q shall be plotted. if false, number of phonons (per site) and kinetic energy are plotted
plotMomentsLarger2 = False and plotQCoordinates

global_t_max = 6.5
dt = 0.025

results = []

for f in ref_file_or_folder:
	if num_sites == 2:
		res = re.findall(r"\d*\.\d+|\d+", f)
		t_max_ref = float(res[2])
		#dt_ref = float(res[3])
		hopping = float(res[4])
		gamma = float(res[5])
		omega = float(res[6])
		
		
	else:
		res = re.findall(r"[-+]?\d*\.\d+|\d+", f)
		if num_sites != int(res[0]) or num_electrons != int(res[1]):
			continue
		
		hopping = float(res[2])
		omega = float(res[3])
		gamma = float(res[4])
	
	
	if selectParameterSet and (omega not in selectedOmega or gamma not in selectedGamma):
		continue
	
	ref = data_manager.RefData(initial_state, psi_initial, num_sites, omega, gamma)
	
	t_max = min(ref.time()[-1], global_t_max)

	sim = data_manager.Simulation(omega, gamma, t_max, dt,
	                              num_samples, num_sites=num_sites,
	                              symmetric_y_init=symmetrize_start,
	                              num_electrons=num_electrons,
	                              fullStorage=plotMomentsLarger2,
	                              printProgress=printProgress,
	                              initial_state=initial_state,
	                              psi_initial=psi_initial,
	                              calcYSquared=True)
	
	ref.set_cut_time(t_max)
	
	"""
	if num_sites == 2:
		suffix = "L2M40t01.0omega%.1fgam-%.1fdt0.10tot50.PB0.bin"%(omega, gamma)
		n0_ref = np.fromfile(path_to_ref+f + "/n0"+suffix, dtype=float)[:cut_ref]
		n1_ref = np.fromfile(path_to_ref+f + "/n1"+suffix, dtype=float)[:cut_ref]
		density_ref = np.empty((n0_ref.size, 2))
		density_ref[:, 0] = n0_ref
		density_ref[:, 1] = n1_ref
	elif num_sites == 13:
		density_ref = np.loadtxt(prefix_dmrg13 + "density_fermions" + suffix_dmrg13)[:cut_ref]
		#n_ph_ref = np.loadtxt(prefix_dmrg13 + "density_phonons" + suffix_dmrg13)[:cut_ref]
		#E_kin_ref = np.loadtxt(prefix_dmrg13 + "kinetic_energy" + suffix_dmrg13)[:cut_ref]
		#x_ref = np.loadtxt(prefix_dmrg13 + "lattice_displacement" + suffix_dmrg13)[:cut_ref]
		#p_ref = np.loadtxt(prefix_dmrg13 + "lattice_momentum" + suffix_dmrg13)[:cut_ref]
		#nk_ref = np.loadtxt(prefix_dmrg13 + "nk" + suffix_dmrg13)[:cut_ref]			# electron density in momentum space
		#nl_nm_ref = np.loadtxt(prefix_dmrg13 + "nl_nm" + suffix_dmrg13)[:cut_ref]		# electron density correlation between mid and site l>= L//2
	"""
	O_CDW = sim.O_CDW()
	O_CDW_ref = ref.O_CDW()
	
	ind_breaked_density = np.argmax(np.abs(O_CDW_ref-O_CDW) > eps_abs_O_CDW)
	if ind_breaked_density == 0:	# not breaked within simulation time
		t_break = -sim.time()[-1]
	else:
		t_break = sim.time()[ind_breaked_density]
	
	print(omega, gamma, t_break)
	results.append((omega, gamma, t_break))
	
	
	if plot_data: # TODO make it also runnable for 13 (and other?) sizes
		fig, axs = plt.subplots(3 if plotQCoordinates else 2,3, figsize=(16, 9))
	
		fig.suptitle("%s: L=%i, $\\omega=%g$, $\\gamma=%g$"%(initial_state, num_sites, omega, gamma))
		
		ax_OCDW = axs[0, 0]
		ax_OCDW.set_xlabel("t*t0")
		ax_OCDW.set_ylabel(r"$\mathcal{O}_{\rm{CDW}}$")
		
		ax_density = axs[0, 1]
		ax_density.set_xlabel("t*t0")
		ax_density.set_ylabel(r"$\langle n \rangle$")
		
		ax_px = axs[0, 2]
		ax_diff_p = None
		ax_xt = None
		ax_pt = None
		
		ax_qt = axs[1, 0] if plotQCoordinates else None
		ax_pqt = axs[1, 1] if plotQCoordinates else None
		ax_pqq = axs[1, 2] if plotQCoordinates else None
		ax_var_qpq = axs[2,0] if plotQCoordinates else None
		
		ax_skew_qpq = axs[2,1] if plotMomentsLarger2 else None
		ax_kurtosis_qpq = axs[2, 2] if plotMomentsLarger2 else None
		
		ax_nph = axs[-1, 1] if not plotMomentsLarger2 else None
		ax_ekin = axs[-1, 2] if not plotMomentsLarger2 else None
		
		
		t_break = abs(t_break)
		
		# O_CDW(t)
		ax_OCDW.cla()
		ax_OCDW.plot(sim.time(), O_CDW, label="ehrenfest")
		ax_OCDW.plot(sim.time(), O_CDW_ref, '--', label="ref")
		ax_OCDW.plot([t_break]*2, [O_CDW.min(), 1], 'k')
		ax_OCDW.set_xlabel("t*t0")
		ax_OCDW.set_ylabel(r"$\mathcal{O}_{\rm{CDW}}$")
		ax_OCDW.legend()
		
		interp_nel = interp1d(ref.time(), ref.nel(), axis= 0, kind=3, copy=False, assume_sorted=True)
		interp_x = interp1d(ref.time(), ref.x(), axis= 0, kind=3, copy=False, assume_sorted=True)
		interp_p = interp1d(ref.time(), ref.p(), axis= 0, kind=3, copy=False, assume_sorted=True)
		
		
		# n(t)
		ax_density.cla()
		lines = ax_density.plot(sim.time(), sim.nel())
		ax_density.plot(sim.time(), interp_nel(sim.time()), '--')
		ax_density.plot([t_break]*2, [0, 1], 'k')
		#ax_density.set_xlim((0, 1.5*t_break))
		ax_density.set_xlabel("t*t0")
		ax_density.set_ylabel(r"$\langle n \rangle$")
		ax_density.legend(lines, ("site 0", "site 1"), loc="lower right")
		
		
		# x(t)
		if ax_xt is not None:
			ax_xt.cla()
			#plt.title("x")
			lines = ax_xt.plot(sim.time(), sim.x())
			ax_xt.plot(sim.time(), interp_x(sim.time()), '--')
			ax_xt.plot([t_break, t_break], [sim.x().min(), sim.x().max()], 'k')
			#ax_xt.set_xlim((0, 1.5*t_break))
			ax_xt.set_xlabel("t*t0")
			ax_xt.set_ylabel("x/l0")
			#ax_xt.legend(lines, ("site 0", "site 1"))

		
		# p(t)
		if ax_pt is not None:
			ax_pt.cla()
			lines = ax_pt.plot(sim.time(), sim.p())
			ax_pt.plot(sim.time(), interp_p(sim.time()), '--')
			ax_pt.plot([t_break, t_break], [sim.p().min(), sim.p().max()], 'k')
			#ax_pt.set_xlim((0, 1.5*t_break))
			
			ax_pt.set_xlabel("t*t0")
			ax_pt.set_ylabel("p*l0")
			ax_pt.legend()
		
		
		# p(x)
		if ax_px is not None:
			ax_px.cla()
			ax_px.plot(sim.x(), sim.p())
			ax_px.plot(interp_x(sim.time()), interp_p(sim.time()), '--')
			
			ax_px.set_xlabel("x/l0")
			ax_px.set_ylabel("p*l0")
			#ax_px.legend()
		
		# p_diff(t)
		if ax_diff_p is not None:
			ax_diff_p.cla()
			ax_diff_p.plot(sim.time(), sim.p()-interp_p(sim.time()))
			ax_diff_p.plot([t_break]*2, [-eps_abs_p, eps_abs_p], 'k')
			ax_diff_p.set_xlabel("t*t0")
			ax_diff_p.set_ylabel("diff p = p*l0-p$_\\mathrm{ref}$")
			#ax_diff_p.legend()
		
		
		if num_sites == 2:			
			interp_q1 = interp1d(ref.time(), ref.q(), kind=3, copy=False, assume_sorted=True)
			interp_q2 = interp1d(ref.time(), ref.q2(), kind=3, copy=False, assume_sorted=True)
			interp_q3 = interp1d(ref.time(), ref.q3(), kind=3, copy=False, assume_sorted=True)
			interp_q4 = interp1d(ref.time(), ref.q4(), kind=3, copy=False, assume_sorted=True)
			
			interp_pq1 = interp1d(ref.time(), ref.pq(), kind=3, copy=False, assume_sorted=True)
			interp_pq2 = interp1d(ref.time(), ref.pq2(), kind=3, copy=False, assume_sorted=True)
			interp_pq3 = interp1d(ref.time(), ref.pq3(), kind=3, copy=False, assume_sorted=True)
			interp_pq4 = interp1d(ref.time(), ref.pq4(), kind=3, copy=False, assume_sorted=True)
		
			if ax_qt is not None:
				ax_qt.cla()
				ax_qt.plot(sim.time(), sim.q())
				#ax_qt.plot(ref.time(), q1_ref, "--")
				ax_qt.plot(sim.time(), interp_q1(sim.time()), "--")
				ax_qt.set_xlabel("time t*t0")
				ax_qt.set_ylabel("(q = x0-x1)/l0")
			
			if ax_pqt is not None:
				ax_pqt.cla()
				ax_pqt.plot(sim.time(), sim.pq())
				#ax_pqt.plot(ref.time(), pq1_ref, "--")
				ax_pqt.plot(sim.time(), interp_pq1(sim.time()), "--")
				ax_pqt.set_xlabel("time t*t0")
				ax_pqt.set_ylabel("(pq = p0-p1)*l0")
			
			if ax_pqq is not None:
				ax_pqq.cla()
				ax_pqq.plot(sim.q(), sim.pq())
				ax_pqq.plot(interp_q1(sim.time()), interp_pq1(sim.time()), "--")
				ax_pqq.set_xlabel("(q = x0-x1)/l0")
				ax_pqq.set_ylabel("(pq = p0-p1)*l0")
			
			var_q_ref  = ref.q2() - ref.q()**2
			var_pq_ref = ref.pq2() - ref.pq()**2
			
			if ax_var_qpq is not None:
				ax_var_qpq.cla()
				if plot_power_instead_var_skew_kurt:
					ax_var_qpq.plot(sim.time(), sim.q2(), label="q")
					ax_var_qpq.plot(sim.time(), sim.pq2(), label="p_q")
					ax_var_qpq.plot(ref.time(), ref.q2(), "--")
					ax_var_qpq.plot(ref.time(), ref.pq2(), "--")
					ax_var_qpq.set_ylabel("mean q^2 or pq^2")
				elif num_samples >= 4:
					ax_var_qpq.plot(sim.time(), sim.get_q_variance())
					ax_var_qpq.plot(sim.time(), sim.get_p_q_variance())
					ax_var_qpq.plot(ref.time(),  var_q_ref, "--")
					ax_var_qpq.plot(ref.time(), var_pq_ref, "--")
					ax_var_qpq.set_ylabel("variance q/l0 or pq*l0")
				
				ax_var_qpq.set_xlabel("time t*t0")
			
			
			if ax_skew_qpq is not None:
				ax_skew_qpq.cla()
				if plot_power_instead_var_skew_kurt:
					ax_skew_qpq.plot(sim.time(), (sim.get_q_full()**3).mean(axis=1), label="q")
					ax_skew_qpq.plot(sim.time(), (sim.get_p_q_full()**3).mean(axis=1), label="p_q")
					ax_skew_qpq.plot(ref.time(),  ref.q3(), "--")
					ax_skew_qpq.plot(ref.time(), ref.pq3(), "--")
					ax_skew_qpq.set_ylabel("mean q^3 or pq^3")
				elif num_samples >= 4:
					ax_skew_qpq.plot(sim.time(), sim.get_q_skew())
					ax_skew_qpq.plot(sim.time(), sim.get_p_q_skew())
					ax_skew_qpq.plot(ref.time(), ( ref.q3() - 3* ref.q()* ref.q2() + 2* ref.q()**3)/var_q_ref**1.5, "--")
					ax_skew_qpq.plot(ref.time(), (ref.pq3() - 3*ref.pq()*ref.pq2() + 2*ref.pq()**3)/var_pq_ref**1.5, "--")
					ax_skew_qpq.set_ylabel("skew q*l0 or pq/l0")
				
				ax_skew_qpq.set_xlabel("time t*t0")
			
			if ax_kurtosis_qpq is not None:
				ax_kurtosis_qpq.cla()
				if plot_power_instead_var_skew_kurt:
					ax_kurtosis_qpq.plot(sim.time(), sim.q4(), label="q")
					ax_kurtosis_qpq.plot(sim.time(), sim.pq4(), label="p_q")
					ax_kurtosis_qpq.plot(ref.time(), ref.q4(), "--")
					ax_kurtosis_qpq.plot(ref.time(), ref.pq4(), "--")
					ax_kurtosis_qpq.set_ylabel("mean q^4 or pq^4")
				elif num_samples >= 4:
					ax_kurtosis_qpq.plot(sim.time(), sim.get_q_kurtosis(), label="q")
					ax_kurtosis_qpq.plot(sim.time(), sim.get_p_q_kurtosis(), label="p_q")
					ax_kurtosis_qpq.plot(ref.time(), ( ref.q4() - 4* ref.q()* ref.q3() + 6* ref.q()**2* ref.q2()-3* ref.q()**4)/(var_q_ref**2) - 3, "--")
					ax_kurtosis_qpq.plot(ref.time(), (ref.pq4() - 4*ref.pq()*ref.pq3() + 6*ref.pq()**2*ref.pq2()-3*ref.pq()**4)/(var_pq_ref**2) - 3, "--")
					ax_kurtosis_qpq.set_ylabel("excess kurtosis of q*l0 or p_q/l0")
				
				#ax_kurtosis_qpq.plot([t_break]*2, [sim.get_q_kurtosis().min(), sim.get_q_kurtosis().max()], 'k')
				ax_kurtosis_qpq.set_xlabel("time t*t0")
				ax_kurtosis_qpq.legend()
		
		if ax_nph is not None:
			lines = ax_nph.plot(sim.time(), sim.nph())
			ax_nph.plot(sim.time(), sim.nph().sum(axis=1)/num_electrons, 'k')
			interp_nph = interp1d(ref.time(), ref.nph(), axis=0, kind=3, copy=False, assume_sorted=True)
			ax_nph.plot(sim.time(), interp_nph(sim.time()), '--')
			ax_nph.plot(sim.time(), interp_nph(sim.time()).sum(axis=1)/num_electrons, '--k')
			ax_nph.legend(lines, ("site 0", "site 1"))
			ax_nph.set_xlabel("time")
			ax_nph.set_ylabel(r"$\langle N_{ph}\rangle$")
			
		if ax_ekin is not None:
			lines = ax_ekin.plot(sim.time(), sim.energy_kinetic())
			interp_ekin = interp1d(ref.time(), ref.energy_kinetic(), kind=3, copy=False, assume_sorted=True)
			ax_ekin.plot(sim.time(), interp_ekin(sim.time()), '--')
			ax_ekin.legend(lines, ("site 0", "site 1"))
			ax_ekin.set_xlabel("time")
			ax_ekin.set_ylabel(r"$\langle H_{kin}\rangle$")
		
		
		fig.tight_layout()
		if savePlots:
			try:
				mkdir(folder_prefix)
			except FileExistsError:
				pass
			try:
				mkdir(folder_prefix+"/gamma=%.2f"%(gamma))
			except FileExistsError:
				pass
			fig.savefig(folder_prefix+"omega=%.2f_gamma=%.2f.png"%(omega, gamma), dpi=300)
			fig.savefig(folder_prefix+"gamma=%.2f/omega=%.2f_gamma=%.2f.png"%(gamma, omega, gamma), dpi=300)
		else:
			plt.show()
		plt.close(fig)
		
if saveResultToFile:
	np.savetxt(folder_prefix+initial_state+"_ehrenfest_L%i_break_times_num_samples=%i.dat"%(num_sites, num_samples),
		       np.array(results),
		       header="#omega gamma t_break")
