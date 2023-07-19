import numpy as np
import os
import re

class DataManager(object):
	def __init__(self):
		if type(self) is DataManager:
			raise Exception("DataManager is an abstract class and cannot be instantiated directly")
		self._export_path = "data_export/"
	
	
	def get_cut_index(self, t_max):
		if self._time[-1] <= t_max:
			cut = self._time.size
		else:
			cut = np.argmax(self._time > t_max)
		return cut
	
	
	def set_cut_time(self, t_max):
		self.cut = self.get_cut_index(t_max)
	
	
	def reset_cut_time(self):
		self.cut = self._time.size
	
	
	def time(self, apply_cut=True, BO_steps=False):
		if apply_cut:
			retVal = self._time[:self.cut]
		else:
			retVal = self._time
		
		if BO_steps and hasattr(self, "BO_nel_save_dist"):
			retVal = retVal[::self.BO_nel_save_dist]
		return retVal
	
	
	def nel(self, site_x=None, site_y=None, squeezeY=True, apply_cut=True):
		return self._site_selection(self.obs_nel, site_x, site_y, squeezeY, apply_cut)
	
	
	def nel_BO(self, state=None, apply_cut=True):
		if self.obs_nel_BO is None:
			raise RuntimeError("Born Oppenheimer occupation number wasn't calculated")
		
		if state is None:
			retVal =  self.obs_nel_BO[:, :]
		else:
			retVal =  self.obs_nel_BO[:, state]
		
		if apply_cut:
			if hasattr(self, "BO_nel_save_dist"):
				return retVal[:1+(self.cut-1)//self.BO_nel_save_dist]
			else:
				return retVal[:self.cut]
		else:
			return retVal
	
	
	#FIXME adopt to changed interface
	def O_CDW(self, retEqVal=False, ensure_pos_start=True):
		o_cdw = np.zeros(self.time().size)
		if not hasattr(self, "num_sites_y") or self.num_sites_y == 1:
			for i in range(self.num_sites):
				if i%2 == 1:
					o_cdw += self.nel(i)
				else:
					o_cdw -= self.nel(i)
			o_cdw_eq = (self.num_sites%2)/self.num_sites
		else:
			for i in range(self.num_sites_y):
				for i in range(self.num_sites_x):
					if (i+j)%2 == 0:
						o_cdw += self.nel(i, j)
					else:
						o_cdw -= self.nel(i, j)
			o_cdw_eq = ((self.num_sites_x*self.num_sites_y)%2)/(self.num_sites_x*self.num_sites_y)
		o_cdw /= self.num_electrons
		if ensure_pos_start and o_cdw[0] < 0:
			o_cdw = -o_cdw
			o_cdw_eq = -o_cdw_eq
		if retEqVal:
			return o_cdw, o_cdw_eq
		else:
			return o_cdw
	
	
	def nph(self, site_x=None, site_y=None, squeezeY=True, apply_cut=True):
		return self._site_selection(self.obs_nph, site_x, site_y, squeezeY, apply_cut)
	
	
	def x(self, site_x=None, site_y=None, squeezeY=True, apply_cut=True):
		return self._site_selection(self.obs_y[:, :, :, 0], site_x, site_y, squeezeY, apply_cut)
	
	
	def p(self, site_x=None, site_y=None, squeezeY=True, apply_cut=True):
		return self._site_selection(self.obs_y[:, :, :, 1], site_x, site_y, squeezeY, apply_cut)
	
	
	def energy_kinetic(self, electron="sum", scaled=True, apply_cut=True):
		if self.obs_el_ekin is None:
			raise RuntimeError("kinetic energy wasn't stored in calculation")
		if electron == "sum":
			retVal = self.obs_el_ekin[:, :].sum(axis=-1)
		elif electron == "all":
			retVal =  self.obs_el_ekin[:, :]
		else:
			retVal =  self.obs_el_ekin[:, electron]
		
		if apply_cut:
			retVal = retVal[:self.cut]
		
		if scaled:
			return retVal*self.hopping
		else:
			return retVal
	
	# TODO doesn't work for ref_data
	def energy_interaction(self, site="sum", scaled=True, apply_cut=True):
		if self.obs_el_ph_inter is None:
			raise RuntimeError("interaction energy wasn't stored in calculation")
		if site == "sum":
			retVal = self.obs_el_ph_inter[:].sum(axis=(-2, -1))
		elif site == "all":
			retVal =  self.obs_el_ph_inter[:]
		else:
			retVal =  self.obs_el_ph_inter[:]
		
		if apply_cut:
			retVal = retVal[:self.cut]
			
		if scaled:
			return retVal*self.gamma
		else:
			return retVal
	
	
	def energy_total(self, apply_cut=True):
		return   self.energy_interaction("sum", True, apply_cut)\
		       + self.energy_kinetic("sum", True, apply_cut)\
		       + self.omega*self.nph(squeezeY=False, apply_cut=apply_cut).sum(axis=(1, 2))
	
	
	def q(self):
		if self.num_sites != 2:
			raise RuntimeError("q coordinate undefined because num_sites == %i != 2"%self.num_sites)
		if not hasattr(self, "phononCoordQ"):
			self.phononCoordQ = self.x(0) - self.x(1)
		return self.phononCoordQ[:self.cut]
	
	
	def pq(self):
		if self.num_sites != 2:
			raise RuntimeError("pq coordinate undefined because num_sites == %i != 2"%self.num_sites)
		if not hasattr(self, "phononCoordPQ"):
			self.phononCoordPQ = self.p(0) - self.p(1)
		return self.phononCoordPQ[:self.cut]
	
	
	# shapes observable into desired shape with subselection
	# obs must have shape (self.num_times, self.num_sites_x, self.num_sites_y)
	def _site_selection(self, obs, site_x, site_y, squeezeY, apply_cut):
		if site_x is not None:
			retVal = obs[:, :, site_x]
		else:
			retVal = obs[:, :, :]
		
		if squeezeY and self.num_sites_y == 1:
			retVal = retVal[:, 0]
		elif site_y is not None:
			retVal = retVal[:, site_y]
		
		if apply_cut:
			retVal = retVal[:self.cut]
		
		return retVal



# TODO check interface
class RefData(DataManager):
	"""class to coordinate loading reference data for Holstein model from files
	
	the class with do it's bets to etsimate whether a dataset from David (ED) or
	Jan (DMRG) should be used, based on the number of sites
	"""
	
	def __init__(self,
	             initial_state,
	             psi_initial,
	             psi_state,
	             num_sites_x,
	             num_sites_y,
	             hopping,
	             omega,
	             gamma):
		"""creates data container with the requested dataset from files or might fail if data not found
		
		
		if there are 2 sites, ED data are used, otherwise DMRG is used
		
		
		initial_state
		state at which the system started.
		must be either bare or dressed
		
		num_sites (= num_sites_x)
		number of sites of the system (only 1d, implicit num_sites_y = 1)
		
		omega
		phonon frequency
		
		gamma
		coupling strength
		"""
		super(RefData, self).__init__()
		
		if initial_state != "bare" and initial_state != "dressed":
			raise ValueError("initial_state must be either bare or dressed")
		if psi_initial not in ["local", "eigen", "CDW"]:
			raise ValueError("reference data only exist for local and eigen psi_initial")
		
		
		self.num_sites = num_sites_x*num_sites_y
		if psi_initial =="CDW":
			self.num_electrons = self.num_sites//2
		else:
			self.num_electrons = 1
		# right now fixed and no support for DMRG2 folder
		self.initial_state = initial_state
		# rescale all parameters to hopping = 1 for path construction and later change time and observables
		self.omega = omega/hopping
		self.gamma = gamma/hopping
		self.hopping = 1.
		
		def parse_hopping(f):
			f_hopping = re.findall(r"t0(\d+)\.(\d*)gam", f)[0]
			return float(f_hopping[0] + "." + f_hopping[1])
		
		def parse_gamma(f):
			f_gamma = re.findall(r"gam-?(\d+)\.(\d*)", f)[0]
			return float(f_gamma[0] + "." + f_gamma[1])
		
		def parse_omega(f):
			f_omega = re.findall(r"omg(\d+)\.(\d*)", f)[0]
			return float(f_omega[0] + "." + f_omega[1])
		
		def filter_folder_for_params(folder_list):
			folder_list_corr_params = []
			
			#folder_list = np.sort(folder_list)
			#print(self.hopping, self.gamma, self.omega)
			for f in folder_list:
				#print(f)
				try:
					f_hopping = parse_hopping(f)
					f_gamma = parse_gamma(f)
					f_omega = parse_omega(f)
					#print(f_hopping, f_gamma, f_omega)
					
					if np.allclose([f_hopping/f_omega, f_gamma/f_omega],
					               [self.hopping/self.omega, self.gamma/self.omega]):
						folder_list_corr_params.append(f)
					#	print("correct params")
					#else:
					#	print("wrong params")
				except Exception as e:
					print("failed to parse", f, "in folder", self.filePrefix)
					print(e)
			return folder_list_corr_params
		
		self.filePrefix = "ref_data/"
		
		if psi_initial == "eigen":
			self.dataSource = "MTE_Michael"
			self.filePrefix = "ref_data/eigen/MTE_%s_psi_initial=%s_hopping=%g_omega=%g_gamma=%g_num_sites_x=%i_num_sites_y=1/"\
			                  %(initial_state, psi_initial + str(psi_state) if psi_initial == "eigen" else psi_initial,\
			                    hopping, omega, gamma, self.num_sites)
			self.fileSuffix = ""
			data = np.loadtxt(self.filePrefix + "electron_density")
			self._time = data[:, 0]
			self.electronDensity = data[:, 1:]
			self.phononCoordX = np.loadtxt(self.filePrefix + "x")[:, 1:]
		else:
			if self.num_sites in [2,3]:
				self.dataSource = "ED"
				
				
				# build filenames
				if self.num_sites == 2:
					self.filePrefix += "dimer_" + self.initial_state + "/"
				else: # three sites
					self.filePrefix += "trimer_exact/"
				
				pattern = "holst_"
				if self.initial_state == "bare":
					pattern += "freeFerm"
				elif self.initial_state == "dressed":
					pattern += "pol"
				pattern += f"_L{self.num_sites}"
				
				if self.num_sites == 3:
					initial_site = psi_state[0]*num_sites_x + psi_state[1]
					pattern += f"site_in{initial_site}"
				
				folder_list = [s for s in os.listdir(self.filePrefix) if len(re.findall(pattern, s)) != 0]
				
				folder_list_corr_params = filter_folder_for_params(folder_list)
				
				best_folder = ""
				for f in folder_list_corr_params:
					#print(f)
					params = re.findall(r"M(\d+)tot(\d+\.\d+)dt(\d+\.\d+)", f)
					M = int(params[0][0])
					if best_folder == "" or M > best_M:
						best_folder = f
						best_M = M
						t_max = float(params[0][1])
						dt = float(params[0][2])
				
				self.filePrefix += best_folder + "/"
				
				self.omega    = parse_omega(best_folder)
				self.gamma    = parse_gamma(best_folder)
				self.hopping  = parse_hopping(best_folder)
				
				self.energy_rescale_fact = self.omega/omega
				#print(f"{best_folder=}, {t_max=}, {dt=}")
				#print(f"{self.energy_rescale_fact=}")
				
				self._time = np.linspace(0, t_max, int(np.round(t_max/dt)), endpoint=False)
				
				self.fileSuffix = ""	# this value is ignored
				
			elif psi_initial =="CDW":
				self.dataSource = "DMRG"
				use_new_files = (self.num_sites != 13 or self.num_electrons != 6) \
				                and self.omega == 0.1 and self.gamma == 0.4 \
				                or self.initial_state != "bare"
				
				self.filePrefix += "DMRG/" + ("BCDW" if self.initial_state == "bare" else "DCDW")
				if use_new_files:
					self.filePrefix += "_omega=%g_gamma=%g/"%(self.omega, self.gamma)
				else:
					self.filePrefix += "_melting/"
				
				self.filePrefix += "Holstein_model_CDW_"
				if use_new_files and self.initial_state == "dressed":
					self.filePrefix += "dressed_"
				self.filePrefix += "decay_"
				
				dmrg_num_phonons = self._get_dmrg_num_phonons(use_new_files)
				if dmrg_num_phonons == 0:
					raise RuntimeError("DMRG data NOT FOUND for num_sites=%i, num_electrons=%i,"\
					                   " omega=%g, gamma=%g in %s files"\
					                   %(self.num_sites, self.num_electrons, self.omega, self.gamma,
					                     "new" if use_new_files else "old"))
				
				if use_new_files and self.initial_state == "dressed" and self.omega == 2. and (self.gamma in [3., 4.]):
					trunc_weight = "1e-08"
				else:
					trunc_weight = "1e-07"
				self.fileSuffix = "_L%i_N%i_t%g_omega%g_gamma%g_tau0.025_time_steps260_bdim2000_trunc%s_cutph%i_truncOM%s_TEBDOM"\
				                  %(self.num_sites,
				                    self.num_electrons,
				                    self.hopping,
				                    self.omega,
				                    self.gamma,
				                    trunc_weight,
				                    dmrg_num_phonons,
				                    trunc_weight)
				
				if self.omega == 0.1 and self.initial_state == "dressed":
					self.fileSuffix += "_li"
				self.fileSuffix += ".txt"
				self._time = np.loadtxt(self.filePrefix + "time" + self.fileSuffix)/self.hopping
				
				# further stored observables:
				#nk_dmrg13 = np.loadtxt(self.filePrefix + "nk" + self.fileSuffix)			# electron density in momentum space
				#nl_nm_dmrg13 = np.loadtxt(self.filePrefix + "nl_nm" + self.fileSuffix)		# electron density correlation between mid and site l>= L//2
			elif psi_initial == "local":
				self.dataSource = "DMRG2"
				
				self.filePrefix += "DMRG2/"
				
				param_string = "t0%.1fgam%.1fomg%.1fdt0.01"%(hopping, gamma, omega)
				
				if num_sites_y != 1:
					M = 20 if num_sites_x ==5 else 30
					initial_site = psi_state[0]*num_sites_x + psi_state[1] + 1
					self.filePrefix += "2dholtdvplbo_single_ebcdwLx{num_sites_x}Ly{num_sites_y}"
					self.filePrefix += f"M{M}site_in{initial_site}Md3000lboMd100" + param_string + "dt0.01tot20.0lbocut1e-09cut1e-09DIR/"
				
				else:
					# parse holtdmrglbo_single_e.cdwL.*M.*site_in.*Md3000lboMd100t0.*gam.*omg.*dt0.01tot20.0lbocut.*cut.*DIR
					
					pattern = "holtdmrglbo_single_e"
					if self.initial_state == "bare":
						pattern += "b"
					else:
						pattern += "d"
					initial_site = psi_state[1]+1
					pattern += f"cdwL{num_sites_x}M(\\d+)site_in{initial_site}Md3000lboMd100"
					
					folder_list = [s for s in os.listdir(self.filePrefix) if len(re.findall(pattern, s)) != 0]
					
					folder_list_corr_params = filter_folder_for_params(folder_list)
					
					best_folder = ""
					
					for f in folder_list_corr_params:
						M  = int(re.findall(r"M(\d+)", f)[0])
						Md = int(re.findall(r"Md(\d+)", f)[0])
						lboMd = int(re.findall(r"lboMd(\d+)", f)[0])
						lbocut = float(re.findall(r"lbocut(.+)cut(.+)DIR", f)[0][0])
						cut = float(re.findall(r"lbocut(.+)cut(.+)DIR", f)[0][1])
						
						if    best_folder == ""\
						   or M == 40 and Md >= best_Md and lboMd >= best_lboMd\
						      and lbocut <= best_lbocut and cut <= best_cut:
							best_folder = f
							best_M = M
							best_Md = Md
							best_lboMd = lboMd
							best_lbocut = lbocut
							best_cut = cut
					
					self.omega    = parse_omega(best_folder)
					self.gamma    = parse_gamma(best_folder)
					self.hopping  = parse_hopping(best_folder)
					
					self.filePrefix += best_folder + "/"
					#print(self.filePrefix)
				
				try:
					self._time = np.loadtxt(self.filePrefix + "time.dat")
				except OSError:
					num_data_points = np.loadtxt(self.filePrefix + "nlocSite1.dat").size
					self._time = np.arange(num_data_points)*0.01
				self.electronDensity = np.empty((self._time.size, num_sites_y*num_sites_x))
				self.phononDensity = np.empty((self._time.size, num_sites_y*num_sites_x))
				
				for y in range(num_sites_y):
					for x in range(num_sites_x):
						ind = y*num_sites_x + x
						self.electronDensity[:, ind] = np.loadtxt(self.filePrefix +   f"nlocSite{ind+1}.dat")
						self.phononDensity  [:, ind] = np.loadtxt(self.filePrefix + f"NphlocSite{ind+1}.dat")
				
				if num_sites_y == 1:
					self.electronDensity = self.electronDensity.reshape(self._time.size, num_sites_x)
					self.phononDensity   = self.phononDensity  .reshape(self._time.size, num_sites_x)
				
				
		if not hasattr(self, "_time"):
			raise ValueError("No reference data found")
		# rescale time to correct parameters
		self._time *= self.omega/omega
		
		# correct simulation parameters to the actual requested
		self.omega = omega
		self.gamma = gamma
		self.hopping = hopping
		
		# additionally, energy_kinetic must be rescaled. happens at loading in that function
		
		self.reset_cut_time()
		# load any data set to ensure that the path is valid
		self.nel()
	
	
	def _load_ed_file(self, obs_name):
		filenames = os.listdir(self.filePrefix)
		#print(filenames)
		pattern = r"\A" + obs_name + r"\.*"
		#print(pattern)
		fname = [f for f in filenames if len(re.findall(pattern, f)) != 0][0]
		#print(fname, obs_name)
		return np.fromfile(self.filePrefix + fname, dtype=float)
	
	
	def nel(self, site=None, apply_cut=True):
		if not hasattr(self, "electronDensity"):
			try:
				if self.dataSource == "ED":
					self.electronDensity = np.empty((self._time.size, self.num_sites))
					for i in range(self.num_sites):
						ni = self._load_ed_file("n%i"%i)
						self.electronDensity[:, i] = ni
				elif self.dataSource == "DMRG":
					self.electronDensity = np.loadtxt(self.filePrefix + "density_fermions" + self.fileSuffix)
				else:
					raise FileNotFoundError("unkown data source")
			except FileNotFoundError as e:
				print("failed to load electronDensity:" + str(e))
				self.electronDensity = np.full((self._time.size, self.num_sites), np.nan)
		
		if site is not None:
			retVal =  self.electronDensity[:, site]
		else:
			retVal =  self.electronDensity[:]
		
		if apply_cut:
			return retVal[:self.cut]
		else:
			return retVal
	
	
	def energy_kinetic(self, electron="sum", scaled=True, apply_cut=True):
		if electron != "sum" or not scaled:
			raise ValueError("kinetic energy is only provided summed over electrons and scaled in reference data")
		if not hasattr(self, "el_kinetic_energy"):
			try:
				if self.dataSource == "ED":
					self.el_kinetic_energy = self._load_ed_file("ek")/self.hopping/self.energy_rescale_fact
				elif self.dataSource == "DMRG":
					self.el_kinetic_energy = np.loadtxt(self.filePrefix + "kinetic_energy" + self.fileSuffix)
				else:
					raise FileNotFoundError("unkown data source")
				self.el_kinetic_energy *= self.hopping
			except FileNotFoundError as e:
				print("failed to load kinetic energy:" + str(e))
				self.el_kinetic_energy = np.full(self._time.size, np.nan)
		
		if apply_cut:
			return self.el_kinetic_energy[:self.cut]
		else:
			return self.el_kinetic_energy
	
	
	def nph(self, site=None, apply_cut=True):
		if not hasattr(self, "phononDensity"):
			try:
				if self.dataSource == "ED":
					self.phononDensity = np.empty((self._time.size, self.num_sites))
					for i in range(self.num_sites):
						ni = self._load_ed_file("nph%i"%i)
						self.phononDensity[:, i] = ni
				elif self.dataSource == "DMRG":
					self.phononDensity = np.loadtxt(self.filePrefix + "density_phonons" + self.fileSuffix)
				else:
					raise FileNotFoundError("unkown data source")
			except FileNotFoundError as e:
				print("failed to load phononDensity:" + str(e))
				self.phononDensity = np.full((self._time.size, self.num_sites), np.nan)
		if site is not None:
			retVal = self.phononDensity[:, site]
		else:
			retVal = self.phononDensity[:]
		
		if apply_cut:
			return retVal[:self.cut]
		else:
			return retVal
	
	
	def x(self, site=None, apply_cut=True):
		if not hasattr(self, "phononCoordX"):
			try:
				if self.dataSource == "ED":
					self.phononCoordX = np.empty((self._time.size, self.num_sites))
					for i in range(self.num_sites):
						xi = self._load_ed_file("x%i"%i)
						self.phononCoordX[:, i] = xi
				elif self.dataSource == "DMRG":
					self.phononCoordX = np.loadtxt(self.filePrefix + "lattice_displacement" + self.fileSuffix)/np.sqrt(2)
				else:
					raise FileNotFoundError("unkown data source")
			except FileNotFoundError as e:
				print("failed to load phononCoordX:" + str(e))
				self.phononCoordX = np.full((self._time.size, self.num_sites), np.nan)
		if site is not None:
			retVal = self.phononCoordX[:self.cut, site]
		else:
			retVal = self.phononCoordX[:self.cut]
		
		if apply_cut:
			return retVal[:self.cut]
		else:
			return retVal
	
	
	def p(self, site=None, apply_cut=True):
		if not hasattr(self, "phononCoordP"):
			try:
				if self.dataSource == "ED":
					self.phononCoordP = np.empty((self._time.size, self.num_sites))
					for i in range(self.num_sites):
						ni = self._load_ed_file("p%i"%i)
						self.phononCoordP[:, i] = ni
				elif self.dataSource == "DMRG":
					self.phononCoordP = np.loadtxt(self.filePrefix + "lattice_momentum" + self.fileSuffix)/np.sqrt(2)
				else:
					raise FileNotFoundError("unkown data source")
			except FileNotFoundError as e:
				print("failed to load phononCoordP:" + str(e))
				self.phononCoordP = np.full((self._time.size, self.num_sites), np.nan)
		if site is not None:
			retVal =  self.phononCoordP[:, site]
		else:
			retVal =  self.phononCoordP[:]
			
		if apply_cut:
			return retVal[:self.cut]
		else:
			return retVal
	
	
	def q(self, apply_cut=True):
		if self.num_sites != 2:
			raise RuntimeError("q coordinate undefined because num_sites == %i != 2"%self.num_sites)
		if not hasattr(self, "phononCoordQ"):
			try:
				if self.dataSource == "ED":
					self.phononCoordQ = self._load_ed_file("q1")
				elif self.dataSource == "DMRG":
					DataManager.q(self)
				else:
					raise FileNotFoundError("unkown data source")
			except FileNotFoundError as e:
				print("failed to load phononCoordQ:" + str(e))
				self.phononCoordQ = np.full((self._time.size, self.num_sites), np.nan)
		
		if apply_cut:
			return self.phononCoordQ[:self.cut]
		else:
			return self.phononCoordQ
	
	
	def q2(self):
		if self.num_sites != 2:
			raise RuntimeError("q coordinate undefined because num_sites == %i != 2"%self.num_sites)
		if not hasattr(self, "phononCoordQ2"):
			self.phononCoordQ2 = self._load_ed_file("q2"%i)
		return self.phononCoordQ2[:self.cut]
	
	
	def q3(self):
		if self.num_sites != 2:
			raise RuntimeError("q coordinate undefined because num_sites == %i != 2"%self.num_sites)
		if not hasattr(self, "phononCoordQ3"):
			self.phononCoordQ3 = self._load_ed_file("q3")
		return self.phononCoordQ3[:self.cut]
	
	
	def q4(self):
		if self.num_sites != 2:
			raise RuntimeError("q coordinate undefined because num_sites == %i != 2"%self.num_sites)
		if not hasattr(self, "phononCoordQ4"):
			self.phononCoordQ4 = self._load_ed_file("q4")
		return self.phononCoordQ4[:self.cut]
	
	
	def pq(self, apply_cut=True):
		if self.num_sites != 2:
			raise RuntimeError("pq coordinate undefined because num_sites == %i != 2"%self.num_sites)
		if not hasattr(self, "phononCoordPQ"):
			try:
				if self.dataSource == "ED":
					self.phononCoordPQ = self._load_ed_file("qp1")
				elif self.dataSource == "DMRG":
					DataManager.pq(self)
				else:
					raise FileNotFoundError("unkown data source")
			except FileNotFoundError as e:
				print("failed to load phononCoordQ:" + str(e))
				self.phononCoordPQ = np.full((self._time.size, self.num_sites), np.nan)
		
		if apply_cut:
			return self.phononCoordPQ[:self.cut]
		else:
			return self.phononCoordPQ
	
	
	def pq2(self):
		if self.num_sites != 2:
			raise RuntimeError("q coordinate undefined because num_sites == %i != 2"%self.num_sites)
		if not hasattr(self, "phononCoordPQ2"):
			self.phononCoordPQ2 = self._load_ed_file("qp2")
		return self.phononCoordPQ2[:self.cut]
	
	
	def pq3(self):
		if self.num_sites != 2:
			raise RuntimeError("q coordinate undefined because num_sites == %i != 2"%self.num_sites)
		if not hasattr(self, "phononCoordPQ3"):
			self.phononCoordPQ3 = self._load_ed_file("qp3")
		return self.phononCoordPQ3[:self.cut]
	
	
	def pq4(self):
		if self.num_sites != 2:
			raise RuntimeError("q coordinate undefined because num_sites == %i != 2"%self.num_sites)
		if not hasattr(self, "phononCoordPQ4"):
			self.self._load_ed_file("qp4")
		return self.phononCoordPQ4[:self.cut]


	def _get_dmrg_num_phonons(self, use_new_files):
		num_phonon = 0
		if use_new_files:
			if self.omega == 0.1 and self.gamma == 0.4:
				if self.initial_state == "bare":
					num_phonon = 20
				elif self.initial_state == "dressed":
					if self.num_sites == 5:
						num_phonon = 50
					else:
						num_phonon = 40
			elif self.initial_state == "dressed" and self.omega == 0.1 and self.gamma == 0.1:
				num_phonon = 30
			# TODO incomplete, e.g. 2,4 not found
			elif self.initial_state == "dressed" and self.omega == 2.:
				if self.gamma == 1.:
					num_phonon = 20
				elif self.gamma == 3.:
					num_phonon = 30
				elif self.gamma == 4.:
					num_phonon = 40
		else:
			if self.omega == 0.1:
				if self.gamma == 0.05 or self.gamma == 0.1:
					num_phonon = 20
				elif self.gamma == 0.2:
					num_phonon = 30
				elif self.gamma == 0.28:
					num_phonon = 40
				elif self.gamma == 0.4:
					num_phonon = 60
			elif self.omega == 0.2:
				if self.gamma == 0.2:
					num_phonon = 20
				elif self.gamma == 0.45:
					num_phonon = 30
				elif self.gamma == 0.6:
					num_phonon = 40
			elif self.omega == 0.5:
				if self.gamma == 0.25:
					num_phonon = 20
				elif self.gamma == 0.5:
					num_phonon = 30
				elif self.gamma == 1.:
					num_phonon = 30
				elif self.gamma == 2.:
					num_phonon = 70
			elif self.omega == 1.:
				if self.gamma == 0.05:
					num_phonon = 10
				elif self.gamma == 0.5:
					num_phonon = 30
				elif self.gamma == 1.:
					num_phonon = 20
				elif self.gamma == 2.:
					num_phonon = 30
			elif self.omega == 2.:
				if self.gamma == 0.05:
					num_phonon = 20
				elif self.gamma == 1:
					num_phonon = 10
				elif self.gamma == 3:
					num_phonon = 30
				elif self.gamma == 4:
					num_phonon = 40
		return num_phonon
	"""
	# for future adaptations to automatically get the available files from the folder:
	ref_file_or_folder = [f for f in listdir(path_to_ref)
		                     if isfile(join(path_to_ref, f)) and f[len(prefix_dmrg13):len(prefix_dmrg13)+4] == "time"]
	
	prefix_dmrg13 = path_to_ref + prefix_dmrg13
	
	phonon_suffix = "_truncOM1e-07_TEBDOM.txt"
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
	"""


class Wigner(DataManager):
	"""class for storage and access to the wigner functions
	
	available fields:
	q_sampled  	sampled q-coordinates
	p_q_sampled	sampled p_q-coordinates
	wigner   	values of the wigner functions. indexed as wigner[p_q, q]
	time    	all available times to load. loading reuquires the index instead of the value
	qi, pqi   	i-th moment of q or pq for all times
	"""
	
	def __init__(self, omega, gamma, initial_state, psi_initial="local", load=0):
		super(Wigner, self).__init__()
		self.initial_state = initial_state
		self.psi_initial = psi_initial
		self.omega = omega
		self.gamma = gamma
			
		self.path = "wigner_functions/%s_%s_omega=%g_gamma=%g/"%("BCDW" if initial_state == "bare" else "DCDW", psi_initial, omega, gamma)

		l0 = omega/gamma#/2#/np.sqrt(2)
		self.wigner_scale_fact_q = 1/l0*np.sqrt(2)
		self.wigner_scale_fact_p_q = l0*np.sqrt(2)


		f = open(self.path+"Q.DAT", 'rb')
		f.seek(4)	# skip four bytes which tell the size of the array
		self.q_sampled = np.fromfile(f, dtype='float64')*self.wigner_scale_fact_q
		f.close()
		
		f = open(self.path+"P.DAT", 'rb')
		f.seek(4)	# skip four bytes which tell the size of the array
		self.p_q_sampled = np.fromfile(f, dtype='float64')*self.wigner_scale_fact_p_q
		f.close()

		self.dq = self.q_sampled[1]-self.q_sampled[0]
		self.dp = self.p_q_sampled[1]-self.p_q_sampled[0]
		
		f = open(self.path+"TIMES.DAT", 'rb')
		f.seek(4)	# skip four bytes which tell the size of the array
		self._time = np.fromfile(f, dtype='float64')
		f.close()
		
		self.loaded_time = -1
		self.load(load)
		
		self.load_moments()
		self.reset_cut_time()
	
	
	def load(self, time_index : int):
		"""loads the wigner functions corresponding to the given index and returns it
		
		if any int value outside of [0, self.time.size) is passed, a ValueError is raised
		"""
		# TODO if time_index is float, search in self.time?
		if time_index < 0 or time_index >= self._time.size:
			raise ValueError("the given time index %i is outside of the available times"%i)
		elif self.loaded_time != time_index:
			f = open(self.path+"%i.DAT"%(time_index), 'rb')
			f.seek(4)	# skip four bytes which tell the size of the array
			self.wigner = np.fromfile(f, dtype=np.float64).reshape(self.q_sampled.size, self.p_q_sampled.size)\
				          /abs(self.wigner_scale_fact_q*self.wigner_scale_fact_p_q)
			f.close()
			self.loaded_time = time_index
		
		return self.wigner
	
	
	def nel(self, apply_cut=True):
		if not hasattr(self, "electronDensity"):
			try:
				self.electronDensity = np.loadtxt(self.path+"PROBABILITIES.DAT")[:, [3,4]]
			except OSError:
				try:
					ref = RefData(self.initial_state, self.psi_initial, 2, self.omega, self.gamma)
					if not np.allclose(ref.time()[1], self.time[1]):
						raise RuntimeError("wigner took nel from RefData, but timestep differs")
					self.electronDensity = ref.nel()[:self.time.size]
				except Exception:
					self.electronDensity = None
		
		if self.electronDensity is None:
			raise RuntimeError("No electron density information found")
				
		if apply_cut:
			return self.electronDensity[:self.cut]
		else:
			return self.electronDensity
	
	
	def nel_BO(self, apply_cut=True):
		if not hasattr(self, "electronDensity_BO"):
			try:
				self.electronDensity_BO = np.loadtxt(self.path+"PROBABILITIES.DAT")[:, [1,2]]
			except OSError:
				self.electronDensity_BO = None
		
		if self.electronDensity_BO is None:
			raise RuntimeError("No electron density information found")
		
		if apply_cut:
			return self.electronDensity_BO[:self.cut]
		else:
			return self.electronDensity_BO
	
	def q_density(self):
		return self.wigner.sum(axis=0)*(self.dp*self.dq)
	
	
	def pq_density(self):
		return self.wigner.sum(axis=1)*(self.dp*self.dq)
	
	
	def q_moment(self, order):
		return (self.q_density()*self.q_sampled**order).sum()
	
	
	def pq_moment(self, order):
		return (self.pq_density()*self.p_q_sampled**order).sum()
	
	def load_moments(self):
		try:
			self.q1 = np.load(self.path+"q1.npy")
			self.q2 = np.load(self.path+"q2.npy")
			self.q3 = np.load(self.path+"q3.npy")
			self.q4 = np.load(self.path+"q4.npy")
			
			self.pq1 = np.load(self.path+"pq1.npy")
			self.pq2 = np.load(self.path+"pq2.npy")
			self.pq3 = np.load(self.path+"pq3.npy")
			self.pq4 = np.load(self.path+"pq4.npy")
		except FileNotFoundError:
			self.recalc_moments()
			
	
	def recalc_moments(self):
		old_loaded_time = self.loaded_time
		print("recalculating moments in %s..."%self.path, end="", flush=True)
		
		self.q1 = np.empty(self._time.size)
		self.q2 = np.empty(self._time.size)
		self.q3 = np.empty(self._time.size)
		self.q4 = np.empty(self._time.size)
		
		self.pq1 = np.empty(self._time.size)
		self.pq2 = np.empty(self._time.size)
		self.pq3 = np.empty(self._time.size)
		self.pq4 = np.empty(self._time.size)
		
		for i in range(self._time.size):
			self.load(i)
			self.q1[i] = self.q_moment(1)
			self.q2[i] = self.q_moment(2)
			self.q3[i] = self.q_moment(3)
			self.q4[i] = self.q_moment(4)
			
			self.pq1[i] = self.pq_moment(1)
			self.pq2[i] = self.pq_moment(2)
			self.pq3[i] = self.pq_moment(3)
			self.pq4[i] = self.pq_moment(4)
		
		np.save(self.path+"q1.npy", self.q1)
		np.save(self.path+"q2.npy", self.q2)
		np.save(self.path+"q3.npy", self.q3)
		np.save(self.path+"q4.npy", self.q4)
		
		np.save(self.path+"pq1.npy", self.pq1)
		np.save(self.path+"pq2.npy", self.pq2)
		np.save(self.path+"pq3.npy", self.pq3)
		np.save(self.path+"pq4.npy", self.pq4)
		
		# restore old loaded wigner_function
		self.load(old_loaded_time)
		print("done.")




# TODO new shared superclass for both Simualation classes with load/store functionality
class Simulation(DataManager):
	""" Contianer class for Simulations with MTE, MCEv1 and MCEv2. superclass of Simulation_MTE and Simulation_MCE"""
	
	def __init__(self, hopping, omega, gamma,
	             num_sites_x, num_sites_y,
	             time_max, time_step,
	             initial_state, psi_initial, psi_state=None,
	             phonon_dressed_displaced_sites=None,
	             fullStorage=False,
	             force_new_simulation=False,
	             load_sampled_state=False,
	             override_export_path=None, extra_subfolder=None, # path options for previous 2 options
	             calc_BO_electron_density=False,
	             BO_nel_save_dist=1):
		super(Simulation, self).__init__()
		
		# sanity checks
		if psi_initial not in ["local", "eigen", "CDW"]:
			raise ValueError("psi init must be one of \"local\", \"eigen\", \"CDW\", where \"CDW\" is not supported by MCE")
		
		
		if int(num_sites_x) != num_sites_x or num_sites_x < 1 or int(num_sites_y) != num_sites_y or num_sites_y < 1:
			raise ValueError("num_sites_x/y must be a positive integer")
		
		# store parameters
		self.hopping = hopping
		self.omega = omega
		self.gamma = gamma
		
		self.num_sites_x = int(num_sites_x)
		self.num_sites_y = int(num_sites_y)
		self.num_sites   = self.num_sites_x*self.num_sites_y
		
		if psi_initial == "CDW" and self.num_sites == 2:
			psi_initial = "local"
		if psi_initial == "CDW":
			if psi_state is not None:
				print("WARNING: psi_state is not None, but is ignored for psi_initial == \"CDW\"")
			
			# every second site, if even, biased to (lower) left, which has smallest number of electrons
			psi_state = []
			for y in range(0, self.num_sites_y, 1):
				for x in range((self.num_sites_x + y)%2, self.num_sites_x, 2):
					psi_state.append((y, x))
		
		elif psi_initial == "local":
			def is_valid_pos(pos):
				return type(psi_state) == tuple and len(psi_state) == 2 and 0 <= pos[0] <= self.num_sites_y and 0 <= pos[1] <= self.num_sites_x
			
			if psi_state is None:
				psi_state = [((self.num_sites_y-1)//2, (self.num_sites_x-1)//2)]
			elif is_valid_pos(psi_state):
				psi_state = [psi_state]
			elif type(psi_state) != list or False in [is_valid_pos(pos) for pos in psi_state]:
				raise ValueError("for psi_initial == \"local\", psi_state must be either a tuple of length 2 or a list of tuples of length 2, "\
				                 "each with values (y,x) in ([0, num_sites_y), [0, num_sites_x)), or None")
		
		else: # "eigen"
			def is_valid_state(state):
				return type(psi_state) == int and 0 <= psi_state < self.num_sites
			
			if psi_state is None:
				psi_state = [0]
			elif is_valid_state(psi_state):
				psi_state = [psi_state]
			elif type(psi_state) != list or False in [is_valid_state(state) for state in psi_state]:
				raise ValueError("for psi_initial == \"eigen\", psi_state must be either an int or a list of int with values in [0, num_sites) or None")
		
		self.num_electrons = len(psi_state)
		
		if initial_state not in ["bare", "dressed"]:
			raise ValueError("initial state must be \"bare\" or \"dressed\"")
		
		if initial_state == "dressed":
			if psi_initial == "CDW":
				if phonon_dressed_displaced_sites is not None:
					import warning
					warning.warn("phonon_dressed_displaced_sites is not None, but is ignored for psi_initial == \"CDW\"")
				phonon_dressed_displaced_sites = psi_state
			else:
				if phonon_dressed_displaced_sites is None:
					phonon_dressed_displaced_sites = []
					if psi_initial == "local":
						phonon_dressed_displaced_sites.extend(psi_state)
					
				elif is_valid_pos(phonon_dressed_displaced_sites):
					phonon_dressed_displaced_sites = [phonon_dressed_displaced_sites]
				elif type(phonon_dressed_displaced_sites) != list or False not in [is_valid_pos(pos) for pos in phonon_dressed_displaced_sites]:
					raise ValueError("phonon_dressed_displaced_sites must be either a tuple of length 2 or a list of tuples of length 2, "\
					                 "each with values (y,x) in ([0, num_sites_y), [0, num_sites_x)), or None")
			
			if len(phonon_dressed_displaced_sites) == 0:
				raise ValueError(f"dressed sites couldn't be determined with psi_initial={psi_initial}. "\
				                 "Use phonon_dressed_displaced_sites to specify at least one site or use initial_state=\"bare\" instead")
		else:
			phonon_dressed_displaced_sites = None
		self.phonon_dressed_displaced_sites = phonon_dressed_displaced_sites
		
		
		self.initial_state = initial_state
		self.psi_initial   = psi_initial
		self.psi_state     = psi_state
		self.BO_nel_save_dist = BO_nel_save_dist
		
		# try to load old simulation data:
		if force_new_simulation or fullStorage:
			self.run_simulation = True
		else:
			succeeded_load = self.load_observables(override_export_path=override_export_path,
			                                       extra_subfolder=extra_subfolder, silent=True)
			self.run_simulation = not succeeded_load or self._time[-1] < time_max\
			                      or calc_BO_electron_density and self.obs_nel_BO is None
		
		if not self.run_simulation or load_sampled_state:
			try:
				self.load_initial_binary(override_export_path=override_export_path,
				                         extra_subfolder=extra_subfolder, silent=True)
				if self.run_simulation:
					print("reusing old sampled initial distribution")
			except Exception as e:
				print("failed to load_initial_binary:", e)
		
		
		if self.run_simulation:
			self.num_times = int(np.round(time_max/time_step)) + 1  # +1 for t=0
			self._time = np.linspace(0, time_max, self.num_times)
		
		self.set_cut_time(time_max)
	
	def set_export_path(self, new_export_path):
		if type(new_export_path) is not str:
			raise TypeError("new_export_path muste be a str, but was "\
			                 + str(type(new_export_path)) + " instead")
		self.export_path
	
	
	def export_path(self):
		return self._export_path
	
	
	def get_initial_state_name(self):
		if self.psi_initial != "CDW" and self.psi_state is None:
			raise RuntimeError("psi_state must be specified to generate the initial_state name if psi_initial != \"CDW\"")
		
		text_init = self.initial_state
		if self.initial_state == "dressed" and self.psi_initial != "CDW"\
		   and self.psi_state != self.phonon_dressed_displaced_sites:
			if len(self.phonon_dressed_displaced_sites) == 1:
				text_init += str(self.phonon_dressed_displaced_sites[0])
			else:
				text_init += str(self.phonon_dressed_displaced_sites)
		
		text_init += "_" + self.psi_initial
		if self.psi_initial != "CDW":
			if len(self.psi_state) == 1:
				simplified_state = self.psi_state[0]
			else:
				simplified_state = self.psi_state
			
			text_init += str(simplified_state).replace(' ', '')
		
		return text_init
	
	
	def get_sim_name(self):
		raise NotImplementedError()
	
	
	def get_storage_path(self, override_export_path=None, extra_subfolder=None):
		# setup the path to the storage location
		if override_export_path is not None:
			if type(override_export_path) is not str:
				raise TypeError("new_export_path muste be a str, but was "\
			                 + str(type(override_export_path)) + " instead")
			path = override_export_path
		else:
			path = self.export_path()
		
		if len(path) != 0 and path[-1] != "/":
			path += "/"
		
		path += self.get_sim_name() + "/"
		
		if extra_subfolder is not None:
			if type(extra_subfolder) is not str:
				raise TypeError("extra_subfolder muste be a str, but was "\
			                 + str(type(extra_subfolder)) + " instead")
			path += extra_subfolder
		
		if path[-1] != "/":
			path += "/"
		
		return path.replace(" ", "")
	
	
	def save_observables(self, override_export_path=None, extra_subfolder=None,
	                     delimiter="\t", save_fmt="%-16.14g", silent=False):
		path = self.get_storage_path(override_export_path, extra_subfolder)
		if not silent:
			print("saving observables to:", path)
		# create the folder (if not yet existing)
		os.system("mkdir -p " + path.replace("(", "\(").replace(")", "\)"))
		
		# save the observables as (time, obs) in a txt file
		shaped_time = self.time(apply_cut=False).reshape(-1, 1)
		header_start = "time  "
		
		site_list_2d = ""
		for y in range(self.num_sites_y):
			for x in range(self.num_sites_x):
				site_list_2d += delimiter+"(%2i, %2i)"%(x, y)
		
		if hasattr(self, "finalInitParams"):
			try:
				import json
				with open(path+"finalInitParams.json", "w") as f:
					json.dump(self.finalInitParams, f, indent=2)
			except Exception as e:
				print("saving finalInitParams failed:", e)
		
		try:
			header_obs = "electron density in real space on site (x, y)\n"\
			             + header_start + site_list_2d
			
			shaped_nel = self.nel(apply_cut=False).reshape(shaped_time.size, -1)
			np.savetxt(path + "electron_density.txt", np.hstack((shaped_time, shaped_nel)),
			           save_fmt, delimiter, header=header_obs)
		except Exception as e:
			print("saving local electron_density failed:", e)
		
		if self.obs_nel_BO is not None:
			try:
				header_obs = "electron density in Born-Oppenheimer Basis (BO), "\
				             "state number starting at 0 = groundstate, sorted by energy\n"\
				             + header_start + "\t"\
				             + str(np.arange(self.num_sites))[1:-1].replace(" ", "       "+delimiter)
				
				shaped_time_BO = self.time(apply_cut=False, BO_steps=True).reshape(-1, 1)				
				shaped_nel_BO = self.nel_BO(apply_cut=False).reshape(shaped_time_BO.size, -1)
				np.savetxt(path + "electron_density_BO.txt", np.hstack((shaped_time_BO, shaped_nel_BO)),
				           save_fmt, delimiter, header=header_obs)
			except Exception as e:
				print("saving electron_density_BO failed:", e)
		
		try:
			header_obs = "phonon number on site (x, y)\n"\
			             + header_start + site_list_2d
			
			shaped_nph = self.nph(apply_cut=False).reshape(shaped_time.size, -1)
			np.savetxt(path + "phonon_number.txt", np.hstack((shaped_time, shaped_nph)),
			           save_fmt, delimiter, header=header_obs)
		except Exception as e:
			print("saving local phonon_number failed:", e)
		
		try:
			header_obs = "lattice displacement of oscillator on site (x, y)\n"\
			             + header_start + site_list_2d
			
			shaped_x = self.x(apply_cut=False).reshape(shaped_time.size, -1)
			np.savetxt(path + "x.txt", np.hstack((shaped_time, shaped_x)),
			           save_fmt, delimiter, header=header_obs)
		except Exception as e:
			print("saving lattice displacement failed:", e)
		
		try:
			header_obs = "lattice momentum of oscillator on site (x, y)\n"\
			             + header_start + site_list_2d
			
			shaped_p = self.p(apply_cut=False).reshape(shaped_time.size, -1)
			np.savetxt(path + "p.txt", np.hstack((shaped_time, shaped_p)),
			           save_fmt, delimiter, header=header_obs)
		except Exception as e:
			print("saving lattice momentum failed:", e)
		
		try:
			header_obs = "kinetic energy of electron(s)\n"\
			             + header_start + "\t"\
			             + str(np.arange(1, self.num_electrons+1))[1:-1].replace(" ", "       "+delimiter)
			
			shaped_el_ekin = self.energy_kinetic("all", True, apply_cut=False).reshape(shaped_time.size, -1)
			np.savetxt(path + "energy_kinetic.txt", np.hstack((shaped_time, shaped_el_ekin)),
			           save_fmt, delimiter, header=header_obs)
		except Exception as e:
			print("saving lattice displacement failed:", e)
		
		try:
			header_obs = "interaction energy between electron(s) (summed over electrons) and oscillator on site (x, y)\n"\
			             + header_start + site_list_2d
			
			shaped_energy_interaction = self.energy_interaction("all", True, apply_cut=False).reshape(shaped_time.size, -1)
			np.savetxt(path + "energy_interaction.txt", np.hstack((shaped_time, shaped_energy_interaction)),
			           save_fmt, delimiter, header=header_obs)
		except Exception as e:
			print("saving interaction energy failed:", e)
		
		if type(self) == Simulation_MCE:
			try:
				header_obs = header_start + delimiter + "total energy"
				
				shaped_energy_total = self.energy_total(apply_cut=False).reshape(shaped_time.size, -1)
				np.savetxt(path + "energy_total.txt", np.hstack((shaped_time, shaped_energy_total)),
					       save_fmt, delimiter, header=header_obs)
			except Exception as e:
				print("saving interaction energy failed:", e)
			
			if hasattr(self, "fullStorage_num_configs") and self.fullStorage_num_configs is not None:
				try:
					header_obs = header_start + delimiter + "number of configurations"
					
					shaped_num_configs = self.fullStorage_num_configs.reshape(shaped_time.size, -1)
					np.savetxt(path + "num_configs.txt", np.hstack((shaped_time, shaped_num_configs)),
							   save_fmt, delimiter, header=header_obs)
				except Exception as e:
					print("saving num_configs failed:", e)
	
	
	def load_observables(self, override_export_path=None, extra_subfolder=None, silent=False):
		path = self.get_storage_path(override_export_path, extra_subfolder)
		
		# remove all old contents
		self._time = None
		self.obs_nel = None
		self.obs_nel_BO = None
		self.obs_nph = None
		self.obs_y = None
		self.obs_el_ekin = None
		self.obs_el_ph_inter = None
		
		try:
			import json
			with open(path+"finalInitParams.json", "r") as f:
				self.finalInitParams = json.load(f)
		except Exception as e:
			if not silent:
				print("loading finalInitParams failed:", e)
		
		try:
			data = np.loadtxt(path + "electron_density.txt", ndmin=2)
			self._time = data[:, 0]
			self.obs_nel = data[:, 1:].reshape(self._time.size, self.num_sites_y, self.num_sites_x)
		except Exception as e:
			if not silent:
				print("failed loading electron_density:", e)
		
		try:
			data = np.loadtxt(path + "phonon_number.txt", ndmin=2)
			if self._time is None:
				self._time = data[:, 0]
			self.obs_nph = data[:, 1:].reshape(self._time.size, self.num_sites_y, self.num_sites_x)
		except Exception as e:
			if not silent:
				print("loading phonon_number failed:", e)
		
		try:
			data_x = np.loadtxt(path + "x.txt", ndmin=2)
			data_p = np.loadtxt(path + "p.txt", ndmin=2)
			if self._time is None:
				self._time = data_x[:, 0]
			
			x = data_x[:, 1:].reshape(self._time.size, self.num_sites_y, self.num_sites_x)
			p = data_p[:, 1:].reshape(self._time.size, self.num_sites_y, self.num_sites_x)
			self.obs_y = np.stack((x, p), axis=-1)
		except Exception as e:
			if not silent:
				print("loading lattice displacement and momentum failed:", e)
		
		try:
			data = np.loadtxt(path + "energy_kinetic.txt", ndmin=2)
			if self._time is None:
				self._time = data[:, 0]
			
			self.obs_el_ekin = data[:, 1:].reshape(self._time.size, -1)\
			                   *(1/self.hopping if self.hopping != 0 else 1)
		except Exception as e:
			if not silent:
				print("loading lattice energy_kinetic failed:", e)
		
		try:
			data = np.loadtxt(path + "energy_interaction.txt", ndmin=2)
			if self._time is None:
				self._time = data[:, 0]
			self.obs_el_ph_inter = data[:, 1:].reshape(self._time.size, self.num_sites_y, self.num_sites_x)\
			                       *(1/self.gamma if self.gamma != 0 else 1)
		except Exception as e:
			if not silent:
				print("loading lattice energy_interaction failed:", e)
		
		try:
			data = np.loadtxt(path + "electron_density_BO.txt", ndmin=2)
			if self._time is None:
				self._time = data[:, 0]
				self.BO_nel_save_dist = 1
			else:
				self.BO_nel_save_dist = int(np.round(data[1, 0]/self._time[1]))
			self.obs_nel_BO = data[:, 1:]
		except Exception as e:
			self.BO_nel_save_dist = 1
			if not silent:
				print("loading electron_density_BO failed:", e)
		
		if type(self) == Simulation_MCE and self.num_configs != self.max_num_configs:
			try:
				data = np.loadtxt(path + "num_configs.txt", ndmin=2)
				self.fullStorage_num_configs = data[:, 1]
			except Exception as e:
				print("loading num_configs failed:", e)
		
		
		succeeded = True
		if self.obs_y is None or self.obs_nph is None or self.obs_nel is None\
			   or self.obs_el_ekin is None or self.obs_el_ph_inter is None\
			   or self._time is None:
			succeeded = False
		
		if not silent:
			if succeeded:
				print("loaded data from", path)
			else:
				print("failed to load data from", path)
		return succeeded
	
	def save_initial_binary(self, override_export_path=None, extra_subfolder=None, silent=False):
		""" saves all sampled initial values as raw binary. Endianess of the system is never checked,
		so when transfering the binary output files between different machines, compare the endianess!
		All data are stored as np.complex/double complex (16 Byte) and to Stefan's knowledge,
		just little-endian machines have been used with this code.
		To check your endianess, run in python:
		from sys import byteorder
		print(byteorder) # prints "little" on Stefan's Laptop and on ds10 of rocks cluster
		
		in case you have different values between the two machines in the transfer,
		you'll need to convert the data by flipping the byteorder in every block of 16 bytes
		(code not included, but should be easy in C by using a char* reference)
		"""
		raise NotImplementedError()
		
		
	def load_initial_binary(self, override_export_path=None, extra_subfolder=None, silent=False):
		"""loads the data saved by save_initial_binary, READ THAT DOCUMENTATION for notice about data layout
		"""
		raise NotImplementedError()


from simulate import simulate_system_MTE
from src.ehrenfest_init import initial_state_MTE

class Simulation_MTE(Simulation):
	"""container class for MTE simulation data. The Simulation is run at object creation"""
	
	
	def get_sim_name(self):
		return "MTE_%s_Lx=%lu_Ly=%lu_hopping=%g_omega=%g_gamma=%g_numtraj=%i"\
		       %(self.get_initial_state_name(), self.num_sites_x, self.num_sites_y,\
		         self.hopping, self.omega, self.gamma,
		         self.num_samples)
	
	
	def save_initial_binary(self, override_export_path=None, extra_subfolder=None, silent=False):
		path = self.get_storage_path(override_export_path, extra_subfolder)
		if not silent:
			print("saving initial conditions to:", path)
		os.system("mkdir -p " + path.replace("(", "\(").replace(")", "\)"))
		
		with open(path+"MTE_raw_start.bin", "wb") as file_bin_raw_start:
			self.start_y  .tofile(file_bin_raw_start)
			self.start_psi.tofile(file_bin_raw_start)
	
	
	def load_initial_binary(self, override_export_path=None, extra_subfolder=None, silent=False):
		path = self.get_storage_path(override_export_path, extra_subfolder)
		
		with open(path+"MTE_raw_start.bin", "rb") as file_bin_raw_start:
			if not silent:
				print("loading initial conditions from:", path)
			item_count = self.num_samples*self.num_sites
			
			self.start_y   = np.fromfile(file_bin_raw_start, dtype=np.complex, count=item_count, offset=0)\
			                 .reshape((self.num_samples, self.num_sites_y, self.num_sites_x))
			self.start_psi = np.fromfile(file_bin_raw_start, dtype=np.complex, count=self.num_electrons*item_count)\
			                 .reshape((self.num_samples, self.num_electrons, self.num_sites_y, self.num_sites_x))
			
	
	
	#new signature:
	def __init__(self, hopping, omega, gamma,
	             num_sites_x, num_sites_y,
	             time_max, time_step,
	             num_samples,
	             initial_state, psi_initial, psi_state=None,
	             phonon_dressed_displaced_sites=None,
	             calcYSquared=True, fullStorage=False,
	             symmetric_y_init=False, tol_energy_error_per_site=1e-6,
	             printProgress=True,
	             force_new_simulation=False,
	             load_sampled_state=False,
	             store_results=True,
	             override_export_path=None, extra_subfolder="0"): # path options for previous 2 options
	    # TODO rewrite documentation
	    # TODO allow data save and load?
		"""run a new simulation of the 2d Holstein model. The created object contains the results.
		Some methods are only available if fullStorage=True.
		
		omega, gamma and hopping=1 are the system parameters.
		All classical coordinates (x, p, q, p_q) are stored in natural units.
		x -> x/l0, p->p*l0 where l0=sqrt((\\hbar=1)/mass_ho/omega) = 1
		
		time_max last time (inclusive) to which shall be integrated
		time_step the time step at which the data are stored
		tol_energy_error_per_site is used as upper limit of inetgration timestep as 
		time_step = min(0.05, sqrt(10*tol_energy_error_per_site*hopping/omega/gamma**2))
		otherwise the given time_step is used and tol_energy_error_per_site ignored.
		Every time step is stored (might change in future as option)
		
		num_samples number of classical trajectories which shall be generated (>=0).
		Multiplied by 4 if symmetric_y_init (default: True).
		if num_samples==0, a single trajectory at x_i=0, p_i=0 is used,
		otherwise, chosen from groundstate of quantum harmonic oscillator (gaussian).
		This parameter is ignored if an explicit initial state is given in y_initial.
		
		symmetric_y_init whether generated inital state is symmetrized.
		see num_samples
		
		num_sites number of sites in the system.
		may be omitted, if either y_initial or psi_initial are given.
		In that case, it is inferred from them.
		Must be consistent with shapes of y_initial and psi_initial (if they are given)
		
		y_initial  initial state of the classical particles to use.
		data layout: y_initial[sample, site_y, site_x, x or p]
		may be omitted, if num_samples and (num_sites or psi_initial) are given,
		in which case it is generated.
		
		num_electrons number of (non-interacting) electrons (>=0) in the system.
		may be ommited if psi_initial is given. No electrons are allowed.
		
		psi_initial initial wavefunction of the non-interacting electrons.
		currently only accepts strings "local" and "eigen"
		local: localized as CDW (charge density wave)
		eigen: eigenstate of Born-Oppenheimer Hamiltonian, to select one (or more) use psi_state
		
		initial_state must be either (if y_initial is not None)
		   "bare" (bare charge density wave)
		or "dressed" (dressed charge density wave).
		For "bare", the generated y_initial is centered around 0 (no phonons)
		For "dressed", the generated y_initial is centered around 0 at unoccupied sites
		and around gamma/omega*sqrt(2) at occupied sites
		if y_initial is not None, the argument is ignored
		
		fullStorage whether all classical trajectories shall be stored completly.
		some functions (like histograms in q space) are only available if fullStorage==True.
		defaults to False
		
		printProgress whether the progress of the simulation shall be printed.
		default: True
		
		"""
		if int(num_samples) != num_samples or num_samples <= 0:
			raise ValueError("num_samples must be a positive integer")
		self.num_samples = num_samples
		super(Simulation_MTE, self).__init__(hopping, omega, gamma,
		                                     num_sites_x, num_sites_y,
		                                     time_max, time_step,
		                                     initial_state, psi_initial, psi_state,
		                                     phonon_dressed_displaced_sites,
		                                     fullStorage,
		                                     force_new_simulation, load_sampled_state,
		                                     override_export_path, extra_subfolder,
		                                     calc_BO_electron_density=False, BO_nel_save_dist=1)
		
		
		if self.run_simulation:
			# TODO shift to ehrenfest_init
			x_shift = np.zeros((num_sites_y, num_sites_x))
			el_grid = np.zeros((num_sites_y, num_sites_x), dtype=np.bool)
			for p in self.psi_state:
				el_grid[p] = True
			
			if initial_state == "dressed":
				x_shift[el_grid] = gamma/omega*np.sqrt(2)
			
			y_initial = initial_state_MTE(self.num_sites_x, self.num_sites_y, self.num_samples, symmetric_y_init, self.omega, x_shift)
			# num_samples might have changed by symmetrization
			self.num_samples = y_initial.shape[0]
			
			if type(psi_initial) is str:
				start_psi = np.zeros((self.num_samples, self.num_electrons, self.num_sites_y, self.num_sites_x), dtype=np.complex)
				if psi_initial in ["local", "CDW"]:
					X, Y = np.meshgrid(np.arange(num_sites_x), np.arange(num_sites_y))
					x = X[el_grid].flatten()
					y = Y[el_grid].flatten()
					
					for i in range(len(y)):
						start_psi[0, i, y[i], x[i]] = 1
					start_psi[:] = start_psi[0].reshape(1, self.num_electrons, self.num_sites_y, self.num_sites_x)
				else:
					# FIXME
					# TODO psi_state for selection of eigenstate. might be an array
					# self.num_sites
					if psi_initial != "eigen":
						raise RuntimeError("unkown value of psi_initial encountered")
					if num_sites_y != 1:
						raise Exception("eigenstates only supported in 1d (num_sites_y == 1)")
					
					if psi_state is None:
						raise ValueError("psi_state is None, although an eigenstate is requested")
					
					# TODO also allow 2d version, which is not tridiagonal
					# construct Hamiltonian
					coupling_strength = self.gamma*np.sqrt(2)
					H = np.zeros((self.num_samples, self.num_sites, self.num_sites))
					diag_ind = np.arange(self.num_sites)
					
					H[:, diag_ind[:-1], diag_ind[1:]] = -self.hopping
					H[:, diag_ind[1:], diag_ind[:-1]] = -self.hopping
					H[:, diag_ind, diag_ind] = (-coupling_strength)*y_initial[:, :, :, 0].reshape(self.num_samples, self.num_sites)
					
					# find eigenvalues and eigenstates
					eigval, eigvecs = np.linalg.eigh(H)
					eigval = eigval.reshape(self.num_samples, self.num_sites)
					order = np.argsort(eigval, axis=1)
					
					# use eigenstates according to psi_state as start_psi
					for n in range(self.num_samples):
						start_psi[n] = eigvecs[n, :, order[n, psi_state]]
			# TODO explicit setup is missing?
			else:
				raise Exception("still need to code explicit starting...")
			
			
			self.start_psi = start_psi.copy()
			self.start_y   = y_initial.copy()
			res = simulate_system_MTE(y_initial, start_psi,
			                          self.hopping, self.omega, self.gamma,
			                          self._time.size, time_step, tol_energy_error_per_site,
			                          returnNumPhonons=False,
			                          returnYSquared=calcYSquared,
			                          printProgress=printProgress,
			                          fullSampleReturn=fullStorage,
			                          returnElEkin=True,
			                          returnInteractionEnergy=True)
			
			self.obs_y = res["y_mean"]
			self.obs_nel = res["density"]
			self.y_squared = res["y_squared"]
			self.y_full = res["y_full"]
			self.obs_el_ekin = res["el_kinetic_energy"]
			self.obs_el_ph_inter = res["energy_interaction"]
			
			self.obs_nel_BO = None
		
			try:
				self.obs_nph = ((  self.x2(squeezeY=False, apply_cut=False) \
				                 + self.p2(squeezeY=False, apply_cut=False))\
				                - (1. if self.num_samples > 1 else 0.))/2.
			except RuntimeError:
				self.obs_nph = np.full((self.num_times, self.num_sites_y, self.num_sites_x), np.nan)
			
			if store_results:
				self.save_observables(override_export_path=override_export_path,
				                      extra_subfolder=extra_subfolder, silent=not printProgress)
				self.save_initial_binary(override_export_path=override_export_path,
				                         extra_subfolder=extra_subfolder, silent=not printProgress)
		
		self.set_cut_time(time_max)
	
	
	def _calc_y_squared_from_y_full(self):
		if self.y_full is None:
			raise RuntimeError("calculating y_squared from y_full=None is impossible")
		self.y_squared = (self.y_full**2).mean(axis=1)
	
	
	def x2(self, site_x=None, site_y=None, squeezeY=True, apply_cut=True):
		if self.y_squared is None and self.y_full is not None:
			self._calc_y_squared_from_y_full()
		if self.y_squared is None:
			raise RuntimeError("no mean of squared x or p data available")
		
		return self._site_selection(self.y_squared[:, :, :, 0], site_x, site_y, squeezeY, apply_cut)
	
	
	def p2(self, site_x=None, site_y=None, squeezeY=True, apply_cut=True):
		if self.y_squared is None and self.y_full is not None:
			self._calc_y_squared_from_y_full()
		if self.y_squared is None:
			raise RuntimeError("no mean of squared x or p data available")
		
		return self._site_selection(self.y_squared[:, :, :, 1], site_x, site_y, squeezeY, apply_cut)
	
	
	def q2(self):
		if self.num_sites_x != 2 or self.num_sites_y != 1:
			raise RuntimeError("q coordinate undefined because num_sites == %i != 2"%self.num_sites)
		return self.x2(0)-self.x2(1)
	
	
	def q3(self):
		if self.num_sites_x != 2 or self.num_sites_y != 1:
			raise RuntimeError("q coordinate undefined because num_sites == %i != 2"%self.num_sites)
		if not hasattr(self, "phononCoordQ3"):
			self.phononCoordQ3 = (self.get_q_full()**3).mean(axis=1)
		return self.phononCoordQ3[:self.cut]
	
	
	def q4(self):
		if self.num_sites_x != 2 or self.num_sites_y != 1:
			raise RuntimeError("q coordinate undefined because num_sites == %i != 2"%self.num_sites)
		if not hasattr(self, "phononCoordQ4"):
			self.phononCoordQ4 = (self.get_q_full()**4).mean(axis=1)
		return self.phononCoordQ4[:self.cut]
	
	
	def pq2(self):
		if self.num_sites_x != 2 or self.num_sites_y != 1:
			raise RuntimeError("q coordinate is only defined for 2 sites")

		return self.p2(0)-self.p2(1)
	
	
	def pq3(self):
		if self.num_sites_x != 2 or self.num_sites_y != 1:
			raise RuntimeError("q coordinate is only defined for 2 sites")
		if not hasattr(self, "phononCoordPQ3"):
			self.phononCoordPQ3 = (self.get_p_q_full()**3).mean(axis=1)
		return self.phononCoordPQ3[:self.cut]
	
	
	def pq4(self):
		if self.num_sites_x != 2 or self.num_sites_y != 1:
			raise RuntimeError("q coordinate is only defined for 2 sites")
		if not hasattr(self, "phononCoordPQ4"):
			self.phononCoordPQ4 = (self.get_p_q_full()**4).mean(axis=1)
		return self.phononCoordPQ4[:self.cut]
	
	
	def get_q_full(self):
		"""returns the full storage of q coordinates (without p_q)
		
		The data are only calculated once and then stored.
		Calling this function if self.num_sites != 2 causes a RuntimeError
		only available if all y_full has been stored
		
		result shape: (num_times, num_samples)
		"""
		
		if self.num_sites_x != 2 or self.num_sites_y != 1:
			raise RuntimeError("q coordinate is only defined for 2 sites")
		if self.y_full is None:
			raise RuntimeError("full data are only available if the simulation stored all trajectories")
		
		if not hasattr(self, "q_full"):
		#num_samples, num_times, num_sites_y, num-sites_x
			self.q_full = self.y_full[:, 0, :, 0] - self.y_full[:, 0, :, 1]
		return self.q_full
	
	
	def get_p_q_full(self):
		"""returns the full storage of p_q coordinates (without q)
		
		The data are only calculated once and then stored.
		Calling this function if self.num_sites != 2 causes a RuntimeError
		only available if all y_full has been stored
		
		result shape: (num_times, num_samples)
		"""
		
		if self.num_sites_x != 2 or self.num_sites_y != 1:
			raise RuntimeError("q coordinate is only defined for 2 sites")
		if self.y_full is None:
			raise RuntimeError("full data are only available if the simulation stored all trajectories")
		
		if not hasattr(self, "p_q_full"):
			self.p_q_full = self.y_full[:, 0, :, 1] - self.y_full[:, 0, :, 1]
		return self.p_q_full
	
	
	def _update_statistics_q(self):
		if self.num_sites_x != 2 or self.num_sites_y != 1:
			raise RuntimeError("q coordinate is only defined for 2 sites")
		self.q_variance, self.q_skew, self.q_kurtosis = _get_statistic_estimators_from_sample(self.get_q_full(), self.q(truncate=False))
	
	
	def get_q_variance(self):
		"""returns the estimated variance of the q coordinate
		
		The data are only calculated once and then stored.
		Calling this function if self.num_sites != 2 causes a RuntimeError
		
		result shape: (num_times,)
		"""
		
		if not hasattr(self, "q_variance"):
			self._update_statistics_q()
		return self.q_variance[:self.cut]
	
	
	def get_q_skew(self):
		"""returns the estimated skew of the q coordinate
		
		The data are only calculated once and then stored.
		Calling this function if self.num_sites != 2 causes a RuntimeError
		
		result shape: (num_times,)
		"""
		
		if not hasattr(self, "q_skew"):
			self._update_statistics_q()
		return self.q_skew[:self.cut]
	
	
	def get_q_kurtosis(self):
		"""returns the estimated excess kurtosis of the q coordinate
		
		The data are only calculated once and then stored.
		Calling this function if self.num_sites != 2 causes a RuntimeError
		
		result shape: (num_times,)
		"""
		
		if not hasattr(self, "q_kurtosis"):
			self._update_statistics_q()
		return self.q_kurtosis[:self.cut]
	
	
	def _update_statistics_pq(self):
		if self.num_sites_x != 2 or self.num_sites_y != 1:
			raise RuntimeError("q coordinate is only defined for 2 sites")
		self.pq_variance, self.pq_skew, self.pq_kurtosis = _get_statistic_estimators_from_sample(self.get_p_q_full(), self.pq())
	
	
	def get_p_q_variance(self):
		"""returns the estimated variance of the p_q coordinate
		
		The data are only calculated once and then stored.
		Calling this function if self.num_sites != 2 causes a RuntimeError
		
		result shape: (num_times,)
		"""
		
		if not hasattr(self, "pq_variance"):
			self._update_statistics_pq()
		return self.pq_variance[:self.cut]
	
	def get_p_q_skew(self):
		"""returns the estimated skew of the p_q coordinate
		
		The data are only calculated once and then stored.
		Calling this function if self.num_sites != 2 causes a RuntimeError
		
		result shape: (num_times,)
		"""
		
		if not hasattr(self, "pq_skew"):
			self._update_statistics_pq()
		return self.pq_skew[:self.cut]
	
	def get_p_q_kurtosis(self):
		"""returns the estimated excess kurtosis of the p_q coordinate
		
		The data are only calculated once and then stored.
		Calling this function if self.num_sites != 2 causes a RuntimeError
		result shape: (num_times,)
		"""
		if not hasattr(self, "pq_kurtosis"):
			self._update_statistics_pq()
		return self.pq_kurtosis[:self.cut]
	
	# FIXME histograms don't have the new data layout yet
	def get_histogram_xp(self, histogram_bins=0):
		if self.y_full is None:
			raise RuntimeError("histograms are only available if the simulation stored all trajectories")
		
		if type(histogram_bins) == int and histogram_bins >= 0:
			if histogram_bins == 0:
				histogram_bins = int(5*self.num_samples**.25)
		else:
			raise ValueError("histogram_bins must be an non-negative integer")
		"""
			num_bins = (histogram_bins, histogram_bins)
			
		elif type(histogram_bins) == tuple and len(histogram_bins) == 2:
			num_bins = (0,0)
			if type(histogram_bins[0]) == int:
				num_bins[0] = histogram_bins[0]
				auto_bin_range_x = True
			else:
				num_bins[0] = histogram_bins[0][0]
				auto_bin_range_x = False
			
			if type(histogram_bins[1]) == int:
				num_bins[1] = histogram_bins[1]
				auto_bin_range_p = True
			else:
				num_bins[1] = histogram_bins[1][0]
				auto_bin_range_p = False
		
		if hasattr(self, histogram_xp_edges) and hasattr(self, histogram_xp):
			# check whether the old histogram can be returned
			
		
		"""
		# TODO store histogram? check for old one before recalc?
		# storage currently needed for hist_phase_space.py, but could be changed to function
		# define the edges of the histogram
		self.histogram_xp_edges = np.empty((2, histogram_bins+1))
		self.histogram_xp_edges[0] = np.linspace(self.y_full[:, :, 0::2].min(), self.y_full[:, :, 0::2].max(), histogram_bins+1)
		self.histogram_xp_edges[1] = np.linspace(self.y_full[:, :, 1::2].min(), self.y_full[:, :, 1::2].max(), histogram_bins+1)
		
		# calculate all histograms
		self.histogram_xp = np.empty((self.cut, self.num_sites, histogram_bins, histogram_bins))
		for t in range(self.cut):
			for s in range(self.num_sites):
				self.histogram_xp[t, s] = np.histogram2d(self.y_full[t, :, 2*s], self.y_full[t, :, 2*s+1],
				                                         bins=self.histogram_xp_edges,
				                                         density=True)[0]
		
		return self.histogram_xp, self.histogram_xp_edges
	
	
	def get_histogram_qpq(self, histogram_bins=0):
		if self.y_full is None:
			raise RuntimeError("histograms are only available if the simulation stored all trajectories")
		
		if type(histogram_bins) == int and histogram_bins >= 0:
			if histogram_bins == 0:
				histogram_bins = int(5*self.num_samples**.25)
		else:
			raise ValueError("histogram_bins must be an non-negative integer")
		
		if not hasattr(self, "q_full"):
			self.q_full   = self.y_full[:, :, 0] - self.y_full[:, :, 2]
			self.p_q_full = self.y_full[:, :, 1] - self.y_full[:, :, 3]
		
		
		# define the edges of the histogram
		self.histogram_q_edges = np.empty((2, histogram_bins+1))
		self.histogram_q_edges[0] = np.linspace(self.get_q_full().min(),   self.get_q_full().max(),   histogram_bins+1)
		self.histogram_q_edges[1] = np.linspace(self.get_p_q_full().min(), self.get_p_q_full().max(), histogram_bins+1)
		
		self.q_histogram = np.empty((self.cut, histogram_bins, histogram_bins))
		for t in range(self.cut):
			self.q_histogram[t] = np.histogram2d(self.get_q_full()[t, :], self.get_p_q_full()[t, :],
			                                     bins=self.histogram_q_edges,
			                                     density=True)[0]
		return self.q_histogram, self.histogram_q_edges
	
	
	def get_histogram_pqMinusq(self, histogram_bins=0):
		if self.y_full is None:
			raise RuntimeError("histograms are only available if the simulation stored all trajectories")
		
		self.data = self.get_p_q_full() - self.get_q_full()
		self.histogram_edges = np.linspace(data.min(), data.max(), histogram_bins+1)
		
		self.hist = np.empty((self._time.size, histogram_bins))
		for t in range(self._time.size):
			self.hist[t], _ = np.histogram(data[t, :], bins=histogram_bins, density=True)
		
		return self.hist, self.histogram_edges
	
def _get_statistic_estimators_from_sample(data_full, data_mean):
	""" calculates the unbiased estimators of variance, skew and excess kurtosis"""
	
	n = data_full.shape[1]
	full_centered = data_full - data_mean.reshape(-1, 1)
	full_centered_square = full_centered**2
	variance = full_centered_square.sum(axis=1)/(n-1)
	skew = (full_centered_square*full_centered).sum(axis=1)/variance**1.5*(n/(n-1)/(n-2))
	kurtosis = (full_centered_square**2).sum(axis=1)/variance**2*(n*(n+1)/(n-1)/(n-2)/(n-3))
	
	return variance, skew, kurtosis - 3*(n-1)**2/(n-2)/(n-3)



from simulate import simulate_system_MCE
from src.ehrenfest_init import initial_state_MCE

class Simulation_MCE(Simulation):
	"""container class for Multi-Configurational Ehrenfest (MCE) simulation data.
	   Both Versions (1 and 2) are supported.
	   at creation, MCE_binary_results is checked for simulation data, if none are found, a new simulation is run"""
	
	
	def get_sim_name(self):
		name = "MCEv%i_%s_Lx=%lu_Ly=%lu_hopping=%g_omega=%g_gamma=%g"\
		       %(self.MCE_version, self.get_initial_state_name(), self.num_sites_x, self.num_sites_y,\
		         self.hopping, self.omega, self.gamma)
		
		if self.local_psi_init_select_eigenstate:
			name += "_initEigen"
		if self.flip_initial_state_for_half_configs:
			name += "_initFlip"
		if self.scale_el == 0:
			name += "_scale_el=0"
		
		if self.pancakeSampling:
			name += "_pancake%g"%self.pancakeScale
		
		
		name += "_numconf=%i-%i"%(self.num_configs, self.max_num_configs)
		
		return name
	
	
	# may be used during creation as long as get_sim_name() works and self.num_sites_x,y is set
	def load_binary_observable(self, name):
		input_path = "MCE_binary_results/"+ self.get_sim_name() + "/"
		output_path = input_path
		
		target_shape = None
		if name == "obs_y":
			target_shape = (-1, self.num_sites_y, self.num_sites_x, 2)
		elif name == "obs_nel":
			target_shape = (-1, self.num_sites_y, self.num_sites_x)
		elif name == "obs_nph":
			target_shape = (-1, self.num_sites_y, self.num_sites_x)
		elif name == "obs_el_ekin":
			target_shape = (-1, 1)
		elif name == "obs_el_ph_inter":
			target_shape = (-1, self.num_sites_y, self.num_sites_x)
		elif name != "time":
			raise RuntimeError("unknown observable \""+name+"\" requested")
		
		# TODO also check for already exsisting npy file?
		f = open(input_path + name + ".bin", "rb")
		arr = np.fromfile(f, dtype=np.float64)
		f.close()
		
		if target_shape != None:
			arr = arr.reshape(*target_shape)
		
		return arr
	
	def save_initial_binary(self, override_export_path=None, extra_subfolder=None, silent=False):
		path = self.get_storage_path(override_export_path, extra_subfolder)
		if not silent:
			print("saving initial conditions to:", path)
		os.system("mkdir -p " + path.replace("(", "\(").replace(")", "\)"))
		
		with open(path+"MCE%i_raw_start.bin"%self.MCE_version, "wb") as file_bin_raw_start:
			self.start_hoz.tofile(file_bin_raw_start)
			self.start_el_coeff.tofile(file_bin_raw_start)
			if self.MCE_version == 2:
				self.start_coupling_coeff.tofile(file_bin_raw_start)
	
	
	def load_initial_binary(self, override_export_path=None, extra_subfolder=None, silent=False):
		path = self.get_storage_path(override_export_path, extra_subfolder)
		
		with open(path+"MCE%i_raw_start.bin"%self.MCE_version, "rb") as file_bin_raw_start:
			if not silent:
				print("loading initial conditions from:", path)
			item_count = self.num_configs*self.num_sites
			
			self.start_hoz = np.fromfile(file_bin_raw_start, dtype=np.complex, count=item_count, offset=0)
			self.start_el_coeff = np.fromfile(file_bin_raw_start, dtype=np.complex, count=item_count)
			
			if self.MCE_version == 1:
				self.start_hoz = self.start_hoz.reshape(self.num_sites_y, self.num_sites_x, self.num_configs)
				self.start_el_coeff = self.start_el_coeff.reshape(self.num_sites_y, self.num_sites_x, self.max_num_configs)
				self.start_coupling_coeff = None
				
			elif self.MCE_version == 2:
				self.start_hoz = self.start_hoz.reshape(self.num_configs, self.num_sites_y, self.num_sites_x)
				self.start_el_coeff = self.start_el_coeff.reshape(self.num_configs, self.num_sites_y, self.num_sites_x)
				self.start_coupling_coeff = np.fromfile(file_bin_raw_start, dtype=np.complex)
				
			else:
				raise ValueError("invalid MCE_version encountered in load_initial_binary")
	
	
	def __init__(self, version, hopping, omega, gamma,
	             num_sites_x, num_sites_y,
	             time_max, time_step,
	             num_configs,
	             initial_state, psi_initial, psi_state=None,
	             max_num_configs=None, allow_config_delete=False,
	             allow_config_clone=True,
	             fullStorage=False, printProgress=True,
	             phonon_dressed_displaced_sites=None,
	             force_new_simulation=False,
	             load_sampled_state=False,
	             store_results=True,
	             override_export_path=None, extra_subfolder="0", # path options for previous 2 options
	             calc_BO_electron_density=False,
	             BO_nel_save_dist=1,
	             local_psi_init_select_eigenstate=None,
	             seed=None, allowCompression=True, allowExpansion=None,
	             scale_ph=None, scale_el=None,
	             pancakeSampling=False, pancakeScale=1.,
	             flip_initial_state_for_half_configs=False,
	             calcCondition=False):
		if version != 1 and version != 2:
			raise ValueError("Version of MCE simulation must be either 1 or 2")
		
		self.MCE_version = version
		
		if max_num_configs is None:
			max_num_configs = num_configs
		
		if int(num_configs) != num_configs or int(max_num_configs) != max_num_configs or num_configs <= 0 or max_num_configs < num_configs:
			raise ValueError("num_configs and max_num_configs must be positive integers with max_num_configs >= num_configs")
		
		self.num_configs = int(num_configs)
		self.max_num_configs = int(max_num_configs)
		
		if psi_initial == "CDW":
			raise ValueError("MCE doesn't support CDW initial state, because multiple electrons would be required")
		
		if local_psi_init_select_eigenstate is None:
			local_psi_init_select_eigenstate = (version == 2 and initial_state == "bare")
		
		self.scale_el = scale_el
		
		self.local_psi_init_select_eigenstate = local_psi_init_select_eigenstate
		self.pancakeSampling = pancakeSampling
		self.pancakeScale = float(pancakeScale)
		self.flip_initial_state_for_half_configs = flip_initial_state_for_half_configs
		
		super(Simulation_MCE, self).__init__(hopping, omega, gamma,
		                                     num_sites_x, num_sites_y,
		                                     time_max, time_step,
		                                     initial_state, psi_initial, psi_state,
		                                     phonon_dressed_displaced_sites,
		                                     fullStorage,
		                                     force_new_simulation, load_sampled_state,
		                                     override_export_path, extra_subfolder,
		                                     calc_BO_electron_density, BO_nel_save_dist)
		
		if len(self.psi_state) != 1:
			raise ValueError("MCE only supports 1 electron")
		
		print(self.get_sim_name())
		if self.run_simulation:
			if    not hasattr(self, "start_hoz") or self.start_hoz is None\
			   or not hasattr(self, "start_el_coeff") or self.start_el_coeff is None\
			   or self.MCE_version == 2 and (not hasattr(self, "start_coupling_coeff") or self.start_coupling_coeff is None):
				# no stored simulation found, start a new one
				
				hoz, el_coeff, coupling_coeff, self.finalInitParams = \
				    initial_state_MCE(version, self.num_sites_x, self.num_sites_y,
				                      num_configs, max_num_configs,
				                      hopping, omega, gamma,
				                      psi_initial, self.psi_state[0], # [0] because only one electron allowed
				                      self.phonon_dressed_displaced_sites,
				                      local_psi_init_select_eigenstate=local_psi_init_select_eigenstate,
				                      printProgress=printProgress,
				                      returnFinalParams=True,
				                      seed=seed, allowCompression=allowCompression,
				                      allowExpansion=allowExpansion,
				                      scale_ph=scale_ph, scale_el=scale_el,
				                      pancakeSampling=pancakeSampling, pancakeScale=pancakeScale,
				                      flip_initial_state_for_half_configs=flip_initial_state_for_half_configs,
				                      calcCondition=calcCondition)
				
				# save initial placements in case they shall be stored later
				if self.MCE_version == 1:
					self.start_hoz = hoz[:, :, :self.num_configs].copy()
					self.start_el_coeff = el_coeff[0, :, :, :self.num_configs].copy()
					self.start_coupling_coeff = None
				else:
					self.start_hoz = hoz[:self.num_configs, :, :].copy()
					self.start_el_coeff = el_coeff[:self.num_configs, 0, :, :].copy()
					self.start_coupling_coeff = coupling_coeff[:self.num_configs].copy()
			else:
				if self.MCE_version == 1:
					hoz = np.empty((self.num_sites_y, self.num_sites_x, self.max_num_configs), dtype=np.complex)
					np.copyto(hoz[:, :, :self.num_configs], self.start_hoz, casting="no")
					
					el_coeff = np.empty((self.num_electrons, self.num_sites_y, self.num_sites_x, self.max_num_configs), dtype=np.complex)
					np.copyto(el_coeff[0, :, :, :self.num_configs], self.start_el_coeff, casting="no")
					
					coupling_coeff = None
				else:
					hoz = np.empty((self.max_num_configs, self.num_sites_y, self.num_sites_x), dtype=np.complex)
					np.copyto(hoz[:self.num_configs], self.start_hoz, casting="no")
					
					el_coeff = np.empty((self.max_num_configs, self.num_electrons, self.num_sites_y, self.num_sites_x), dtype=np.complex)
					np.copyto(el_coeff[:self.num_configs, 0], self.start_el_coeff, casting="no")
					
					coupling_coeff = self.start_coupling_coeff.copy()
			
			
			self.obs_y, self.obs_nph, self.obs_nel, self.obs_el_ekin, self.obs_el_ph_inter,\
			self.fullStorage_num_configs, self.fullStorage_el_coeff, self.fullStorage_hoz,\
			self.obs_nel_BO = \
			    simulate_system_MCE(hoz, el_coeff, coupling_coeff,\
			                        hopping, omega, gamma,\
			                        self.num_times, time_step, num_configs, allow_config_delete, allow_config_clone,\
			                        printProgress, fullStorage, calc_BO_electron_density, self.BO_nel_save_dist)
			
			if self.fullStorage_num_configs is not None and self.fullStorage_el_coeff is not None and self.fullStorage_hoz is not None:
				for t in range(self.fullStorage_num_configs.size):
					self.fullStorage_hoz     [t, :, :, self.fullStorage_num_configs[t]:] = np.nan
					self.fullStorage_el_coeff[t, :, :, self.fullStorage_num_configs[t]:] = np.nan
			
			if store_results:
				self.save_observables(override_export_path=override_export_path,
				                      extra_subfolder=extra_subfolder, silent=not printProgress)
				self.save_initial_binary(override_export_path=override_export_path,
				                         extra_subfolder=extra_subfolder, silent=not printProgress)
		else:	# not self.run_simulation
			if not hasattr(self, "fullStorage_num_configs"):
				self.fullStorage_num_configs = None
			self.fullStorage_el_coeff = None
			self.fullStorage_hoz = None
		
	
	
	def get_Wignerfunction_qp(self, time_index=None,
	                          q=None, pq=None, num_q=101, num_pq=101, symmetric_pq=True,
	                          autorange_sigma=2, autorange_all_times=False):
		""" calculates the Wigner function in q-pq phasespace at a given time index
		only works in the dimer (num_sites_x==2 and num_sites_y == 1) and with fullStorage at init
		only num_configs configurations are considered, cloning during simulation will lead to incorrect results
		
		time index: index of the time at which the Wigner function shall be calculated.
		            obtain the time as time()[time_index] (assuming time_index < self.cut).
		            If None (default), self.time() is used instead.
		
		q, pq: genarting edges (1d array) of full meshgrid at which the Wigner function shall be evaluated, defaults to None.
		       if None, automatic limits are determined, which can be controlled by following arguments
		       automatic grid is equally spaced.
		       full grid can be obtained by Q, PQ = np.meshgrid(q, pq)
		
		num_q, num_pq: number of points on each axis of generated meshgrid.
		               Only used if corresponding q or pq is None, otherwise overwritten with q.size (pq.size)
		
		symmetric_pq: if True, the generated intererval of pq is symmteric around 0 (|pq_max| = |pq_min|)
		
		autorange_sigma: number of standard deviations (sigma) padded around the
		                 convex hull of stored locations from phonons
		
		autorange_all_times: if True, the positions+momenta at all times are considered for the range.
		
		returns wigner_function, np.ndarray
		                         if len(time_index) == 1: shape (time_index.size, pq.size, q.size), dtype = np.float
		                         else: shape (pq.size, q.size), dtype = np.float
		"""
		
		if time_index is None:
			time_index = self.time()
		
		if self.num_sites_x != 2 or self. num_sites_y != 1:
			raise ValueError("q-pq space is only defined for the dimer(num_sites_x = 2, num_sites_y = 1). Selected at creation of this Simulation.")
		
		if self.fullStorage_hoz is None or self.fullStorage_el_coeff is None:
			raise RuntimeError("simulation wasn't run with fullStorage!")
		
		if self.fullStorage_num_configs is not None and len(np.unique(self.fullStorage_num_configs)) > 1:
			raise NotImplementedError("get_Wignerfunction_qp currently doesn't support a variable number of configurations over runtime (likely from cloning)")
		
		time_index = np.array(time_index, ndmin=1, dtype=np.int)
		if np.any(time_index < 0) or np.any(time_index >= self._time.size):
			raise ValueError("invalid time_index passed. (Every value) must be in range of [0, num_times).")
		
		
		if q is None or pq is None:
			# TODO allow for cloning?
			if not hasattr(self, "_qrange"):
				qrange_conf = self.fullStorage_hoz[:, 0, 0, :self.num_configs].real - self.fullStorage_hoz[:, 0, 1, :self.num_configs].real
				self._qrange = np.sqrt(2)*(np.hstack((qrange_conf.min(axis=1, keepdims=True), qrange_conf.max(axis=1, keepdims=True))))
			if not hasattr(self, "_pqrange"):
				pqrange_conf = self.fullStorage_hoz[:, 0, 0, :self.num_configs].imag - self.fullStorage_hoz[:, 0, 1, :self.num_configs].imag
				self._pqrange = np.sqrt(2)*(np.hstack((pqrange_conf.min(axis=1, keepdims=True), pqrange_conf.max(axis=1, keepdims=True))))
		
		if q is None:
			# determine range
			if autorange_all_times:
				q_min = np.min(self._qrange[:, 0])
				q_max = np.max(self._qrange[:, 1])
			else:
				q_min = np.min(self._qrange[time_index, 0])
				q_max = np.max(self._qrange[time_index, 1])
			
			q = np.linspace(q_min - autorange_sigma, q_max + autorange_sigma, num_q)
		else:
			q = np.array(q).flatten()
			num_q = q.size
		
		if pq is None:
			# determine range
			if autorange_all_times:
				pq_min = np.min(self._pqrange[:, 0])
				pq_max = np.max(self._pqrange[:, 1])
			else:
				pq_min = np.min(self._pqrange[time_index, 0])
				pq_max = np.max(self._pqrange[time_index, 1])
			
			if symmetric_pq:
				pq_max = max(abs(pq_min), abs(pq_max))
				pq_min = -pq_max
			
			pq = np.linspace(pq_min - autorange_sigma, pq_max + autorange_sigma, num_pq)
		else:
			pq = np.array(pq).flatten()
			num_pq = pq.size
		
		
		Q, PQ = np.meshgrid(q, pq)
		wigner = np.zeros((time_index.size, pq.size, q.size))
		
		#prefact[k, j], j-> shape(1, num_configs), k-> shape(num_configs, 1)
		prefact = np.empty((self.num_configs, self.num_configs), dtype=np.complex)
		
		
		# TODO calculate wigner
		for ti in range(len(time_index)):
			t = time_index[ti]
			prefact *= 0
			
			prefact += ( (  self.fullStorage_hoz[t, 0, 0, :self.num_configs].real
			              + self.fullStorage_hoz[t, 0, 1, :self.num_configs].real).reshape(1, self.num_configs)
			            -(  self.fullStorage_hoz[t, 0, 0, :self.num_configs].real
			              + self.fullStorage_hoz[t, 0, 1, :self.num_configs].real).reshape(self.num_configs, 1))**2\
			          +( (  self.fullStorage_hoz[t, 0, 0, :self.num_configs].imag
			              + self.fullStorage_hoz[t, 0, 1, :self.num_configs].imag).reshape(1, self.num_configs)
			            -(  self.fullStorage_hoz[t, 0, 0, :self.num_configs].imag
			              + self.fullStorage_hoz[t, 0, 1, :self.num_configs].imag).reshape(self.num_configs, 1))**2
			prefact /= -4
			
			# all 4 phase terms:
			#prefact += 1j*(-self.fullStorage_hoz[t, 0, 0, :self.num_configs].real.reshape(1, self.num_configs)*self.fullStorage_hoz[t, 0, 1, :self.num_configs].imag.reshape(self.num_configs, 1)\
			#               -self.fullStorage_hoz[t, 0, 1, :self.num_configs].real.reshape(1, self.num_configs)*self.fullStorage_hoz[t, 0, 0, :self.num_configs].imag.reshape(self.num_configs, 1)\
			#               +self.fullStorage_hoz[t, 0, 0, :self.num_configs].imag.reshape(1, self.num_configs)*self.fullStorage_hoz[t, 0, 1, :self.num_configs].real.reshape(self.num_configs, 1)\
			#               +self.fullStorage_hoz[t, 0, 1, :self.num_configs].imag.reshape(1, self.num_configs)*self.fullStorage_hoz[t, 0, 0, :self.num_configs].real.reshape(self.num_configs, 1))
			
			# simpler phase terms by using symmetry:
			prefact += 1j*(+self.fullStorage_hoz[t, 0, 0, :self.num_configs].imag.reshape(1, self.num_configs)*self.fullStorage_hoz[t, 0, 1, :self.num_configs].real.reshape(self.num_configs, 1)\
			               +self.fullStorage_hoz[t, 0, 1, :self.num_configs].imag.reshape(1, self.num_configs)*self.fullStorage_hoz[t, 0, 0, :self.num_configs].real.reshape(self.num_configs, 1))
			prefact.imag -= prefact.imag.T
			
			prefact = np.exp(prefact)
			
			# multiply with electronic overlap coefficient
			prefact *= np.sum(np.conj(self.fullStorage_el_coeff[t, :, :, :self.num_configs]).reshape(self.num_sites_y, self.num_sites_x, self.num_configs, 1)*\
			                          self.fullStorage_el_coeff[t, :, :, :self.num_configs] .reshape(self.num_sites_y, self.num_sites_x, 1, self.num_configs), axis=(0, 1))
			
			prefact /= 2*np.pi
			
			for k in range(self.num_configs):
				wigner[ti] += prefact[k, k].real*np.exp(-.5*(  ( Q-np.sqrt(2)*(self.fullStorage_hoz[t, 0, 0, k].real - self.fullStorage_hoz[t, 0, 1, k].real))**2\
				                                             + (PQ-np.sqrt(2)*(self.fullStorage_hoz[t, 0, 0, k].imag - self.fullStorage_hoz[t, 0, 1, k].imag))**2))
				for j in range(0, k):
					exp_real  = ( Q- ((  self.fullStorage_hoz[t, 0, 0, j].real - self.fullStorage_hoz[t, 0, 1, j].real
					                   + self.fullStorage_hoz[t, 0, 0, k].real - self.fullStorage_hoz[t, 0, 1, k].real)/np.sqrt(2)))**2
					exp_real += (PQ- ((  self.fullStorage_hoz[t, 0, 0, j].imag - self.fullStorage_hoz[t, 0, 1, j].imag
					                   + self.fullStorage_hoz[t, 0, 0, k].imag - self.fullStorage_hoz[t, 0, 1, k].imag)/np.sqrt(2)))**2
					exp_real *= -.5
					
					exp_imag  =  Q*((  self.fullStorage_hoz[t, 0, 0, j].imag - self.fullStorage_hoz[t, 0, 1, j].imag
					                 - self.fullStorage_hoz[t, 0, 0, k].imag + self.fullStorage_hoz[t, 0, 1, k].imag)/np.sqrt(2))
					exp_imag -= PQ*((  self.fullStorage_hoz[t, 0, 0, j].real - self.fullStorage_hoz[t, 0, 1, j].real
					                 - self.fullStorage_hoz[t, 0, 0, k].real + self.fullStorage_hoz[t, 0, 1, k].real)/np.sqrt(2))
					
					wigner[ti] += 2*(prefact[k, j]*np.exp(exp_real + 1j*exp_imag)).real
		
		
		if time_index.size == 1:
			return q, pq, wigner[0]
		return q, pq, wigner
	
	
	def get_num_configs_over_time(self):
		if self.fullStorage_num_configs is not None:
			return self.fullStorage_num_configs
		else:
			return np.full(self._time.size, self.num_configs)
	
	
	def get_config_highest_occupied_eigenstate(self, time_index):
		if self.fullStorage_hoz is None or self.fullStorage_el_coeff is None:
			raise RuntimeError("simulation wasn't run with fullStorage!")
		
		if self.num_sites_x != 1 and self.num_sites_y != 1:
			raise NotImplementedError("only 1D systems BO occupation of each state supported")
		
		num_configs = self.get_num_configs_over_time()[time_index]
		x = self.fullStorage_hoz[time_index, : , :, :num_configs].real.reshape(self.num_sites, num_configs).T
		el_coeff = self.fullStorage_el_coeff[time_index, : , :, :num_configs].reshape(self.num_sites, num_configs).T
		
		# TODO only works in 1D
		H = np.zeros((num_configs, self.num_sites, self.num_sites))
		diag_ind = np.arange(self.num_sites)
		
		H[:, diag_ind[:-1], diag_ind[1:]] = -self.hopping
		H[:, diag_ind[1:], diag_ind[:-1]] = -self.hopping
		H[:, diag_ind, diag_ind] = (-np.sqrt(2)*self.gamma)*x
		
		BO_eigvals, BO_eigvecs = np.linalg.eigh(H)
		BO_coeffcients = np.conjugate(BO_eigvecs.transpose(0, 2, 1)) @ el_coeff.reshape(num_configs, self.num_sites, 1)
		BO_coeffcients = BO_coeffcients.reshape(num_configs, self.num_sites)
		# sqrt(occupation probability) along each axis, nut normalized
		BO_prob = np.abs(BO_coeffcients)
		
		return np.argmax(BO_prob, axis=1)
		

