from constants import *
import numpy as np
import soc

eps=1.e-12 # accuracy for tensor train rounding

# The number of electroc DoFs in the system part:
M_laps  = 1
N_ions_per_lap = 3
Dim_elecDOF = M_laps * N_ions_per_lap   
Spin = True
if Spin: 
    N_spin=2
else:
    N_spin=1
Spin_alt_order = True
# True:  spin-up (1st DoF), spin-down (1st DoF), spin-up (2nd DoF), spin-down (2nd DoF), spin-up (3rd DoF), spin-down (3rd DoF), ...
# False: spin-up (1st DoF), spin-up (2nd DoF), spin-up (3rd DoF),..., spin-down (1st DoF), spin-down (2nd DoF), spin-down (3rd DoF), ...

# The number of vibrational DoFs in the system part
Dim_vibDOF = Dim_elecDOF    
# The size of basis set for vibrational DoFs in the system part. For example, [20]*Dim_vibDOF means taking into account 20 vibrational states for each DoF. It can also be [20, 30, 25,...] for the first, second, third, ... DoF. 
vibstates=  5
vibbasis = [vibstates]*Dim_vibDOF  

# Number of Pade poles for Fermi function 
N_fermipoles = 1
# Number of leads
N_leads = 2

head_tail_coup = True
# the type of environments: fermi
# (head-tail-only) means only the head and tail sites are coupled to left and right lead, respectively.
Outputs="current, n_vib" 
#Outputs="current,elec_pop[pop_el_states],vib_pop{pop_vib_states}, n_vib" 
# observables to be calcualted and printed out: 'current' or 'Current', 
#                                               'elec_pop[all/file_name]' or 'Pop_el[all/file_name]', 
#                                               'n_vib' or 'N_vib', for average vibrational excitation
#                                               'vib_pop{all/file_name}' or 'Pop_vib{all/file_name}' 
# for example, for Outputs="elec_pop[all], current, n_vib, vib_pop{pop_vib_states}" 
# we choose to print the population of all electronic states, electronic current, average vibrational excitation, population of some vibrational states(specified in the file pop_vib_states.out) 
# format of output files: 
# 'current': time, I_L, I_R
# 'elec_pop': time, pop_of electronic states...
# 'vib_pop': time, pop_of vibrational states...
# 'n_vib': time, <b_1^+ b_1>, <b_2^+ b_2>,...

Init_state="initial_state"
# "ground_state" or ground_state: 0,0,0,0,0,0,.....
# "initial_state:  specified initial state in system subspace (in Hilbert not in twin space) in the file initial_state.in 
# for example: when Dim_elecDOF = 2, N_spin =2, Dim_vibDOF=2, 
#  1, 0, 0, 0, 2, 0 in initial_state.in means the first electronic DOF with spin up is occupied at the beginning, other electronic dofs are unoccupied. The first vibrational DoF is excited initially to the second excited state and the second vib DoF remains in the ground state.
# "equilibrium": the equilibrium state of the system+environment complex, obtained through imaginary time propagtion.  

Dt = 0.001
Tend = 200
Nsteps = int(Tend/Dt)
dt = Dt/au2fs# in units of femtosecond
Nsections = 10000

Init_propagation = True
t_in=0./au2fs
# 'Init_propagation = True' if we start from t_0=0, 'False' if we resume the simulation starting from t_0 = t_in (given below)


t0 = 1.0 
# Vibrational frequencies
vib_freq = np.zeros((Dim_vibDOF))
for i_vib in range(Dim_vibDOF):
    vib_freq[i_vib] =  t0*1.4/au2eV  # in units of eV

# Vibronic Coupling
vib_elec_coup=np.zeros( (Dim_elecDOF, N_spin, Dim_vibDOF) )
for i_elec in range(Dim_elecDOF):
    for i_spin in range(N_spin):
        for i_vib in range(Dim_vibDOF):
            if (i_vib==i_elec): vib_elec_coup[i_elec][i_spin][i_vib]=0.0/au2eV  
# Here, it is assumed that every electronic DoF is coupled to one vibrational DoF. But it is not necessary the case. The code also allows the cross couplings. 


# Electronic Hamiltonian
# H = \sum_{i} onsite_energy[i] d_i^{dagger} d_i 
#   + \sum_{i < j} (exchange_energy[i,j] d_i^{dagger} d_j +h.c.) 
#   + \sum_{i < j} Couloumb_energy[i,j] d_i^{dagger} d_i d_j^{dagger} d_j
onsite_energy   = np.array([[1.0,-1.0],[1.0,-1.0],[1.0,-1.0]])#np.full((Dim_elecDOF, N_spin), e_onsite/au2eV, dtype=float) 
Couloumb_energy = np.zeros((Dim_elecDOF, N_spin, Dim_elecDOF, N_spin),dtype=float) 
#Couloumb_energy[0,0,1,0] =  15/au2eV

exchange_energy = np.zeros((Dim_elecDOF, N_spin, Dim_elecDOF, N_spin),dtype=float) 
for i_elec in range(Dim_elecDOF-1):
    for i_spin in range(N_spin):
        exchange_energy[i_elec,i_spin,i_elec+1,i_spin] = -t0/au2eV

# Fransson does not want to tell us what parameters he is using.
# Therefore, I am setting them arbitrarily.  
# chirality = "+" or "-" for Right or Left-handedness,
# r_helix = radius of helix,
# l_helix = length of the helix,
chirality = "+"
r_helix = 1.0 
l_helix = 1.0
so_coup = soc.spin_orbital_coupling(r_helix, l_helix, M_laps, Dim_elecDOF, chirality)

# Fransson's next-nearest-neighbor spin-orbit coupling
spin_orbit_coupling = np.zeros((Dim_elecDOF,N_spin,N_spin),dtype=complex)
# lambda0 = SOC-strength.
lambda0 = t0/40/au2eV
for i_elec in range(Dim_elecDOF-2):
    vsi = so_coup.VSI(i_elec)
    for i_spin in range(N_spin): 
        for j_spin in range(N_spin):
            spin_orbit_coupling[i_elec,i_spin,j_spin] = 1j * lambda0 * vsi[i_spin,j_spin]

# nearest-neighbor electron-phonon coupling
nnhopping_coupling_ph = np.zeros((Dim_elecDOF,N_spin,N_spin),dtype=float)
tPh = 0.1*t0/au2eV
for i_elec in range(Dim_elecDOF-1):
    for i_spin in range(N_spin):  
        for j_spin in range(N_spin):
            nnhopping_coupling_ph[i_elec,i_spin,j_spin] = tPh

# next-nearest-neighbor spin-dependent electron-phonon coupling
spin_orbit_coupling_ph = np.zeros((Dim_elecDOF,N_spin,N_spin),dtype=complex)
lambdaPh = 0.0/au2eV
for i_elec in range(Dim_elecDOF-2):
    vsi = so_coup.VSI(i_elec)
    for i_spin in range(N_spin):  
        for j_spin in range(N_spin):
            spin_orbit_coupling_ph[i_elec,i_spin,j_spin] = 1j * lambdaPh * vsi[i_spin,j_spin]

temperature = 0.1/au2eV # in units of eV
beta = 1./temperature
Ntaus = 50
dtau = beta /Ntaus
temperature_fermibath = np.full((N_leads, N_spin), temperature)

Phi = 0 * 0.1/au2eV   # Bias Voltage, in units of V
Chem_potentials = np.full( (N_leads, N_spin), 0.)
for i_lead in range(N_leads):
    Chem_potentials[i_lead,0] = (-1)**(i_lead)*Phi/2
    Chem_potentials[i_lead,1] = (-1)**(i_lead)*Phi/2

WBL=False
broadenings = 10.0/au2eV
broadenings_fermibath = np.full((N_leads, N_spin), broadenings)
# 'WBL=True' to invoke wide-band approximation. If it is 'False', need to set the broadenings (in units of eV) 

# Coupling between electronic DOFs and leads (in units of eV)
coupling=0.0
if head_tail_coup:
    elec_lead_couplings=np.zeros((N_spin, N_leads))
    for i_spin in range(N_spin):
        for i_lead in range(N_leads):
            elec_lead_couplings[i_spin][i_lead]=np.sqrt(coupling)
else:
    #elec_lead_couplings=np.full((Dim_elecDOF, N_spin, N_leads),np.sqrt(coupling))
    elec_lead_couplings=np.zeros((Dim_elecDOF, N_spin, N_leads))
    for i_elec in range(Dim_elecDOF):
        for i_spin in range(N_spin):
            for i_lead in range(N_leads):
                if (i_lead==i_elec):
                    elec_lead_couplings[i_elec,i_spin,i_lead]=np.sqrt(coupling)
                else:
                    elec_lead_couplings[i_elec,i_spin,i_lead]=np.sqrt(coupling)

