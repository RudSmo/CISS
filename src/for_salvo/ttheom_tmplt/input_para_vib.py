from constants import *
import numpy as np

eps=1.e-12 # accuracy for tensor train rounding

# The number of electroc DoFs in the system part:
n_dummy = 0
Dim_elecDOF = 1 
pole_scheme = 'pade' # choose 'pade' or 'barycentric' 
Spin = False
if Spin: 
    N_spin=2
else:
    N_spin=1
Spin_alt_order = True
N_sign = 2
# True:  spin-up (1st DoF), spin-down (1st DoF), spin-up (2nd DoF), spin-down (2nd DoF), spin-up (3rd DoF), spin-down (3rd DoF), ...
# False: spin-up (1st DoF), spin-up (2nd DoF), spin-up (3rd DoF),..., spin-down (1st DoF), spin-down (2nd DoF), spin-down (3rd DoF), ...

# The number of vibrational DoFs in the system part
Dim_vibDOF = 1
# The size of basis set for vibrational DoFs in the system part. For example, [20]*Dim_vibDOF means taking into account 20 vibrational states for each DoF. It can also be [20, 30, 25,...] for the first, second, third, ... DoF. 
vibstates=  15
vibbasis = [vibstates]*Dim_vibDOF  
# Number of Pade poles for Fermi function 
N_fermipoles = 6
# Number of leads
N_leads = 2

head_tail_coup = True
# the type of environments: fermi
# (head-tail-only) means only the head and tail sites are coupled to left and right lead, respectively.
Outputs="current, n_vib" 
#Outputs="current,elec_pop[pop_el_states],vib_pop{pop_vib_states}, n_vib" 

#                                               'elec_pop[all/file_name]' or 'Pop_el[all/file_name]', 
#                                               'n_vib' or 'N_vib', for average vibrational excitation
#                                               'vib_pop{all/file_name}' or 'Pop_vib{all/file_name}' 
# for example, for Outputs="elec_pop[all], current, n_vib, vib_pop{pop_vib_states}" 

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

# Subotniks parameters
# Delta = 0
# A = B = 1
# kappa_x = 0
# chi = 1
# kappa_y = 0.8
# Gamma = 1
G0 = 0.01

lambda_1 = 0.1/np.sqrt(2)
lambda_2 = 0.2/np.sqrt(2)
t_hop = 0.0
Ei = np.zeros((2,N_spin))
Tij = np.zeros((2,2,N_spin))
for i in range(N_spin):
    Ei[0,i] = 0.15
#    Ei[1,i] = -0.15
#    Tij[0,1,i] = t_hop

vib_elec_coup = np.zeros((2,2,N_spin,Dim_vibDOF),dtype=complex)
for i in range(N_spin):
    vib_elec_coup[0,0,i,0] = lambda_1
 #   vib_elec_coup[1,1,i,0] = -lambda_1
 #   vib_elec_coup[0,0,i,1] = lambda_2
 #   vib_elec_coup[1,1,i,1] = -lambda_2


#vib_freq = np.array([0.1,0.1])
vib_freq = np.array([0.1])
#kappas = np.array([0.0,0.0])
kappas = np.array([0.0])


Dt = 0.05
Tend = 0.05

Nsteps = int(Tend/Dt)
dt = Dt# in units of femtosecond
Nsections = Nsteps

Init_propagation = True
t_in=0 # ./au2fs
# 'Init_propagation = True' if we start from t_0=0, 'False' if we resume the simulation starting from t_0 = t_in (given below)

temperature = 0.0258 # in units of eV
beta = 1./temperature
Ntaus = 50
dtau = beta /Ntaus
temperature_fermibath = np.full((N_leads, N_spin), temperature)

Phi = 0 * 0.01 # * 0.05    # Bias Voltage, in units of V
Chem_potentials = np.full( (N_leads, N_spin), 0.)
for i_lead in range(N_leads):
    Chem_potentials[i_lead] = (-1)**(i_lead)*Phi/2


W = 10.0
WBL = True
fermi_range = 1.5                                                         # Set range of Fermi-Dirac function to decompose (Crucial for large energies!)i
analytic_spectral_function_decomposition = True                            # Leave as true for Lorentzians (we can decompose them analytically)
if WBL == 0:                                                               # If not under the wide-band limit, they are assumed to have a Lorentzian density of states
    specwidth = np.array([W],dtype=float)                                  # Set spectral width of bath Lorentzians (dummy variable in WBL)
    symmetrized_fermi_specwidth = np.array([fermi_range],dtype=float)/au2eV
    Ncutoff = 1                                                            # Leave as 1
elif WBL == 1:
    specwidth = np.array([W,W],dtype=float)
    symmetrized_fermi_specwidth = np.array([fermi_range],dtype=float)/au2eV
    Ncutoff = 1
# Barycentric decomposition parameters
Nsupport_points_barycentric = 2000                                         # Number of support points used in barycentric decomposition. Can safely leave as 2000
tol_F = 1e-6                                                               # Minimum spacing between support points. Leave as 1e-6
tol_Gamma_barycentric = 1e-6                                               # Tolerance of decomposition of Gamma spectral functions - leave as 1e-5 because
                                                                           # we do this part analytically for Lorentzians
tol_fermi_symmetrized_barycentric = 1e-3

broadenings = 1000.0 #/au2eV
broadenings_fermibath = np.full((N_leads, N_spin), broadenings)
# 'WBL=True' to invoke wide-band approximation. If it is 'False', need to set the broadenings (in units of eV) 

# Coupling between electronic DOFs and leads (in units of eV)
P = 0.0
coupling = G0

#elec_lead_couplings=np.zeros((Dim_elecDOF, N_spin, N_leads))
#elec_lead_couplings[0,0,0]=np.sqrt(coupling)
#elec_lead_couplings[1,0,1]=np.sqrt(coupling)

if head_tail_coup:
    elec_lead_couplings=np.zeros((N_spin, N_leads),dtype=complex)
    for i_spin in range(N_spin):
        for i_lead in range(N_leads):
            elec_lead_couplings[i_spin][i_lead]=np.sqrt(coupling)
#else:
#    #elec_lead_couplings=np.full((Dim_elecDOF, N_spin, N_leads),np.sqrt(coupling))
#    elec_lead_couplings=np.zeros((Dim_elecDOF, N_spin, N_leads))
#    for i_elec in range(Dim_elecDOF):
#        for i_spin in range(N_spin):
#            for i_lead in range(N_leads):
#                if (i_lead==i_elec):
#                    elec_lead_couplings[i_elec,i_spin,i_lead]=np.sqrt(coupling)
#                else:
#                    elec_lead_couplings[i_elec,i_spin,i_lead]=np.sqrt(coupling)
#
