module InputHamiltonian
include("hamiltonian.jl")
using .Hamiltonian


# Convention for constructing Hamiltonian
# Hilbert space: H = H_1u ⊗ H_1d ⊗ ... H_Lu ⊗ H_Ld ⊗

num_sites = 3                   # Number of sites 
num_laps  = 1                   # Number of laps
vibStates = 1                   # Number of vibrational states per mode
eps       = [0.,0.]             # On-site energy [E_up,E_dw]
t0        = 1.0                 # Elastic hopping-parameter
tPh       = 0.0                 # Inelastic hopping, coupling to vibrational modes
lambda0   = 0.0                 # Spin-Orbit-coupling  
lambdaPh  = 0.0                 # Non-static SOC 
wPh       = 0.0                 # Vibrational frequency for each mode
r_helix   = 1.0                 # Radius of helix 
c_helix   = 1.0                 # Length of helix

H = Fransson(num_sites,num_laps,vibStates,eps,t0,tPh,lambda0,lambdaPh,wPh,r_helix,c_helix)
display(H)

end