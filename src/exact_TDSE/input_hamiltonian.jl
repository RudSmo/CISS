module InputHamiltonian
include("hamiltonian.jl")
using .Hamiltonian
using PyPlot
using LinearAlgebra

# Convention for constructing Hamiltonian
# Hilbert space: H = H_1u ⊗ H_1d ⊗ ... H_Lu ⊗ H_Ld ⊗

num_sites = 3                   # Number of sites 
num_laps  = 1                  # Number of laps
vibStates = 3                   # Number of vibrational states per mode
T0 	  = 0.04
eps       = [-6*T0,-6*T0]             # On-site energy [E_up,E_dw]
t0        = -T0                 # Elastic hopping-parameter
tPh       = 0.1*T0               # Inelastic hopping, coupling to vibrational modes
lambda0   = T0/40                 # Spin-Orbit-coupling  
#lambdaPh  = 0.1*lambda0                 # Non-static SOC 
wPh       = 0.002                 # Vibrational frequency for each mode
r_helix   = 1.0                 # Radius of helix 
c_helix   = 1.0                 # Length of helix



E = -40*T0:0.0005*T0:5*T0

function Dos(E_i,E,b=1e-4)
    D = [1/(sqrt(2*pi)*b)*sum([exp(-(e-i)^2/(2*b^2)) for i in E_i]) for e in E]
    return D
end

for lambdaPh in [0.1*lambda0,10.0*lambda0]
	

	H = Fransson(num_sites,num_laps,vibStates,eps,t0,tPh,lambda0,lambdaPh,wPh,r_helix,c_helix)
	ew,ev = eigen(H)



	D = Dos(ew,E)
	plt.plot(E,D,label=L"$\lambda_1=\,$"*"$(lambdaPh)")
	plt.legend()
	println(ev[:,1])
end
plt.show()	




end
