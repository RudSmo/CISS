module InputFile

include("rk4.jl")
include("hamiltonian.jl")
using .RK4
using .Hamiltonian


el_state = [1,1,0,0,0,0]
ph_state = [0,0,0]
vibStates = 1
num_sites = Int(length(el_state)/2)
num_laps = 1
sites_lap = num_sites
t0 = 1.0
eps = [1.0,-1.0]
tPh = 0.0
lambda0 = 0.0
lambdaPh = 0.0*lambda0
wPh = 1.4*t0
r_hel = 1.0
c_hel = 1.0
psi0 = InitialState(el_state,ph_state,vibStates,num_sites)
H = Hamiltonian.Fransson(num_sites,num_laps,vibStates,eps,t0,tPh,lambda0,lambdaPh,wPh,r_hel,c_hel)


dt = 0.001
path = "test_data/"
if !isdir(path)
    mkdir(path)
end
filname = "$(path)pops"
imagfilname = "$(path)test_1site_im"
pathIn = "$(path)final_Quantities"

nsteps = Int(100/dt)
insteps = Int(100/dt)
nsave = 100
nflush= 1000

Propagate(psi0,num_sites,vibStates,H,dt,filname,nsteps,nsave,nflush)
#PropagateImag(num_sites,vibStates,H,dt,filname,imagfilname,nsteps,insteps,nsave,nflush)




end
