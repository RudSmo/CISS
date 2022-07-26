module InputFile

include("rk4.jl")
include("hamiltonian.jl")
using .RK4
using .Hamiltonian

el_state = [1,0,0,0]
ph_state = [1,0]
vibStates = 5
num_sites = Int(length(el_state)/2)
num_laps = 1
sites_lap = 2
eps = [0.0,0.0]
t0 = 1.0
tPh = 0.1
lambda0 = 1.0
lambdaPh = 0.1
r_hel = 1.0
c_hel = 1.0
psi0 = InitialState(el_state,ph_state,vibStates,num_sites)

H = Hamiltonian.Fransson(num_sites,num_laps,vibStates,eps,t0,tPh,lambda0,lambdaPh,r_hel,c_hel)
dt = 0.001
path = "data/"
if !isdir(path)
    mkdir(path)
end
filname = "$(path)test_1site.dat"
imagfilname = "$(path)test_1site_im.dat"

nsteps = Int(50/dt)
insteps = Int(10/dt)
nsave = 20
nflush= 400

#Propagate(psi0,num_sites,vibStates,H,dt,filname,nsteps,nsave,nflush)
PropagateImag(num_sites,vibStates,H,dt,filname,imagfilname,nsteps,insteps,nsave,nflush)




end