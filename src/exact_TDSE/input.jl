module InputFile

include("rk4.jl")
include("hamiltonian.jl")
using .RK4
using .Hamiltonian


el_state = [1,0,0,0]
ph_state = [0,0]
vibStates = 3
num_sites = Int(length(el_state)/2)
num_laps = 1
sites_lap = num_sites
t0 = 0.04
eps = -6*t0 * ones(2)
tPh = 0.1*t0
lambda0 = t0/40
lambdaPh = 0.1*lambda0
wPh = 1.4*t0
r_hel = 1.0
c_hel = 1.0
psi0 = InitialState(el_state,ph_state,vibStates,num_sites)

H = Hamiltonian.Fransson(num_sites,num_laps,vibStates,eps,t0,tPh,lambda0,lambdaPh,wPh,r_hel,c_hel)
dt = 0.001
path = "data/"
if !isdir(path)
    mkdir(path)
end
filname = "$(path)test_1site_1.dat"
imagfilname = "$(path)test_1site_im.dat"
pathIn = "$(path)test_1site.dat_final_Quantities"

nsteps = Int(50/dt)
insteps = Int(100/dt)
nsave = 200
nflush= 1000

Propagate(psi0,num_sites,vibStates,H,dt,filname,nsteps,nsave,nflush,pathIn)
#PropagateImag(num_sites,vibStates,H,dt,filname,imagfilname,nsteps,insteps,nsave,nflush)




end
