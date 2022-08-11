module InputFile
src_path=expanduser("~/PhD/projects/CISS/src/exact_TDSE/")
include("$(src_path)rk4.jl")
include("$(src_path)hamiltonian.jl")
using .RK4
using .Hamiltonian

el_state = [1,1,0,0,0,0]
ph_state = [0,0,0]
vibStates = 5 
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
#path = ""
#if !isdir(path)
#    mkdir(path)
#end
filname = "pops_vib_2.dat"
imagfilname = "pop_vib_im_2.dat"
pathIn = "pops_vib.dat_final_Quantities"

nsteps = Int(950/dt)#Int(200/dt)
insteps = Int(10/dt)
nsave = 100
nflush= 600

@timev Propagate(psi0,num_sites,vibStates,H,dt,filname,nsteps,nsave,nflush,pathIn)
#PropagateImag(num_sites,vibStates,H,dt,filname,imagfilname,nsteps,insteps,nsave,nflush)




end
