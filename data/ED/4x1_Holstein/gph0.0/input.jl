module InputFile
src_path=expanduser("~/PhD/projects/CISS/src/exact_TDSE/")
include("$(src_path)rk4.jl")
include("$(src_path)hamiltonian.jl")
using .RK4
using .Hamiltonian

el_state = [1,1,0,0]
ph_state = [0,0]
vibStates = 1
num_sites = Int(length(el_state)/2)
num_laps = 1
sites_lap = num_sites
t0 = 1.0
eps = zeros(2)
gPh = 0.0
tPh = 0.0
lambda0 = 0.0
lambdaPh = 0.0
wPh = 1.0
r_hel = 1.0
c_hel = 1.0
psi0 = InitialState(el_state,ph_state,vibStates,num_sites)

H = Hamiltonian.Holstein(num_sites,vibStates,eps,t0,gPh,wPh)
dt = 0.001
#path = ""
#if !isdir(path)
#    mkdir(path)
#end
filname = "pops_vib.dat"
imagfilname = "pop_vib_im.dat"
pathIn = ""

nsteps = Int(50/dt)#Int(200/dt)
insteps = Int(10/dt)
nsave = 10
nflush= 60

@timev Propagate(psi0,num_sites,vibStates,H,dt,filname,nsteps,nsave,nflush,pathIn)
#PropagateImag(num_sites,vibStates,H,dt,filname,imagfilname,nsteps,insteps,nsave,nflush)




end
