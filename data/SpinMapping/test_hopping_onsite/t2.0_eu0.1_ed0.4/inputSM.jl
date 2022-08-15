module Input

src_path =  expanduser("~/PhD/projects/CISS/src/SpinMapping/NLS_Phonons/")
include("$(src_path)hamiltonian.jl")
include("$(src_path)structures.jl")
include("$(src_path)mapping_NLS_phonons.jl")
using .HamiltonianLib
using .Structures
using .MappingNLSPhonons

L      = 4
Nlvl   = 2*L
t0     = 2.0
ES1    = 0.1
ES2    = 0.4
E      = [ES1,ES2]
tso    = 0.0
a      = 1.0
c      = 1.0
num_laps = 1
tph    = 0.0
tsoph  = 0.0
gamma  = zeros(L)
omega  = ones(L)


H      = FranssonHamiltonian(L,E,t0,tso,a,c,num_laps,tph,tsoph,omega,omega)
beta   = 100
sphere = "W"
state  = [1,1,0,0,0,0,0,0]

ntraj  = 1000
dt     = 0.001
tend   = 200
Nstp   = Int(tend/dt)
Nsave  = 1
Nflush = 2
fn = "dataSM"

if !isdir(fn)
    mkdir(fn)
end


fn = "$(fn)/$(sphere)_S1$(ES1)_S2$(ES2)_ntraj$(ntraj)_dt$(dt)"
fnE = "$(fn)_electrons.dat"
fnP = "$(fn)_phonons.dat"

SimPar = SimulationParameters(Nlvl,sphere,ntraj,dt,Nstp,Nsave,Nflush,fnE,fnP)
HamPar = HamiltonianParameters(H,gamma,omega,beta)
PSVar  = FocusedInitialCondition(SimPar,HamPar,state)

NSimulate(HamPar,SimPar,state)

end
