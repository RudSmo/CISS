module InputFile

include("rk4.jl")
include("hamiltonian.jl")
using .RK4
using .Hamiltonian


state = [0,1,1,0,0,0,0,0,0,0]
num_sites = Int(length(state)/2-1)
t0 = 1.0
JK = 1.0
U = 0.0
psi0 = RK4.ElectronicStateNI(state)
H = Hamiltonian.KondoImpurityNI(num_sites,t0,JK)

dt = 0.001
path = "data_Kondo_compare/NI/"
if !isdir(path)
    mkpath(path)
end
filname = "$(path)L$(num_sites)_U$(U)_JK$(JK)"
imagfilname = "$(path)test_1site_im"
pathIn = "$(path)final_Quantities"

nsteps = Int(100/dt)
insteps = Int(100/dt)
nsave = 25
nflush= 1000

PropagateKondoNI(psi0,num_sites,H,dt,filname,nsteps,nsave,nflush)

end
