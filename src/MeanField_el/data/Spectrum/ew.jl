module InputFile

#src_path=expanduser("~/Libraries/CISS/src/MeanField/")
#include("$(src_path)hamiltonian.jl")
#include("$(src_path)iosave.jl")
#include("$(src_path)multiopenint.jl")
#include("$(src_path)classical_phonons.jl")
include("../../hamiltonian.jl")
include("../../iosave.jl")
include("../../multiopenint.jl")

using .Hamiltonian
using .IOsave
using LinearAlgebra
using .MultiOpenIntegrator
using DelimitedFiles
# -------------- Parameters --------------------------

# Simulation Parameters
dt = 0.001
tend = 1000.0
numstep = Int(tend/abs(dt))
writeMode = "minimal"
initProp  = true
# true = propagate as usuals
# false = use initial density matrix from previous computation

# Geometry Parameters
t0 = 1.0
eu = 0.0
ed = 0.0
E = [eu,ed]
a = 1.0
c = 1.0
num_laps = 20
ions = 4
l = Int(num_laps * ions)
tph = 0.0
tsoph = 0.0
ws = ones(l)*1.4*t0
qs = zeros(l)


Tso = [1.0]#Array(0.0:0.01:20.0)
fil = open("eigvals_laps$(num_laps)_ions$(ions)_eu$(eu)_ed$(ed)_a$(a)_c$(c)_tsoVar.dat","w")
for tso in Tso
	H = Hamiltonian.FranssonHamiltonian(l,E,t0,tso,a,c,num_laps,tph,tsoph,ws,qs)
	print(H == H')
	ew = eigvals(H)
	write(fil,"$(tso) ")
	for j in 1:length(ew)
		write(fil,"$(ew[j]) ")
	end
	write(fil,"\n")
end
close(fil)


end

