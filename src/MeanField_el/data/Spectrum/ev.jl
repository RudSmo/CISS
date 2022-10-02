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
using PyPlot
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
ions = 8
l = Int(num_laps * ions)
tph = 0.0
tsoph = 0.0
ws = ones(l)*1.4*t0
qs = zeros(l)
Tsos = Array(0.0:0.005:3.0)
Nd = zeros(length(Tsos),l)
ev_ind = Int(l)
for (ind,tso) in enumerate(Tsos)
	H = Hamiltonian.FranssonHamiltonian(l,E,t0,tso,a,c,num_laps,tph,tsoph,ws,qs)
	ew,ev = eigen(H)
	ds = [ev[i:i+1,ev_ind]'*ev[i:i+1,ev_ind] for i in 1:2:length(ew)]
	Nd[ind,:] = ds[:]
#	plt.plot(ds .- 0.5*ind,label=L"$\lambda = $"*"$(tso)")
end
plt.imshow(Nd,aspect="auto")
#plt.legend()
plt.show()

end

