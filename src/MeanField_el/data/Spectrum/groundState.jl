include("../../hamiltonian.jl")
include("../../electrons.jl")
using .Hamiltonian
using .Electrons
using PyPlot

t0 = 1.0
eu = 0.0
ed = 0.0
E = [eu,ed]
a = 1.0
c = 1.0
tso = t0/4
num_laps = 3
ions = 4
l = Int(num_laps * ions)
tph = 0.0
tsoph = 0.0
ws = ones(l)*1.4*t0
qs = zeros(l)

H = Hamiltonian.FranssonHamiltonian(l,E,t0,tso,a,c,num_laps,tph,tsoph,ws,qs)

rho = DensityFEx(l,1,H,l,0)
plt.plot([rho[i,i] for i in 1:2:2*l-1])
plt.show()


