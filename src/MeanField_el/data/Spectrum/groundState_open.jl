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
tso = 4*t0
num_laps = 10
ions = 4
l = Int(num_laps * ions)
tph = 0.0
tsoph = 0.0
ws = ones(l)*1.4*t0
qs = zeros(l)

H = Hamiltonian.FranssonHamiltonian(l,E,t0,tso,a,c,num_laps,tph,tsoph,ws,qs)
nP = 50

beta = 1/0.025
kappa,xi= PadeFermi(nP)
chi = 1im .* xi ./beta
chiLP = [chi,chi] 
tL = ones(2)
tR = tL
Gamma1,Gamma2   = GammaMatrices1d(l,tL,tR)
GammaT          = Gamma1+Gamma2

rho = DensityFExOpen(l,1,H,kappa,chiLP,GammaT,beta)
#plt.plot([rho[i,i]-rho[i+1,i+1] for i in 1:2:2*l-1])
plt.plot([rho[i,i] for i in 1:2:2*l],color="red")
plt.plot([rho[i,i] for i in 2:2:2*l],color="blue")
plt.show()


