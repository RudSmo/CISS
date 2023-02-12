p = expanduser("~/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_NEGF/spinResolved/")
include("$(p)hamiltonian.jl")
include("$(p)observables.jl")
include("$(p)structures.jl")

using PyPlot
using LinearAlgebra
using .Hamiltonian
using .Observables
using .CustomStructs

t0 = 0.04
tu = t0
td = 0.75*t0
tud = 0.0
tdu = 0.0
TD = [tu,td]
TO = [tud,tdu]
eu = -6*t0
ed = eu
E = [eu,ed]
a = 1.0
c = 1.0
tso = t0/4
num_laps = 2
ions = 4
N = Int(num_laps * ions)
tph = 0.0
tsoph = 0.0
ws = ones(N)*1.4*t0
qs = zeros(N)

#HS = Hamiltonian.FranssonHamiltonian(N,E,t0,tso,a,c,num_laps,tph,tsoph,ws,qs)
HS = Hamiltonian.SpinHopping(N,E,TD,TO)

En = Array(-10*t0:0.01*t0:t0)
sys_par = SystemParams(N,HS,En)

P = -0.5
gLu = t0/4*(1+P)
gLd = t0/4*(1-P)
gLud = 0.0#1im*20.0*1e-1
gLdu = gLud

gRu = t0/4
gRd = t0/4
gRud = 0.0#1im*20.0*1e-1
gRdu = gRud

muLu = 0.0
muLd = 0.0
muRu = 0.0
muRd = 0.0
tLu  = 0.001#0.025
tLd  = tLu
tRu  = tLu
tRd  = tLu

lead_par = LeadParamsSF(muLu,muLd,muRu,muRd,gLu,gLd,gLud,gLdu,gRu,gRd,gRud,gRdu,tLu,tLd,tRu,tRd)


Vnu = Array(5.0*t0:0.1*t0:24*t0)
Vnd = Array(5.0*t0:0.1*t0:24*t0)

#Tu,Td,Tud,Tdu = LB_TransmissionE(sys_par,lead_par)
#plt.plot(En,Tu,label="u")
#plt.plot(En,Td,label="d")
#plt.plot(En,Tud,label="ud")
#plt.plot(En,Tdu,label="du")
#plt.legend()
#plt.show()



filname = "fS$(ions)_nL$(num_laps)_tu$(tu)_$(td)_P$(P)_currents_spin.dat"
savestp = 1
LB_CurrentV(sys_par,lead_par,Vnu,Vnd,filname,savestp)
#end

