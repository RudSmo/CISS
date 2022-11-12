include("hamiltonian.jl")
include("observables.jl")
include("structures.jl")
using PyPlot
using LinearAlgebra
using .Hamiltonian
using .Observables
using .CustomStructs

t0 = 1.0
eu = 0.0
ed = 0.0
E = [eu,ed]
a = 1.0
c = 1.0
tso = t0/4
num_laps = 3
ions = 4
N = Int(num_laps * ions)
tph = 0.0
tsoph = 0.0
ws = ones(N)*1.4*t0
qs = zeros(N)

HS = Hamiltonian.FranssonHamiltonian(N,E,t0,tso,a,c,num_laps,tph,tsoph,ws,qs)


En = Array(-4.0:0.01:4.0)
sys_par = SystemParams(N,HS,En)

gLu = 0.25
gLd = 0.25
gLud = 0.0#1im*20.0*1e-1
gLdu = gLud

gRu = gLu
gRd = gLd
gRud = 0.0#1im*20.0*1e-1
gRdu = gRud

muLu = 0.0
muLd = 0.0
muRu = 0.0
muRd = 0.0
tLu  = 0.0001
tLd  = 0.0001
tRu  = 0.0001
tRd  = 0.0001

lead_par = LeadParamsSF(muLu,muLd,muRu,muRd,gLu,gLd,gLud,gLdu,gRu,gRd,gRud,gRdu,tLu,tLd,tRu,tRd)


Vnu = Array(-5.0:0.03:5.0)
Vnd = Array(-5.0:0.03:5.0)

#Tu,Td,Tud,Tdu = LB_TransmissionE(sys_par,lead_par)
#plt.plot(En,Tu,label="u")
#plt.plot(En,Td,label="d")
#plt.plot(En,Tud,label="ud")
#plt.plot(En,Tdu,label="du")
#plt.legend()
#plt.show()


path = "data/Fransson_L$(N)"
mkpath(path)

filname = "$(path)/N$(N)_currents_spin.dat"
savestp = 1
LB_CurrentV(sys_par,lead_par,Vnu,Vnd,filname,savestp)
#end

