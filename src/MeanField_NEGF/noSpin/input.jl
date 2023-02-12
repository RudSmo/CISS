include("hamiltonian.jl")
include("observables.jl")
include("structures.jl")
using PyPlot
using LinearAlgebra
using .Hamiltonian
using .Observables
using .CustomStructs

t0 = 0.04
tu = t0
td = 0.75*t0
T0 = [tu,td]
eu = -6*t0
ed = -6*t0
E = [eu,ed]
a = 1.0
c = 1.0
tso = t0/40
num_laps = 2
ions = 4
N = Int(num_laps * ions)
tph = 0.1*t0
tsoph = 0.1*tso
ws = ones(N)*0.01*t0

tL  = 0.025
tR  = 0.025
beta = 1/tL
qs = [SampleWigner(beta,ws[i],"v2")[1] for i in 1:N] 

#HS = Hamiltonian.FranssonHamiltonian(N,E,t0,tso,a,c,num_laps,tph,tsoph,ws,qs)
HS = Hamiltonian.SpinHopping(N,E,T0)

En = Array(-10*t0:0.01*t0:2*t0)
sys_par = SystemParams(N,HS,En)

P = 0.5
gLu = 0.25*(1+P)
gLd = 0.25*(1-P)
gLud = 0.0
gLdu = 0.0

gRu = gLu
gRd = gLd
gRud = 0.0#1im*20.0*1e-1
gRdu = gRud

muL = 0.0
muR = 0.0


lead_par = LeadParamsSF(muL,muR,gLu,gLd,gLud,gLdu,gRu,gRd,gRud,gRdu,tL,tR)


Vnu = Array(0.0*t0:0.1*t0:25*t0)
Vnd = Array(0.0*t0:0.1*t0:25*t0)

T = LB_TransmissionE(sys_par,lead_par)
D = DOS(sys_par.En,sys_par.HS,0.02)

path = "SpinHopping_$(num_laps)x$(ions)"
if !isdir(path)
	mkpath(path)
end

filname = "$(path)/N$(N)_ions$(ions)_tu$(tu)_td$(td)_P$(P)"

filnameT = open("$(filname)_transmission.dat","w")

filnameD = open("$(filname)_dos.dat","w")
for i in 1:length(T)
	write(filnameT,"$(En[i]) $(T[i]) \n")
	write(filnameD,"$(En[i]) $(D[i]) \n")
end

close(filnameT)
close(filnameD)

filname = "$(filname)_currents.dat"
savestp = 1
LB_CurrentV(sys_par,lead_par,Vnu,Vnd,filname,savestp)
