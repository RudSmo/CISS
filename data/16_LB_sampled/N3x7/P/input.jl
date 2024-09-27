p = expanduser("~/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_NEGF/noSpin/") 
include("$(p)hamiltonian.jl")
include("$(p)observables.jl")
include("$(p)structures.jl")
using .Hamiltonian
using .Observables
using .CustomStructs
using DelimitedFiles
using LinearAlgebra
using PyPlot

num_laps = 3
num_ions = 7
N = num_laps * num_ions
ws  = 0.2*0.002*ones(N)
T0  = 0.04/ws[1]  # in eV
t0  = -T0
eu  = -6*T0
ed  = -6*T0
E   = [eu,ed]
gPh = zeros(2)
a   = 1.0
c   = 1.0
tso = T0/40
tPh = 0.1*T0
tsoPh = 0.1*tso
model = FranssonHamiltonian
modelStr = "Fransson"

Xdat = readdlm("init_X_traj1.dat")

for i in 1:1

	qs = Xdat[i,:]
	HS = Hamiltonian.FranssonHamiltonian(N,E,t0,tso,a,c,num_laps,tPh,tsoPh,ws,qs)
	println(eigvals(HS))

	En = Array(2*eu+3*T0:1:-300)
	sys_par = SystemParams(N,HS,En)

	P = 0.5
	Gamma = T0/4
	gLu = Gamma*(1 + P)
	gLd = Gamma*(1 - P) 
	gLud = 0.0
	gLdu = gLud

	gRu = Gamma
	gRd = Gamma
	gRud = 0.0
	gRdu = gRud

	Temp = ws[1]/0.025
	muL = 0.0
	muR = 0.0
	tL  = Temp
	tR  = Temp

	lead_par = LeadParamsSF(muL,muR,gLu,gLd,gLud,gLdu,gRu,gRd,gRud,gRdu,tL,tR)


	#tfunc = LB_TransmissionE(sys_par,lead_par)
	#plt.plot(En,tfunc)
	#plt.show()
	Vn = Array(0.0:25:-2*eu+2*T0)


	filname = "fS$(num_ions)_nL$(num_laps)_tso$(tso)_P$(P)_currents_$(i).dat"
	savestp = 1
	LB_CurrentV(sys_par,lead_par,Vn,filname,savestp)
end
