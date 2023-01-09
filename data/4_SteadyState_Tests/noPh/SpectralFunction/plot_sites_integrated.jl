using LinearAlgebra
using PyPlot
using DelimitedFiles
src_path=expanduser("/Users/rudolfsmorka/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v2_weighting/")
include("$(src_path)hamiltonian.jl")
include("$(src_path)structures.jl")

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

using .Hamiltonian
using .CustomStructures

num_laps = 2
num_ions = 4
N = num_laps * num_ions
t0  = 0.04
eu  = -6*t0
ed  = -6*t0
E   = [eu,ed]
ws  = 0.05*ones(N)
gPh = [0.0,0.0]
a   = 1.0
c   = 1.0
tso = t0/40
tPh = 0.0
tsoPh = 0.0
#model = HolsteinHopping
#modelStr = "Holstein"
model = FranssonHamiltonian
modelStr = "Fransson"
xnIn = zeros(N)
pnIn = zeros(N)

P = 0.5
tL = t0/4 * (ones(2) .+ P * [1,-1])
tR = t0/4 * ones(2)
beta = 1/0.025
mu = 0.0
HamPar = HamiltonianParams( num_laps,
                            num_ions,
                            N,
                            model,
                            modelStr,
                            E,
                            a,
                            c,
                            t0,
                            tPh,
                            tso,
                            tsoPh,
                            gPh,
                            ws,
                            xnIn,
                            pnIn,
                            tL,
                            tR)


H = model(HamPar)

Gamma = zeros(2*N,2*N)*(1+0*1im)
Gamma[1,1] = tL[1]
Gamma[2,2] = tL[2]
Gamma[end-1,end-1] = tR[1]
Gamma[end,end] = tR[2]

Id = Matrix(I,2*N,2*N)
En = Array(-12*t0:0.0005*t0:2.0*t0)
Nu = []
Nd = []

function Fermi(E,beta,mu)
	return 1/(exp(beta*(E-mu))+1)
end

for site in 1:N
	s = 2*site - 1

	specFuncS = []
	specFuncTu = []
	specFuncTd = []
	for E in En
		Gr = inv(Id * E - H + 1im/2 * Gamma + Id * 1im * 1e-12)
		Ga = Gr'

		fe = Fermi(E,beta,mu)
		Lambda = fe * Gr * Gamma * Ga
		push!(specFuncS,Lambda[s,s]-Lambda[s+1,s+1])
		push!(specFuncTu,Lambda[s,s])
		push!(specFuncTd,Lambda[s+1,s+1])
	end

	function Integrate(f,dx)
		return 1/(2*pi) * sum([0.5*(f[i+1]+f[i])*dx for i in 1:length(f)-1])
	end

	push!(Nu,Integrate(specFuncTu,En[2]-En[1]))
	push!(Nd,Integrate(specFuncTd,En[2]-En[1]))
end

lw = 3
plt.plot(Array(1:N),Nu,color="red",linewidth=lw,label=L"$N_\uparrow$"*" analytical NEGF")
plt.plot(Array(1:N),Nd,color="blue",linewidth=lw,label=L"$N_\downarrow$"*" analytical NEGF")

dat = readdlm("../V0.0/L8_tPh0.0_SOC0.001_SOCPh0.0_V0.0_w00.002_P0.5_X0_-1.0_dt0.05_1_electrons.dat")
println(dat[end-100,1])
plt.plot(Array(1:N),dat[end-100,2:2:end],color="red",linewidth=lw,linestyle="dashed",label=L"$N_\uparrow$"*" MF-HEOM")
plt.plot(Array(1:N),dat[end-100,3:2:end],color="blue",linewidth=lw,linestyle="dashed",label=L"$N_\downarrow$"*" MF-HEOM")
plt.legend(loc="best")
plt.xlabel("sites")
plt.ylabel("population")
plt.savefig("Comparison_populations_V0.pdf")

