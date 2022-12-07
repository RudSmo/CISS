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
t0  = 1.0
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
En = Array(-12:0.001:2.0)

function TransmissionFunc(En,P,beta,mu,num_laps,
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
                            pnIn) 
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
	GammaL = zeros(2*N,2*N)*(1+0*1im)
	GammaR = GammaL
	GammaL[1,1] = tL[1]
	GammaL[2,2] = tL[2]
	GammaR[end-1,end-1] = tR[1]
	GammaR[end,end] = tR[2]
	Gamma = GammaL + GammaR

	Id = Matrix(I,2*N,2*N)

	transFunc = []
	for E in En
		Gr = inv(Id * E - H + 1im/2 * Gamma - Id * 1im * 1e-12)
		Ga = Gr'
		Trans = tr(Gr * GammaL * Ga * GammaR)
		push!(transFunc,Trans)
	end
	return transFunc
end

beta = 1/0.025
mu = 0.0
TfuncP = TransmissionFunc(En,0.5,beta,mu,num_laps,
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
                            pnIn)
TfuncN = TransmissionFunc(En,0.0,beta,mu,num_laps,
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
                            pnIn)

lw = 3
fig = plt.figure(figsize = (10,5))
plt.plot(En,TfuncP-TfuncN,color="k",linewidth=lw-1)
plt.fill_between(En,TfuncP-TfuncN,y2=0,color="grey",alpha=0.6)
println(sum(TfuncP-TfuncN))
plt.xlabel(L"$\varepsilon \, [t_0]$",labelpad=-1)
plt.ylabel(L"$T(\varepsilon,P=0.5)-T(\varepsilon,P=0)$")

plt.savefig("Transmission_function_Difference.pdf")

