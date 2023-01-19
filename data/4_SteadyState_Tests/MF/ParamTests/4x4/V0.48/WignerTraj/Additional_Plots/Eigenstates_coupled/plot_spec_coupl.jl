using LinearAlgebra
using PyPlot
src_path=expanduser("/Users/rudolfsmorka/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v2_weighting/")
include("$(src_path)hamiltonian.jl")
include("$(src_path)structures.jl")

rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 15
rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

using .Hamiltonian
using .CustomStructures

num_laps = 4
num_ions = 4
N = num_laps * num_ions
t0  = 0.04
eu  = -6*t0
ed  = -6*t0
E   = [eu,ed]
ws  = 0.05*t0*ones(N)
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
V = 0.48
Gamma = zeros(2*N,2*N)
Gamma[1,1] = tL[1]
Gamma[2,2] = tL[2]
Gamma[end-1,end-1] = tR[1]
Gamma[end,end] = tR[2]
Heff = H - 1im*Gamma


ew = eigvals(H)
ewEff = eigvals(Heff)

X = Array(-1:0.01:1)
for i in 1:length(ew)
	plt.plot(X,ew[i]*ones(length(X)),color="k",linewidth=1)
	plt.plot(X,ewEff[i]*ones(length(X)), color="red",linewidth=1)
end

plt.plot([-4,-2],[V/2,V/2],color="k")
plt.plot([2,4],[-V/2,-V/2],color="k")
plt.plot([-2,-2],[-9,V/2],color="k")
plt.plot([2,2],[-9,-V/2],color="k")
plt.ylabel(L"$\varepsilon \, [\mathrm{eV}]$")
plt.tick_params(direction="in")
plt.xticks([])
plt.annotate(L"$\mu_L$",(-0.05,0.92),xycoords="axes fraction",fontsize=fs)
plt.annotate(L"$\mu_R$",(1.013,0.135),xycoords="axes fraction",fontsize=fs)
plt.xlim(-4,4)
plt.ylim(-8.5*t0,7*t0)

plt.savefig("Spectrum_leads_$(N)_tso$(tso)_$(V).pdf")



