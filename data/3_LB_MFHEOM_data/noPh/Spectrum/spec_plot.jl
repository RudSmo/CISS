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
num_ions = 5
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

P = 0.0
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
ew = eigvals(H)


X = Array(-1:0.01:1)
for i in 1:length(ew)
	plt.plot(X,ew[i]*ones(length(X)),color="k",linewidth=1)
end

plt.plot([-4,-2],[0,0],color="k")
plt.plot([2,4],[0,0],color="k")
plt.plot([-2,-2],[-9,0],color="k")
plt.plot([2,2],[-9,0],color="k")
plt.ylabel(L"$\varepsilon \, [t_0]$")
plt.tick_params(direction="in")
plt.xticks([])
plt.annotate(L"$\mu_L$",(-0.1,0.89),xycoords="axes fraction",fontsize=fs)
plt.annotate(L"$\mu_R$",(1.02,0.89),xycoords="axes fraction",fontsize=fs)
plt.xlim(-4,4)
plt.ylim(-8.5,1)

plt.savefig("Spectrum_leads.pdf")



