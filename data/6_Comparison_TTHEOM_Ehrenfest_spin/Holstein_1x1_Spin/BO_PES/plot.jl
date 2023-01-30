using LinearAlgebra
using PyPlot
src_path=expanduser("/Users/rudolfsmorka/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v2_weighting/")
include("$(src_path)hamiltonian.jl")
include("$(src_path)structures.jl")

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 2
rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

using .Hamiltonian
using .CustomStructures

num_laps = 1
num_ions = 1
N = num_laps * num_ions
t0  = 0.04
eu  = -6*t0
ed  = -6*t0
E   = [eu,ed]
ws  = 0.05*t0*ones(N)
gPh1 = [0.1,0.1]*t0
gPh2 = [0.1,-0.1]*t0
a   = 1.0
c   = 1.0
tso = 0.0*t0
tPh = 0.0
tsoPh = 0.0*tso
model = HolsteinHopping
modelStr = "Holstein"
#model = FranssonHamiltonian
#modelStr = "Fransson"
xnIn = zeros(N)
pnIn = zeros(N)

P = 0.5
tL = t0/4 * (ones(2) .+ P * [1,-1])
tR = t0/4 * ones(2)
function PhononHamiltonian(HamPar)
        H = zeros(2*HamPar.N,2*HamPar.N)
        for i in 1:HamPar.N
                H[2*i-1,2*i-1] = 0.5*HamPar.ws[i]^2 * HamPar.X[i]^2
                H[2*i,2*i] = 0.5*HamPar.ws[i]^2 * HamPar.X[i]^2
        end
        return H
end


fig = plt.figure(figsize=(9,6))

px = 0.1
py = 0.1
sx = 0.415
sy = 0.87
dx = 0.05

ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px+sx+dx,py,sx,sy])
axes = [ax1,ax2]
parameters = [gPh1, gPh2]
title_labels = [L"$g_\uparrow = g_\downarrow$", L"$g_\uparrow = -g_\downarrow$"]
for (ax,params,title) in zip(axes,parameters,title_labels)

	Xrng = -200:1:200
	nx = length(Array(Xrng))
	colors = ["blue","red"]
	for (col,state) in zip(colors,1:2*N)
		PES = zeros(nx)
		for (i,X) in enumerate(Xrng)
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
				    params,
				    ws,
				    [X],
				    pnIn,
				    tL,
				    tR)


			Hel = model(HamPar)
			Hph = PhononHamiltonian(HamPar)
			H = Hel .+ Hph
			ew = eigvals(H)
			PES[i] = ew[state]
		end
		ax.plot(Xrng,PES,color=col,linewidth=lw)
		ax.set_xlim(Xrng[1],Xrng[end])
		
		ax.annotate(title,(0.4,0.9),xycoords="axes fraction",fontsize=fs)
		ax.set_xlabel(L"$X$",labelpad=10)
		if ax == ax1
			ax.set_ylabel(L"$\varepsilon\,[\mathrm{eV}]$")
		end
	end
end
ax2.set_yticklabels([])
savefig("PES_Holstein_$(num_laps)x$(num_ions)_Spin.pdf")
