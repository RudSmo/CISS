src_path=expanduser("~/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v3_typed_PVar/")
include("$(src_path)hamiltonian.jl")
include("$(src_path)structures.jl")

using PyPlot
using LinearAlgebra
using DelimitedFiles
using HDF5

# Hamiltonian Parameters
num_laps = 2
num_ions = 4
N = num_laps * num_ions
ws  = 0.002*ones(N)
Hop = 0.04
T0  = Hop/ws[1]  # in eV
t0  = -T0
eu  = -6*T0
ed  = -6*T0
E   = [eu,ed]
gPh = zeros(2)
a   = 1.0
c   = 1.0
tso = T0/40
tPh = 0.1*T0*ones(2)
tsoPh = 0.1*tso
model = FranssonHamiltonian
modelStr = "Fransson"
# Simulation Parameters
dt = 0.0002
tend = 3
numstep = Int(round(tend/abs(dt)))
writeMode = "minimal"


popElIn = zeros(2*N)
popElIn[1] = 0.0
popPhIn = zeros(N)

xDat = readdlm("phonon.dat")
xnIn,pnIn = xDat[1,2:2:end],zeros(N)
#xnIn,pnIn = zeros(N),zeros(N)
#xDat[end,:], zeros(N)
## Leads
P = 0.0
tl = T0/4
tL = tl * [1,0]#(ones(2) .+ P * [1,-1])
tR = tl * [1,0]#ones(2)

v = 75
V = v * 0.2 * T0
nPoles = 15
beta = ws[1]/0.025

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

fac_t = 1/(0.002 * 1e3 * 4.135668)
time = xDat[2:end,1]
state = N 
ES_pop = zeros(length(time),N)
px = 0.085
py = 0.13
sx = 0.73
sy = 0.83

fig = plt.figure(figsize=(6,4))
ax = plt.axes([px,py,sx,sy])
cax = plt.axes([px+sx+0.02,py,0.02,sy])
cp = plt.cm.viridis
cmap = plt.cm.ScalarMappable(cmap=cp)

for i in 2:length(time)+1
	HamPar.X = xDat[i,2:2:end]
	H = model(HamPar)
	ev = eigvecs(H)[state,:]
	for j in 1:N-1
		k = 2*j-1
		pop = real(ev[k:k+1]' * [1 0; 0 1.] * ev[k+2:k+3])
		ES_pop[i-1,j] = pop
	end
end
vmin = minimum(ES_pop)
vmax = maximum(ES_pop)
image = ax.imshow(ES_pop,aspect="auto",extent=[1,N-1,time[end],time[1]],vmin=vmin,vmax=vmax,cmap="viridis")
cbar = fig.colorbar(image,ax=ax,cax = cax,aspect=90,pad=.001)
cbar.set_label(L"$\langle \psi_{j}^{(k)} \vert \psi_{j+1}^{(k)}\rangle $")

ax.set_xlabel("site "*L"$j$")
ax.set_ylabel(L"$t\,[\mathrm{ps}]$")

plt.savefig("L$(num_laps)x$(num_ions)_tso_$(tso)_tsoPh_$(tsoPh)_XVar_spacetime_hop_es_$(state).pdf")

