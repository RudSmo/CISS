src_path=expanduser("~/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v3_typed_PVar/")
include("$(src_path)hamiltonian.jl")
include("$(src_path)structures.jl")

src=expanduser("~/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v3_typed_swarm_PVar/")
include("$(src)classical_phonons.jl")

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
tPh = 0.0*T0*ones(2)
tsoPh = 0.1*tso
model = FranssonHamiltonian
modelStr = "Fransson"
# Simulation Parameters
dt = 0.0002
tend = 50
numstep = Int(round(tend/abs(dt)))
writeMode = "minimal"


popElIn = zeros(2*N)
popElIn[1] = 0.0
popPhIn = zeros(N)

xDat = readdlm("phonon.dat")
#xnIn,pnIn = xDat[end,2:2:end],zeros(N)
#xDat[end,:], zeros(N)
num_avg = 100000
beta = ws[1]/0.025
xnIn,pnIn = SampleWigner_Multi(beta,ws,N,1)
Xsmpl,Psmpl = SampleWigner_Multi(beta,ws,N,num_avg)
println(size(Xsmpl))
## Leads
P = 0.0
tl = T0/4
tL = tl * [1,0]#(ones(2) .+ P * [1,-1])
tR = tl * [1,0]#ones(2)

v = 75
V = v * 0.2 * T0
nPoles = 15

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

function Dos(E_i,E,b=1e-4)
    D = [1/(sqrt(2*pi)*b)*sum([exp(-(e-i)^2/(2*b^2)) for i in E_i]) for e in E]
    return D
end

tso_rng = [0.0,0.0001,0.0005,0.001,0.005,0.01] ./ ws[1]
#tso_rng = 0.0:0.1*tso:50*tso 

Erng = -0.48:0.0005:0
DOS = zeros(length(Erng),length(tso_rng))
for (i,tsos) in enumerate(tso_rng)
	HamPar.SOCPh = tsos
	for j in 1:num_avg
		HamPar.X = Xsmpl[j,:]
		H = model(HamPar)
		ew = eigvals(H) * ws[1]
		dos = Dos(ew,Erng)
		DOS[:,i] .+= dos
	end
	DOS[:,i] = DOS[:,i] ./ num_avg
	writedlm("L$(num_laps)x$(num_ions)_tPh_$(tPh)_tsoPh_var_avg_dos.dat",hcat(Erng,DOS))
end
image = plt.imshow(DOS[end:-1:1,:] , aspect="auto",  extent = [tso_rng[1],tso_rng[end],Erng[1]*100,Erng[end]*100])
plt.xlabel(L"$\lambda_1 \, [\mathrm{meV}]$")
plt.ylabel(L"$E\,[\mathrm{meV}]$")
cbar = plt.colorbar(image,aspect=90,pad=.001)
cbar.set_label("DOS")

plt.savefig("L$(num_laps)x$(num_ions)_tPh_$(tPh)_tsoPhVar_avg.pdf")

