using LinearAlgebra
using PyPlot
using3D()
src_path=expanduser("/Users/rudolfsmorka/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v3_typed_PVar/")
include("$(src_path)hamiltonian.jl")
include("$(src_path)structures.jl")

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 1.5
rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)


num_laps = 1
num_ions = 2
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
tso = 0.0*T0/40
tPh = 1.0*T0*ones(2)
tsoPh = 0.0*tso
#model = HolsteinHopping
#modelStr = "Holstein"
model = FranssonHamiltonian
modelStr = "Fransson"
xnIn = zeros(N)
pnIn = zeros(N)

P = 0.5
tL = t0/4 * (ones(2) .+ P * [1,-1])
tR = t0/4 * ones(2)
function PhononHamiltonian(HamPar)
        H = zeros(2*HamPar.N,2*HamPar.N)
        for i in 1:HamPar.N
                H[2*i-1,2*i-1] = 0.5* HamPar.X[i]^2
                H[2*i,2*i] = 0.5*HamPar.X[i]^2
        end
        return H
end


fig = plt.figure(figsize=(10,7))

px = 0.09
py = 0.1
sx = 0.9
sy = 0.9
dx = 0.01
dy = 0.015


Xrng = -20:0.1:20
nx = length(Array(Xrng))
colors = ["blue","red","green","gray","black","orange","purple","turquoise"]
X2s  = [-15,-12,-10,-7,-3]
X3s  = [0,5,10,12,15]

ax = plt.axes([px,py,sx,sy])
PES = []
for X1 in Xrng
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
			    [X1,0.0],
			    pnIn,
			    tL,
			    tR)


	Hel = model(HamPar)
	Hph = PhononHamiltonian(HamPar)
	H = Hel .+ Hph
	ew = eigvals(H)
	push!(PES,ew)
end
PES = [[PES[a][b] for a in 1:length(PES)] for b in 1:length(PES[1])]
for (pes,col) in zip(PES,colors)
	ax.plot(Xrng,pes,linewidth=lw,color=col)
end

plt.savefig("PES_$(num_laps)x$(num_ions)_SOC$(tso)_tPh$(tPh[1])_SOCPh$(tsoPh)_multiState_cuts_v2.pdf")
