src_path=expanduser("~/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v3_typed_PVar/")
include("$(src_path)hamiltonian.jl")
include("$(src_path)structures.jl")

using PyPlot
using LinearAlgebra



# Hamiltonian Parameters
num_laps = 4
num_ions = 4
N = num_laps * num_ions
ws  = 0.002*ones(N)
Hop = 0.04
T0  = Hop/ws[1]  # in eV
tso = T0/40
a   = 1.0
c   = 1.0
Hso = SOC_A_nn(N,tso,a,c,num_laps)

ew = eigvals(Hso)

plt.plot(ew)
plt.show()
