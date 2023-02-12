using LinearAlgebra
using PyPlot
using SparseArrays
include("../../../src/exact_TDSE/hamiltonian.jl")
using .Hamiltonian
num_sites = 2
num_laps = 1
vibStates = 3
T0 = 0.04
t0 = T0
E =  0.0 #-6*t0
eps = [E,E]
tPh = 0.0 # * T0
λ0 = 0.0
λPh = 0.0 #  0.1 * T0 / 40
wPh = 0.05 * T0
r_helix = 1.0
c_helix = 1.0

H = Hamiltonian.Fransson(num_sites,num_laps,vibStates,eps,t0,tPh,λ0,λPh,wPh,r_helix,c_helix)
display(sparse(H))

