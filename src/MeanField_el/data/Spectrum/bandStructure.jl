
include("../../hamiltonian.jl")

using .Hamiltonian
using LinearAlgebra
using FFTW
using PyPlot
# -------------- Parameters --------------------------
t0 = 1.0
eu = 0.0
ed = 0.0
E = [eu,ed]
a = 1.0
c = 1.0
num_laps = 20
ions = 4
l = Int(num_laps * ions)
tph = 0.0
tsoph = 0.0
ws = ones(l)*1.4*t0
qs = zeros(l)

#Tso = Array(0.0:0.01:20.0)
tso = 1.0


#for tso in Tso
H = Hamiltonian.FranssonHamiltonian(l,E,t0,tso,a,c,num_laps,tph,tsoph,ws,qs)
ew = eigvals(H)

plt.plot(ew)
plt.show()




