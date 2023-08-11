include("hamiltonian.jl")
using PyPlot
using LinearAlgebra
using SparseArrays
using Combinatorics 
using StatsBase


function LevelSpacings(sites,laps,vibSt,epsilon,t,t1,tso,tsoPh,w,r_hel,c_hel,PBC)

	H = Fransson(sites,laps,vibSt,epsilon,t,t1,tso,tsoPh,w,r_hel,c_hel,PBC)
	ew = eigvals(H)

	Si = [ew[i+1]-ew[i] for i in 1:length(ew)-1]
	Ri = [min(Si[i],Si[i-1])/max(Si[i],Si[i-1]) for i in 2:length(Si)]
	return Ri
end
	#Sm = 1/length(Si) * sum(Si)
	#Si = Si ./ Sm
	#print(1/length(Ri) * sum(Ri))

	#println(Ri)



num_sites = 3
num_laps = 1
N = num_sites * num_laps
vibStates = 1
T0 = 0.04
t0 = -T0
E =  -6*T0
eps = [E,E]
tPh = 0.1*t0
λ0 = T0/40
λPh = 0.1 * λ0
wPh = 0.002
r_helix = 1.0
c_helix = 1.0

P = 0.5
tl = T0/4
tL = tl * (ones(2) .+ P * [1,-1])
tR = tl * ones(2)
Ri1 = LevelSpacings(num_sites,num_laps,vibStates,eps,t0,tPh,λ0,λPh,wPh,r_helix,c_helix,true)
Ri2 = LevelSpacings(num_sites,num_laps,vibStates,eps,t0,tPh,λ0,λPh,wPh,r_helix,c_helix,true)
plt.hist(Ri1,bins=30,color="k")
plt.hist(Ri2,bins=30,color="red",alpha=0.5)
plt.show()
