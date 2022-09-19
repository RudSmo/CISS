using ITensors
using LinearAlgebra
using TimeEvoMPS
using PyPlot

ITensors.op(on::OpName"F", st::SiteType"Qudit", ds::Int...) = op(OpName"Id"(), st, ds...)

function HamiltonianMPO(sites,N,Eon,t0,gPh)
	ampo = OpSum()
	for i in 1:2:N-1
		ampo += (Eon,"Cdag",i,"C",i)
	end
	for i in 1:2:N-3
		ampo += (-t0,"Cdag",i,"C",i+2)
                ampo += (-t0,"Cdag",i+2,"C",i)
	end
	for i in 1:2:N-1
		ampo += (gPh,"n",i,"adag",i+1)
                ampo += (gPh,"n",i,"a",i+1)
	end
	println(ampo)
	return MPO(ampo,sites)
end



fS = 3
N = 2*fS
dt = 0.01
tend = 5
Eon = 0.0
t0 = 1.0
gPh = 0.0

sites = siteinds(n->isodd(n) ? "Fermion" : "Boson",N)

state = []
for i in 1:2:N
	if i == 1
		push!(state,"1")
		push!(state,"0")
	else
		push!(state,"0")
		push!(state,"0")
	end
end

psi = productMPS(sites,state)
println(psi)
H = HamiltonianMPO(sites,N,Eon,t0,gPh)

cb = LocalMeasurementCallback(["N"], sites[1:2:end-1], 0.2)
tdvp!(psi,H,dt,tend,maxdim=maxdim,callback=cb)
ts = measurement_ts(cb)
fil = open("pops.dat","w")

nm = [getindex.(measurements(cb)["N"],j) for j in 1:2:N-1]

for i in 1:length(ts)
	
	write(fil,"$(ts[i]) ")
	
	for j in 1:fS
		write(fil,"$(nm[j][i]) ")
	end
	write(fil,"\n")
end
close(fil)
#



