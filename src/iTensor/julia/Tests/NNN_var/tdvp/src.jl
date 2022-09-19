using ITensors
using LinearAlgebra
using TimeEvoMPS
using PyPlot


function HamiltonianMPO(sites,N,tns)
	ampo = OpSum()
	for j in 1:length(tns)
		for i in 1:N-j
			ampo += (-tns[j],"Cdagup",i,"Cup",i+j)
			ampo += (-tns[j],"Cdagdn",i,"Cdn",i+j)
        	        ampo += (-tns[j],"Cdagup",i+j,"Cup",i)
			ampo += (-tns[j],"Cdagdn",i+j,"Cdn",i)
		end
	end
	return MPO(ampo,sites)
end


N = 4
dt = 0.01
tend = 50
Eon = 0.0
t1 = 1.0
t2 = 1.0
tns = [t1,t2]

sites = siteinds("Electron",N)

state = []
for i in 1:N
	if i == 1
		push!(state,"UpDn")
	else
		push!(state,"0")
	end
end

psi = productMPS(sites,state)
H = HamiltonianMPO(sites,N,tns)

cb = LocalMeasurementCallback(["Nup","Ndn"], sites, 0.2)
tdvp!(psi,H,dt,tend,maxdim=maxdim,callback=cb)
ts = measurement_ts(cb)
fil = open("L$(N)_n$(length(tns))_dt$(dt).dat","w")

nup = [getindex.(measurements(cb)["Nup"],j) for j in 1:N]
ndn = [getindex.(measurements(cb)["Ndn"],j) for j in 1:N]

for i in 1:length(ts)
	
	write(fil,"$(ts[i]) ")
	
	for j in 1:N
		write(fil,"$(nup[j][i]) $(ndn[j][i])")
	end
	write(fil,"\n")
end
close(fil)




