using ITensors
using LinearAlgebra
using TimeEvoMPS
using PyPlot

function HamiltonianMPO(sites,N,t0,Eu,Ed)
	ampo = AutoMPO()
	for i in 1:N
		add!(ampo,Eu,"Cdagup",i,"Cup",i)
		add!(ampo,Ed,"Cdagdn",i,"Cdn",i)
	end
	for i in 1:N-1
		add!(ampo,-t0,"Cdagup",i,"Cup",i+1)
		add!(ampo,-t0,"Cdagup",i+1,"Cup",i)
		add!(ampo,-t0,"Cdagdn",i,"Cdn",i+1)
                add!(ampo,-t0,"Cdagdn",i+1,"Cdn",i)
	end
	return MPO(ampo,sites)
end





N = 3
dt = 0.01
tend = 5
t0 = 1.0
Eu = 1.0
Ed = -1.0
indEE = 2
maxdim = 50
sites = siteinds("Electron",N)
state = []
for i in 1:N
	if i == 1
		push!(state,"UpDn")
	else
		push!(state,"Emp")
	end
end
psi = productMPS(sites,state)
H = HamiltonianMPO(sites,N,t0,Eu,Ed)
cb = LocalMeasurementCallback(["Nup","Ndn"], sites, 0.2)
tdvp!(psi,H,dt,tend,maxdim=maxdim,callback=cb)
ts = measurement_ts(cb)
fil = open("pops.dat","w")

Nu = [getindex.(measurements(cb)["Nup"],j) for j in 1:N]
Nd = [getindex.(measurements(cb)["Ndn"],j) for j in 1:N]

for i in 1:length(ts)
	
	write(fil,"$(ts[i]) ")
	
	for j in 1:N
		write(fil,"$(Nu[j][i]) $(Nd[j][i]) ")
	end
	write(fil,"\n")
end
close(fil)




