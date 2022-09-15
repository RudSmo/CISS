using ITensors
using LinearAlgebra

function HamiltonianMPO(sites,N,t0,Eu,Ed)
	os = OpSum()
	for i in 1:N
		os += Eu*("Cdagup",i,"Cup",i)
		os += Ed*("Cdagdn",i,"Cdn",i)
	end
	for i in 1:N-1
		os += -t0*("Cdagup",i,"Cup",i+1)
		os += -t0*("Cdagup",i+1,"Cup",i)
		os += -t0*("Cdagdn",i,"Cdn",i+1)
                os += -t0*("Cdagdn",i+1,"Cdn",i)
	end
	return MPO(os,sites)
end





L = 3
dt = 0.001
tend = 50
t0 = 1.0
Eu = 1.0
Ed = -1.0
indEE = 2
sites = siteind("Electron",N)
state = InitState(sites)
for i in 1:N
	if i == 1
		state.set(i,"UpDn")
	else
		state.set(i,"Emp")
	end
end
psi = MPS(state)
ampo = HamiltonianMPO(sites,N,t0,Eu,Ed)



