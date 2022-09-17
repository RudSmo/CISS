using ITensors
using LinearAlgebra
using TimeEvoMPS
using PyPlot


function VDot(a)
    sx = [0 1; 1 0.]
    sy = [0 -1im; 1im 0.]
    sz = [1 0; 0 -1.]
    return a[1]*sx + a[2]*sy + a[3]*sz
end
function absv(a)
    return sqrt(dot(a,a))
end
#-------------------------- Electronic Hamiltonian-------------------------------------
function UnitVector(a,c,j,num_laps,N)
    phij = 2*pi*(j-1)*num_laps/(N-1)
    cj = c*(j-1)/(N-1)
    return [a*cos(phij),a*sin(phij),cj]
end
function SOC_A_nn(N,tso,a,c,num_laps,k)
 
    rj = UnitVector(a,c,k,num_laps,N)

    rj1p = UnitVector(a,c,k+1,num_laps,N)
    rj2p = UnitVector(a,c,k+2,num_laps,N)
    rjrj1p = rj-rj1p
    rjrj2p = rj-rj2p

    dj1p = rjrj1p/absv(rjrj1p)
    dj2p = rjrj2p/absv(rjrj2p)

    vjp = cross(dj1p,dj2p)
    Vjp = 1im*tso*VDot(vjp)

    return Vjp
end


function HamiltonianMPO(sites,N,Eu,Ed,t0,tso,a,c,num_laps)
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
	for i in 1:N-2
		Vj = SOC_A_nn(N,tso,a,c,num_laps,i)		
		add!(ampo,Vj[1,1],"Cdagup",i,"Cup",i+2)
		add!(ampo,conj(Vj[1,1]),"Cdagup",i+2,"Cup",i)
		add!(ampo,Vj[1,2],"Cdagup",i,"Cdn",i+2)
		add!(ampo,conj(Vj[1,2]),"Cdagdn",i+2,"Cup",i)
		add!(ampo,Vj[2,1],"Cdagdn",i,"Cup",i+2)
                add!(ampo,conj(Vj[2,1]),"Cdagup",i+2,"Cdn",i)
		add!(ampo,Vj[2,2],"Cdagdn",i,"Cdn",i+2)
                add!(ampo,conj(Vj[2,2]),"Cdagdn",i+2,"Cdn",i)
	end	
	return MPO(ampo,sites)
end




num_sites = 3
num_laps = 1
N = Int(num_sites*num_laps)
dt = 0.01
tend = 5
t0 = 1.0
tso = 0.8
a = 1.0
c = 1.0
Eu = 1.0
Ed = -1.0
indEE = 2
maxdim = 100
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
H = HamiltonianMPO(sites,N,Eu,Ed,t0,tso,a,c,num_laps)
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




