using ITensors
using LinearAlgebra
using ITensorTDVP
using PyPlot
include("solvers.jl")

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
function ode_solver(H, time_step, psi; kwargs...)
  return ode_solver(
    -im * H,
    time_step,
    psi;
    solver_alg=ode_alg,
    ode_kwargs...,
    kwargs...,
  )
end
function krylov_solver(H, time_step, psi; kwargs...)
  return krylov_solver(
    -im * H, time_step, psi; krylov_kwargs..., kwargs...
  )
end


H = HamiltonianMPO(sites,N,t0,Eu,Ed)


start_linkdim = 4
outputlevel = 1
time_step = 0.1
time_stop = 1.0
maxdim = 100
cutoff = 1e-8
tol = 1e-15
krylov_kwargs = (; tol=tol, eager=true)
@show krylov_kwargs

#ode_alg = Tsit5()
#ode_kwargs = (; reltol=tol, abstol=tol)


nsite = 2
psi_t = tdvp(krylov_solver, H, time_stop, psi; time_step, maxdim, cutoff, nsite, outputlevel)
print(psi_t)


#cb = LocalMeasurementCallback(["Nup","Ndn"], sites, 0.2)
#tdvp!(psi,H,dt,tend,maxdim=maxdim,callback=cb)
#ts = measurement_ts(cb)


#fil = open("pops.dat","w")
#
#Nu = [getindex.(measurements(cb)["Nup"],j) for j in 1:N]
#Nd = [getindex.(measurements(cb)["Ndn"],j) for j in 1:N]
#
#for i in 1:length(ts)
#	
#	write(fil,"$(ts[i]) ")
#	
#	for j in 1:N
#		Nu = 
#		write(fil,"$(Nu[j][i]) $(Nd[j][i]) ")
#	end
#	write(fil,"\n")
#end
#close(fil)
#



