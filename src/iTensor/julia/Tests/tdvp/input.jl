using ITensors
using ITensorTDVP
using KrylovKit
using LinearAlgebra
using Observers


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
time_stop = 50.0
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
print(H)
#outputlevel = 1
#nsite = 1
#start_linkdim = 4
#maxdim = 100
#cutoff = 1e-8
#tol = 1e-15
#krylov_kwargs = (; tol=tol, eager=true)
#
#function step(; sweep, bond, half_sweep)
#    if bond == 1 && half_sweep == 2
#      return sweep
#    end
#    return nothing
#  end
#function current_time(; current_time, bond, half_sweep)
#    if bond == 1 && half_sweep == 2
#      return current_time
#    end
#    return nothing
#  end
#function measure_nu(; psi, bond, half_sweep)
#    if bond == 1 && half_sweep == 2
#        Nup = []
#        for j in 1:N
#            push!(Nup,expect(psi, "Nup"; sites=j))
#        end
#        return Nup
#    end
#    return nothing
#end
#
#function measure_nd(; psi, bond, half_sweep)
#    if bond == 1 && half_sweep == 2
#        Ndn = []
#        for j in 1:N
#            push!(Ndn,expect(psi, "Ndn"; sites=j))
#        end
#        return Ndn
#    end
#    return nothing
#end
#
#function return_state(; psi, bond, half_sweep)
#    if bond == 1 && half_sweep == 2
#      return psi
#    end
#    return nothing
#  end
#obs = Observer("steps" => step, "times" => current_time, "psis" => return_state, "Nup" => measure_nu, "Ndn" => measure_nd)
#psi_t = tdvp(
#    H,
#    time_stop*1im,
#    psi;
#    time_step=1im*dt,
#    cutoff,
#    normalize=true,
#    maxdim=30,
#    (observer!)=obs,
#    outputlevel=1)    
#res = results(obs)
#steps = res["steps"]
#times = res["times"]
#psis = res["psis"]
#Nup = res["Nup"]
#Ndn = res["Ndn"]
#
#fil = open("TDVP_L$(N)_dt$(dt).dat","w")
#for i in 1:length(Nup)
#    write(fil,"$(imag(times[i])) ")
#    for j in 1:length(Nup[i])
#        write(fil,"$(Nup[i][j]) $(Ndn[i][j]) ")
#    end
#    write(fil,"\n")
#end
##plt.plot(times,Nup[:,1])
##print(psi_t)
