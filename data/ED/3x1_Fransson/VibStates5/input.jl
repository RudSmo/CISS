module InputFile
src_path=expanduser("~/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/exact_TDSE/")
include("$(src_path)rk4.jl")
include("$(src_path)hamiltonian.jl")
using .RK4
using .Hamiltonian
using LinearAlgebra
using PyPlot

el_state = [1,1,0,0,0,0]
ph_state = [0,0,0]
vibStates = 5
num_sites = Int(length(el_state)/2)
num_laps = 1
sites_lap = num_sites
t0 = 0.04
eps = -0.24 * ones(2)
tPh = 0.004
lambda0 = 0.001
lambdaPh = 0.0001
wPh = 0.0004
r_hel = 1.0
c_hel = 1.0
psi0 = InitialState(el_state,ph_state,vibStates,num_sites)

H = Hamiltonian.Fransson(num_sites,num_laps,vibStates,eps,t0,tPh,lambda0,lambdaPh,wPh,r_hel,c_hel)
ew = eigvals(H)
for e in ew
	plt.plot(e*ones(length(0:1)),0:1,color="k")
end
plt.xlabel(L"$E\,[eV]$")
plt.savefig("L$(num_laps)x$(num_sites)_spectrum_w0_$(wPh)_nvib_$(vibStates).pdf")

#dt = 0.001

#function DOS(En,HS,sigma)
#    en = eigvals(HS)
#    dos = []
#    for i in 1:length(En)
#        ds = sum([1/(sigma*sqrt(pi))*exp(-(En[i]-ei)^2/(sigma^2)) for ei in en])
#        push!(dos,ds)
#    end
#    return dos
#end
#En_1 = eps[1]-t0-0.01
#En_2 = 0.04
#En_s = wPh/6
#En = En_1:En_s:En_2
#
#dos = DOS(En,H,1e-7)
#plt.plot(En,dos)
#plt.show()
#

#path = ""
#if !isdir(path)
#    mkdir(path)
#end
##filname = "pops_vib_2.dat"
##imagfilname = "pop_vib_im_2.dat"
##pathIn = "pops_vib.dat_final_Quantities"
##
##nsteps = Int(950/dt)#Int(200/dt)
##insteps = Int(10/dt)
##nsave = 100
##nflush= 600
##
##@timev Propagate(psi0,num_sites,vibStates,H,dt,filname,nsteps,nsave,nflush,pathIn)
#PropagateImag(num_sites,vibStates,H,dt,filname,imagfilname,nsteps,insteps,nsave,nflush)




end
