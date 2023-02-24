module InputFile


include("rk4.jl")
using .RK4
using DelimitedFiles
using LinearAlgebra

w0 = 2.0
V = 20
dat = readdlm("../../../data/5_Peierls_noSpin_1x2/MFHEOM/NewParameters_I_V_avg/V$(V)/electrons_es_avg_n30.dat")
#dat = readdlm("../../../data/5_Peierls_noSpin_1x2/MFHEOM/eps_-240_w0_2_t0_40_t1_4_avg/V$(V)/electrons_es_avg_n30.dat")
dif = dat[2:end,2]-dat[2:end,3]
t1 = 0.2
beta = 1/0.005
alpha = tanh(beta*w0/2)
dx = 0.1
Xrng = Array(-50:dx:50)
psi0 = alpha/pi * exp.(-alpha/pi * w0^2 .* Xrng .^ 2) 
psi0 = psi0/norm(psi0)

dt = dat[3,1]-dat[2,1]

#path = "eps_-240_w0_2_t0_40_t1_4/"
path = "NewParams_I_V/"
if !isdir(path)
    mkdir(path)
end
filname = "$(path)V$(V)_t1_$(t1)_w0_$(w0)_dx$(dx).dat"

Tend = dat[end,1]
nsteps = Int(Tend/dt)
nsave = 100
nflush= 1000

Propagate(psi0,w0,t1,dif,dt,dx,Xrng,filname,nsteps,nsave,nflush)
end
