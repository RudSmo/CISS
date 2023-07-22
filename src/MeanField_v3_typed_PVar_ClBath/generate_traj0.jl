src_path=expanduser("~/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v3_typed_swarm_PVar/")
include("$(src_path)classical_phonons.jl")
using DelimitedFiles

num_laps = 2
num_ions = 2
N = num_laps * num_ions
ws  = 0.001*ones(N)
eta = 0.05
wc = ws[1] # central frequency
maux = 50

w_rng = 1e-7:1e-5:0.015
ws_aux1 = w_rng[1]:(wc/2-w_rng[1])/(maux/4):wc/2
ws_aux2 = (wc/2):(3*wc/2)/(maux/2-1):2*wc
ws_aux3 = (2*wc):(w_rng[end]-3*wc/2)/(maux/4-1):w_rng[end]
ws_aux  = collect(Iterators.flatten([Array(ws_aux1),Array(ws_aux2),Array(ws_aux3)]))
gaux    = [sqrt(pi * eta * w * exp(-w/wc)/2) for w in ws_aux]  # Ohmic spectral density


beta = ws[1]/0.025
num_avg = 50
#X0,P0           = SampleGS_Multi(N,num_avg) #SampleWigner_Multi(beta,ws,N,num_avg)
#Xaux0,Paux0 	= SampleGS_Multi(maux,num_avg) #SampleWigner_Multi(beta,ws_aux,maux,num_avg)
X0,P0           = SampleWigner_Multi(beta,ws,N,num_avg)
Xaux0,Paux0     = SampleWigner_Multi(beta,ws_aux,maux,num_avg)
X_l,P_l         = zeros(num_avg,N),zeros(num_avg,N)
Xaux_l,Paux_l 	= zeros(num_avg,maux),zeros(num_avg,maux)
X_l[:,:]        = X0
P_l[:,:]        = P0
Xaux_l[:,:] 	= Xaux0
Paux_l[:,:] 	= Paux0

version = 0

open("init_X_traj$(version).dat", "w") do io
        writedlm(io, X_l, "\t")
end
open("init_P_traj$(version).dat", "w") do io
        writedlm(io, P_l, "\t")
end
open("init_Xaux_traj$(version).dat", "w") do io
        writedlm(io, Xaux_l, "\t")
end
open("init_Paux_traj$(version).dat", "w") do io
        writedlm(io, Paux_l, "\t")
end
