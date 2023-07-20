src_path=expanduser("~/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v3_typed_swarm_PVar/")
include("$(src_path)classical_phonons.jl")
using DelimitedFiles

num_laps = 2
num_ions = 2
N = num_laps * num_ions
ws  = 0.001*ones(N)
eta = 10.0
wc = ws[1] # central frequency
w_rng = 0.0002:0.0005:0.1
maux = 100 # number of aux oscillators
ws_aux = [w for w in w_rng[1]:(w_rng[end]-w_rng[1])/(maux-1):w_rng[end]]
gaux = [sqrt(pi * eta * w * exp(-w/wc)/2) for w in ws_aux]  # Ohmic spectral density


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
