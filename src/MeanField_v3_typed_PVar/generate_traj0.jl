src_path=expanduser("~/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v3_typed_swarm_PVar/")
include("$(src_path)classical_phonons.jl")
using DelimitedFiles

num_laps = 2
num_ions = 4
N = num_laps * num_ions
ws  = 0.002*ones(N)
beta = ws[1]/0.025
num_avg = 5
X0,P0           = SampleWigner_Multi(beta,ws,N,num_avg)
X_l,P_l         = zeros(num_avg,N),zeros(num_avg,N)
X_l[:,:]        = X0
P_l[:,:]        = P0

version = 0

open("init_X_traj$(version).dat", "w") do io
        writedlm(io, X_l, "\t")
end
open("init_P_traj$(version).dat", "w") do io
        writedlm(io, P_l, "\t")
end
