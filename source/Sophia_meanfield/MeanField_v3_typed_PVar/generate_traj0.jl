src_path=expanduser("~/Bachelorarbeit/code/CISS/source/Sophia_meanfield/MeanField_v3_typed_swarm_PVar/")
include("$(src_path)classical_phonons.jl")
include("$(src_path)classical_phonons.jl")
using DelimitedFiles

num_laps = 1
num_ions = 2
N = 2*num_laps * num_ions
ws  = 0.002*ones(N)
beta = ws[1]/0.025
num_avg = 1
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
