src_path=expanduser("~/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v3_typed_swarm_PVar/")
include("$(src_path)classical_phonons.jl")
using DelimitedFiles

num_laps = 3
num_ions = 7
N = num_laps * num_ions
ws  = 0.2*0.002*ones(N)
beta = ws[1]/0.025
num_avg = 100
version = 1

for n in 1:version
        X0,P0  = SampleWigner_Multi(beta,ws,N,num_avg)
        open("init_X_traj$(n).dat", "w") do io
                for j in 1:num_avg
                        writedlm(io, X0)
                end
        end
        open("init_P_traj$(n).dat", "w") do io
                for j in 1:num_avg
                        writedlm(io, P0)
                end
        end

end

