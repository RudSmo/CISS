include("classical_phonons.jl")
num_laps = 1
num_ions = 3
N = num_laps * num_ions
ws  = 0.002*ones(N)
beta = ws[1]/0.025
num_avg = 100 
