using LinearAlgebra

# Set model parameters
t_left = 0.25
t_right = 0.25
t_imp = 1.0
delta = 1e5
lambda = 0.2
epsilon1 = 0.1
epsilon2 = 0.12
Gamma_left(energy) = 2*t_left^2 / (energy^2 + delta^2)
Gamma_right(energy) = 2*t_right^2 / (energy^2 + delta^2)
Sigma_imp(energy) = [epsilon1+lambda^2*Gamma_left(energy)  lambda*Gamma_left(energy);
                     lambda*Gamma_right(energy)        epsilon2+lambda^2*Gamma_right(energy)]

# Define NEGF functions
function NEGF(energy)
    Sigma = Sigma_imp(energy)
    G = inv(energy*Matrix(I,2,2) - Sigma)
    SigmaG = Sigma * G
    return G, SigmaG
end

function current(energy, mu_left, mu_right)
    G, SigmaG = NEGF(energy)
    GammaL = -2*imag.(SigmaG[1,1] - SigmaG[2,1])
    GammaR = -2*imag.(SigmaG[2,2] - SigmaG[1,2])
    return 2*(GammaL*(fermi(energy, mu_left)-fermi(energy, mu_right)))
end

function fermi(energy, mu)
    return 1/(1+exp((energy-mu)))
end

# Define simulation parameters
n_energy = 100
energy_min = -5.0
energy_max = 5.0
mu_left = -1.0
mu_right = 1.0
V_bias = mu_right - mu_left

# Compute current as a function of energy
energies = range(energy_min, energy_max, length=n_energy)
currents = [current(energy, mu_left, mu_right) for energy in energies]

# Plot results
using PyPlot
println(currents)
plt.plot(energies, currents)
plt.xlabel("Energy")
plt.ylabel("Current")
plt.show()

