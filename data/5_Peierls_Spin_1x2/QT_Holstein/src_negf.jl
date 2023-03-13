using LinearAlgebra

# Define the Hamiltonian for a two-site Peierls-type system
function hamiltonian(epsilon, t, lambda, gamma)
    H = [epsilon -t;
         -t -epsilon]
    V = [0 lambda;
         lambda 0]
    H + gamma * V
end

# Define the current operator
function current_operator()
    [0 1;
     -1 0]
end

# Define the function that computes the current-voltage characteristics
function current_voltage(V, gamma, epsilon, t, lambda, temperature)
    kT = temperature * 8.617333262145e-5 # Boltzmann constant times temperature in eV
    ID = Matrix(I,2,2)
    current = []
    for i in 1:length(V)
        G = inv((1im * V[i] + kT) * ID - hamiltonian(epsilon, t, lambda, gamma))
        f = 1.0 / (1.0 + exp((V[i] - epsilon) / kT))
        I = real(tr(current_operator() * (G - G')))
        push!(current, 2.0 * pi * gamma * f * I)
    end
    current
end

# Compute the current-voltage characteristics for a system with epsilon = 0.5 eV, t = 0.1 eV,
# lambda = 0.05 eV, gamma = 0.01 eV, and temperature = 300 K, over a voltage range of -2 eV to 2 eV
Volts = -2:0.01:2
current = current_voltage(Volts, 0.01, 0.5, 0.1, 0.05, 300)

# Plot the results
using PyPlot
plt.plot(Volts, current)
plt.xlabel("Voltage (eV)")
plt.ylabel("Current (A)")
plt.show()
