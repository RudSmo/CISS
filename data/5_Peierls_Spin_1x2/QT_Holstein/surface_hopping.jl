using LinearAlgebra

# Define the system Hamiltonian
function hamiltonian(ϵ1, ϵ2, γ)
    return [ϵ1 γ; γ ϵ2]
end


# Define the electronic states
ψ1 = [1, 0]
ψ2 = [0, 1]

# Define the electronic coupling matrix
Γ = [0 1; 1 0]

# Define the nuclear degrees of freedom (e.g. vibrational modes)
q = [0, 0]  # initial nuclear position
p = [0, 0]  # initial nuclear momentum

# Define the electron-phonon coupling strength
α = 0.1
global ϵ1 = 0.0
global ϵ2 = 1.0

# Define the time step and simulation time
dt = 0.01
tmax = 10.0

# Initialize the electronic and nuclear states
ψ = ψ1
Q = [q]

# Initialize the electronic and nuclear velocities
v = zeros(2)
P = [p]

# Initialize the electronic energies
E = [ϵ1]

# Initialize the hopping probability
P_hop = 0.0

# Run the simulation
for t in 0:dt:tmax
    # Compute the current electronic energy
    H = hamiltonian(ϵ1, ϵ2, α*q[1])
    E_now = real(ψ' * H * ψ)
    
    # Compute the force on the nuclear degree of freedom
    F = -α * [ψ1'*Γ*ψ, ψ2'*Γ*ψ]'  # electron-phonon coupling
    
    # Compute the nuclear acceleration
    a = F / m  # assume unit mass
    
    # Update the nuclear position and velocity
    q = q + p * dt + 0.5 * a * dt^2
    p = p + 0.5 * (a + F/m) * dt
    
    # Compute the hopping probability
    P_hop = dt * abs(α * ψ' * [q[1]*ψ1, q[1]*ψ2]' * Γ * [q[1]*ψ2, q[1]*ψ1]' * ψ)
    
    # Decide whether to hop or not
    if rand() < P_hop
        # Update the electronic state and energy
        if ψ == ψ1
            ψ = ψ2
            ϵ1, ϵ2 = ϵ2, ϵ1
        else
            ψ = ψ1
            ϵ1, ϵ2 = ϵ1, ϵ2
        end
        push!(E, real(ψ' * H * ψ))
    else
        push!(E, E_now)
    end
    
    # Update the nuclear positions and velocities
    push!(Q, q)
    push!(P, p)
end

# Compute the current-voltage characteristics
V = range(0.0, stop=1.0, length=length(E))
I = zeros(length(V))
for i in 2:length(V)
    I[i] = (E[i] - E[i-1]) / dt
end

# Plot the results
using PyPlot
plt.plot(V, I)
plt.xlabel("Voltage (V)")
plt.ylabel("Current (a.u.)")
plt.show()


