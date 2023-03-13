using LinearAlgebra
using Arpack

# Define system parameters
const t0 = 1.0
const omega = 0.5
const alpha = 0.1
const beta = 0.2
const gamma = 0.3
const V = 0.4

function eye(n)
	return Matrix(I,n,n)
end

# Define Hamiltonian matrix
function hamiltonian_mat(n::Int64)
    H = zeros(n, n)
    for i = 1:n-1
        H[i, i+1] = t0
        H[i+1, i] = t0
    end
    for i = 1:n
        H[i, i] = omega * (i - 0.5) + alpha * i^2
    end
    return H
end

# Define self-energy matrix
function self_energy_mat(n::Int64, energy::Float64)
    Sigma = zeros(n, n)*(1+0*1im)
    for i = 1:n
        Sigma[i, i] = V^2 / (energy - omega * (i - 0.5) - alpha * i^2 + 1im * beta)
    end
    return Sigma
end

# Define Green's function using Krylov subspace method
function greens_function(n::Int64, energy::Float64)
    H = hamiltonian_mat(n)
    Sigma = self_energy_mat(n, energy)
    A = inv(energy * eye(n) - H - Sigma)
    b = zeros(n)*(1+0*1im)
    b[1] = 1.0
    F = A * b
    λ, v = eigs(A, nev=1, which=:SM, v0=b)
    return F, -λ[1]
end

# Compute transmission coefficient using Landauer formula
function transmission_coefficient(n::Int64, energy_range::Tuple{Float64, Float64}, delta::Float64)
    E_min, E_max = energy_range
    T = 0.0
    energy_list = E_min:delta:E_max
    for energy in energy_list
        F, G = greens_function(n, energy)
        T += gamma * abs(F[n])^2 * G / (2 * pi)
    end
    return T * delta
end

# Compute and print transmission coefficient
T = transmission_coefficient(50, (-10.0, 10.0), 0.01)
println("Transmission coefficient: ", T)

