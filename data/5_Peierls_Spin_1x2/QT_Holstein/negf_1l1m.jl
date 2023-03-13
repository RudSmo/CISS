using LinearAlgebra

# Define the Hamiltonian for the cluster
function hamiltonian(e, v, λ, M)
    c = [0 1; 0 0.]
    cd = c'
    b = zeros(M,M)
    for i in 1:M-1
	b[i,i+1] = sqrt(i)
    end
    bd = b'
    IdB = Matrix(I,M,M)
    IdF = Matrix(I,2,2)
    H = e * kron(cd * c,IdB) + λ * kron(cd * c, bd+b) + v * kron(IdF,bd*b)
    return H
end

function Gamma(gammaL,gammaR,M)
	IdB = Matrix(I,M,M)
	G = zeros(2,2)
	G[1,1] = gammaL
	G[2,2] = gammaR
	return kron(G,IdB)
end

# Define the Fermi function
function fermi(E, μ, T)
    return 1.0/(1.0 + exp((E - μ)/T))
end

# Define the retarded Green's function
function green_function(E, H, Σ, eps = 1e-10)
    ID = Matrix(I,size(H))
    return inv(E*ID - H - Σ .+ 1im*eps)
end

# Define the broadening functions
function broadening(Gamma)
    return -1im * Gamma/2
end

# Calculate the current using the Landauer-Büttiker formula
function current(e, v, λ, ΓL, ΓR, T, μL, μR, Erng, M)
    h = 4.135667696*1e-15
    con = (2.0*1.602176634e-19/h)
    H = hamiltonian(e, v, λ, M)
    Cur = 0.0
    for i in 1:length(Erng)-1
	    E1,E2 = Erng[i],Erng[i+1]
	    GL = Gamma(ΓL,0.0,M) 
            GR = Gamma(0.0,ΓR,M)
	    ΣL = broadening(GL)
	    ΣR = broadening(GR)
	    G1 = green_function(E1, H, ΣL + ΣR)
	    T1 = tr(ΓL*G1*ΓR*G1')
	    dF1 = fermi(E1, μL, T) - fermi(E1, μL, T)

            G2 = green_function(E2, H, ΣL + ΣR)
            T2 = tr(ΓL*G2*ΓR*G2')
            dF2 = fermi(E2, μL, T) - fermi(E2, μL, T)
	    Cur += -1im * con * (E2-E1) * 0.5 * (T1 * dF1 + T2 * dF2)
    end
    return Cur
end

# Define the parameters
e = 0.4 # electronic energy
v = 0.1 # vibrational energy
M = 3 # number of vibrational states
λ = 0.01 # electron-vibration coupling strength
ΓL = 0.25 # left lead broadening
ΓR = 0.25 # right lead broadening
T = 300.0 # temperature
μL = 0.0 # left lead chemical potential
μR = 0.0 # right lead chemical potential
V = range(0.0, 3.0, length=100) # voltage range
Erng = range(-1.5, 1.5, length=201) # energy range for integration

# Calculate the current-voltage characteristics
cur = []
for v in V
    push!(cur, current(e, v, λ, ΓL, ΓR, T, μL + v/2, μR - v/2, Erng, M))
end

using PyPlot
plt.plot(V, real.(cur))
plt.xlabel("Voltage (V)")
plt.ylabel("Current (A)")
plt.show()
