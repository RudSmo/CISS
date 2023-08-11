include("hamiltonian.jl")
using PyPlot
using LinearAlgebra
using SparseArrays
using Combinatorics 
using StatsBase



num_sites = 3
num_laps = 1
N = num_sites * num_laps
vibStates = 4
T0 = 0.04
t0 = -T0
E =  -6*T0
eps = [E,E]
tPh = 0.1*t0
λ0 = T0/40
λPh = 0.1 * λ0
wPh = 0.002
r_helix = 1.0
c_helix = 1.0

P = 0.5
tl = T0/4
tL = tl * (ones(2) .+ P * [1,-1])
tR = tl * ones(2)

E = -45*T0:0.005*T0:5*T0

lambda_rng = 0:0.01*λ0:λ0
DOS_v = zeros(length(E),length(lambda_rng))

function Dos(E_i,E,b=1e-3)
    D = [1/(sqrt(2*pi)*b)*sum([exp(-(e-i)^2/(2*b^2)) for i in E_i]) for e in E]
    return D
end
for (ind,λPh) in enumerate(lambda_rng)
        H = Fransson(num_sites,num_laps,vibStates,eps,t0,tPh,λ0,λPh,wPh,r_helix,c_helix,false)
        ew = eigvals(H)

        #GL,GR = Gamma(tL,tR,N,vibStates)

        #Heff = H - 1im/2 * (GL+GR)

        #ew = eigvals(Heff)
        #ew_r = real(ew)
        #ew_i = imag(ew)


        #D = Dos(ew,E)
        #plt.plot(E,D,color="blue")
        #D_r = Dos(ew_r,E)
        #D_i = Dos(ew_i,E)
        #plt.plot(E,D_r,label="vib=$(vibStates)")
        #plt.plot(E,D_i,color="blue")
        DOS_v[:,ind] = Dos(ew,E)
end
plt.imshow(DOS_v,aspect="auto")
#plt.legend()
plt.show()




