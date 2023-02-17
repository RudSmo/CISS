module InputFile
p = expanduser("~/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v2_weighting_eigenstates/")
include("$(p)hamiltonian.jl")
include("$(p)classical_phonons.jl")
include("$(p)structures.jl")
using .Hamiltonian
using LinearAlgebra
using PyPlot
using .ClassicalPhonons
using .CustomStructures
# -------------- Parameters --------------------------

# Hamiltonian Parameters
num_laps = 1
num_ions = 2
N = num_laps * num_ions
t0  = 0.04
eu  = -6*t0
ed  = -6*t0
E   = [eu,ed]
ws  = 0.05*t0*ones(N)
gPh = [0.0,0.0]
a   = 1.0
c   = 1.0
tso = t0/40
tPh = 0.0
tsoPh = 0.0

model = FranssonHamiltonian_H0
modelStr = "Fransson"
xnIn,pnIn = zeros(N),zeros(N)


## Leads
P = 0.5
tL = t0/4 * (ones(2) .+ P * [1,-1])  
tR = t0/4 * ones(2)
V  = 12*t0
nPoles = 10
beta = 1/0.025
HamPar = HamiltonianParams( num_laps,
                            num_ions,
                            N,
                            model,
                            modelStr,
                            E,
                            a,
                            c,
                            t0,
                            tPh,
                            tso,
                            tsoPh,
                            gPh,
                            ws,
                            xnIn,
                            pnIn,
			    tL,
			    tR)

H0 = model(HamPar)
ew0,ev0 = eigen(H0)

model1 = FranssonHamiltonian_HI
tPh1 = 0.0 * t0
tsoPh1 = 10 * tso
j = 1
Xrng = Array(-50:0.5:50)
for i in 1:length(ew0)
	ol = []
	for x in Xrng

		HamPar1 = HamiltonianParams( num_laps,
                	            num_ions,
                        	    N,
        	                    model1,
	                            modelStr,
	                            E,
        	                    a,
                	            c,
                        	    t0,
	                            tPh1,
        	                    tso,
                	            tsoPh1,
                        	    gPh,
	                            ws,
        	                    [x,0],
                	            pnIn,
                        	    tL,
	                            tR)
	
		HI = model1(HamPar1)
		push!(ol,ev0[:,j]' * HI * ev0[:,i])
	end
	plt.plot(Xrng,ol,label="$(j)-$(i)")
end
plt.legend()
plt.show()



end

