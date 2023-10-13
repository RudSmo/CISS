src_path=expanduser("~/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v3_typed_PVar/")
include("$(src_path)hamiltonian.jl")
include("$(src_path)structures.jl")

using PyPlot
using LinearAlgebra
using DelimitedFiles
# Hamiltonian Parameters
num_laps = 2
num_ions = 4
N = num_laps * num_ions
ws  = 0.002*ones(N)
Hop = 0.04
T0  = Hop/ws[1]  # in eV
t0  = -T0
eu  = -6*T0
ed  = -6*T0
E   = [eu,ed]
gPh = zeros(2)
a   = 1.0
c   = 1.0
tso = T0/40
tPh = 0.1*T0*ones(2)
tsoPh = 0.1*tso
model = FranssonHamiltonian
modelStr = "Fransson"
# Simulation Parameters
dt = 0.0002
tend = 3
numstep = Int(round(tend/abs(dt)))
writeMode = "minimal"


popElIn = zeros(2*N)
popElIn[1] = 0.0
popPhIn = zeros(N)
xnIn,pnIn = zeros(N), zeros(N)
## Leads
P = 0.0
tl = T0/4
tL = tl * [1,0]#(ones(2) .+ P * [1,-1])
tR = tl * [1,0]#ones(2)

v = 85
V = v * 0.2 * T0
nPoles = 15
beta = ws[1]/0.025

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

function PhononHamiltonian(HamPar)
        H = zeros(2*HamPar.N,2*HamPar.N)
        for i in 1:HamPar.N
                H[2*i-1,2*i-1] = 0.5* HamPar.X[i]^2
                H[2*i,2*i] = 0.5*HamPar.X[i]^2
        end
        return H
end


function SimulatedAnnealing(X0,model,HamPar,p_init,p_final,trial,cycle,eps,fac)
	HamPar.X = X0
	acc = 0
	accept = false
	T_init = -1/log(p_init)
	T_final = -1/log(p_final)
	inc = (T_init-T_final)/cycle
	ew = eigvals(model(HamPar))
	E_c = sum([e for e in ew[1:Int(end/2)]])	
	T_l = [T_init]
	E_l = [E_c]
	X_l = [X0]
	X_c = X0
	T = T_init
	for i in 1:cycle
		for j in 1:trial
			X_i = X_c + eps*randn(HamPar.N)
			HamPar.X = X_i
			H_el = model(HamPar)
			H_vib = PhononHamiltonian(HamPar)
			H = H_el + H_vib
			ew = eigvals(H)
			E_i = sum([e*f for (e,f) in zip(ew,fac)])
			r = rand()
			k = exp(-(E_i-E_c)/T)
			if E_i > E_c
				if r < min(1,k)
					accept = true
				else
					accept = false
				end
			else
				accept = true
			end
			if accept == true
				E_c = E_i
				X_c = X_i
				acc += 1
			end
		end
		push!(X_l,X_c)
		push!(E_l,E_c)
		T -= inc
		push!(T_l,T)			
	end
	return E_l,T_l,X_l
end

# findmin(E_l); E_min[1]; X_l[E_min[2]]
fac = readdlm("electrons_es_avg_n21_P.dat")[end,2:end]
X0 = randn(N)
eps = 1e-2
p_init = 0.9
p_final = 0.01
trial = 200
cycle = 12500
E_l,T_l,X_l = SimulatedAnnealing(X0,model,HamPar,p_init,p_final,trial,cycle,eps,fac)
data = zeros(length(E_l),2+length(X_l[1]))
data[:,1] = T_l
data[:,2] = E_l
for i in 1:length(X_l)
	data[i,3:end] = X_l[i]
end
writedlm("L$(num_laps)x$(num_ions)_tso_$(tso)_tPh_$(tPh[1])_tsoPh_$(tsoPh)_simAnn.dat",data)

#println(X_l)
#plt.plot(Array(0:cycle),E_l)
#println(X_l[findmin(E_l)[2]])
#plt.show()

	
#H = model(HamPar)
#Ews = eigvals(H)


#function Dos(E_i,E,b=1e-3)
#    D = [1/(sqrt(2*pi)*b)*sum([exp(-(e-i)^2/(2*b^2)) for i in E_i]) for e in E]
#    return D
#end
#
#
#for ew in Ews
#	plt.plot([-1,1],ew*ones(2)*ws[1],color="k")
#end
#plt.ylabel(L"$\varepsilon \, [\mathrm{eV}]$")
#plt.show()
#
#
