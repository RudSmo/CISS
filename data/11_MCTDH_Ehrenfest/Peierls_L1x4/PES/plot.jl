using LinearAlgebra
#using PyPlot
#using3D()
using Plots
using LaTeXStrings
plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing, grid=true)
src_path=expanduser("/Users/rudolfsmorka/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v2_weighting/")
include("$(src_path)hamiltonian.jl")
include("$(src_path)structures.jl")

#plt.style.use("classic")
#rc("text", usetex=true)
#rc("font",family="serif",serif="Times")
#fs = 17
#rc("axes",labelsize=fs)
#rc("xtick",labelsize=fs-2)
#rc("ytick",labelsize=fs-2)

using .Hamiltonian
using .CustomStructures

num_laps = 2
num_ions = 3
N = num_laps * num_ions
t0  = 0.04
eu  = -6*t0
ed  = -6*t0
E   = [eu,ed]
ws  = 0.05*ones(N)
gPh = [0.0,0.0]
a   = 1.0
c   = 1.0
tso = t0/40
tPh = 0.1*t0
tsoPh = 0.1*tso
#model = HolsteinHopping
#modelStr = "Holstein"
model = FranssonHamiltonian
modelStr = "Fransson"
xnIn = zeros(N)
pnIn = zeros(N)

P = 0.5
tL = t0/4 * (ones(2) .+ P * [1,-1])
tR = t0/4 * ones(2)


Xrng = Array(-90:1:60)
nx = length(Xrng)


function PlotPES(num_laps,
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
	p = []
	for state in 1:2*N
		PES = zeros(nx,nx)
		#state = 5
		for (i,X1) in enumerate(Xrng)
			for (j,X2) in enumerate(Xrng)
	
				xnIn = zeros(N)
				xnIn[1] = X1
				xnIn[1] = X2
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


				H = model(HamPar)
				ew = eigvals(H)
				PES[i,j] = ew[state]
			end
		end
		if state == 1
			col = cgrad([:brown,:brown]) 
		elseif state == 2
			col = cgrad([:blue,:blue]) 
		elseif state == 3
			col = cgrad([:red,:red]) 
		elseif state == 4
			col = cgrad([:orange,:orange])
		elseif state == 5
			col = cgrad([:black,:black]) 
		elseif state == 6
			col = cgrad([:green,:green]) 
		end
		push!(p,plot(Xrng,Xrng,PES,st=:surface,camera=(10,10),c=col,fillalpha=0.5, legend= false, size = (400,00), xlabel=L"$X_1$", ylabel=L"$X_2$", zlabel=L"$\varepsilon\, [\mathrm{eV}]$"))
		#wireframe!(Xrng,Xrng,PES)
	end
	return plot(p)
end
PlotPES(num_laps,
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

savefig("PES_$(num_laps)x$(num_ions)_SOC$(tso)_tPh$(tPh)_SOCPh$(tsoPh)_multiState.pdf")
#show()
#plt.plot_surface(Xrng, Xrng, PES,cmap=plt.cm.Blues)
#plt.show()
