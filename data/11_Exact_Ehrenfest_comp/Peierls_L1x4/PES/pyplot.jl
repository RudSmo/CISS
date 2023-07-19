using LinearAlgebra
using PyPlot
using3D() 
src_path=expanduser("/Users/rudolfsmorka/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v3_typed_PVar/")
include("$(src_path)hamiltonian.jl")
include("$(src_path)structures.jl")

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)


num_laps = 2
num_ions = 2
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
tso = 0.0 * T0/40
tPh = 0.1*T0*ones(2)
tsoPh = 0.0 * 0.1*tso
model = FranssonHamiltonian
modelStr = "Fransson"
ws = ws/2
# Simulation Parameters
dt = 0.0002
tend = 16
numstep = Int(round(tend/abs(dt)))
writeMode = "minimal"

popElIn = zeros(2*N)
popElIn[1] = 1.0
popPhIn = zeros(N)
#xnIn,pnIn = [0.76,-3.61,-5.24,2.4],[-7.40,-4.33,1.59,-0.3]
xnIn,pnIn = zeros(N), zeros(N)
## Leads
P = 0.5
tl = 0.0 * T0/4
tL = tl * (ones(2) .+ P * [1,-1])
tR = tl * ones(2)



function PhononHamiltonian(HamPar)
        H = zeros(2*HamPar.N,2*HamPar.N)
        for i in 1:HamPar.N
                H[2*i-1,2*i-1] = 0.5*HamPar.X[i]^2
                H[2*i,2*i] = 0.5*HamPar.X[i]^2
        end
        return H
end


fig = plt.figure(figsize=(13,4))

px = -0.03
py = 0.09
sx = 0.37
sy = 0.95
dx = -0.055
phi1 = 17#12
phi2 = 110#-80

ax1 = plt.axes([px,py,sx,sy],projection="3d")
ax1.view_init(phi1,phi2)
ax2 = plt.axes([px+sx+dx,py,sx,sy],projection="3d")
ax2.view_init(phi1,phi2)
ax3 = plt.axes([px+2*(sx+dx),py,sx,sy],projection="3d")
ax3.view_init(phi1,phi2)
axes = [ax1,ax2,ax3]
parameters = [[tPh,0.0],[0.0*tPh,0.1*tso],[tPh,0.1*tso]]
title_labels = [L"$t_1=4\, \mathrm{meV}, \ \lambda_1=0\, \mathrm{meV}$", L"$t_1=0\, \mathrm{meV}, \ \lambda_1=0.1\, \mathrm{meV}$", L"$t_1=4\, \mathrm{meV}, \ \lambda_1=0.1\, \mathrm{meV}$"]
for (ax,params,title) in zip(axes,parameters,title_labels)

	Xrng = -15:0.05:15
	nx = length(Array(Xrng))
	colors = ["blue","red","green","gray","black","orange"]
	for (col,state) in zip(colors,1:2*N)
		PES = zeros(nx,nx)
		#state = 5
		for (i,X1) in enumerate(Xrng)
			for (j,X2) in enumerate(Xrng)

				HamPar = HamiltonianParams( num_laps,
							    num_ions,
							    N,
							    model,
							    modelStr,
							    E,
							    a,
							    c,
							    t0,
							    params[1],
							    tso,
							    params[2],
							    gPh,
							    ws,
							    [X1,X2,0,0.],
							    pnIn,
							    tL,
							    tR)			


				Hel = model(HamPar)
				Hph = PhononHamiltonian(HamPar)
				H = Hel .+ Hph
				ew = eigvals(H)
				PES[i,j] = ew[state]
			end
		end
		ax.plot_surface(Xrng' .* ones(nx),ones(nx)' .* Xrng,PES,linewidth=0.0,color=col,antialiased=false,alpha = 0.5)
		ax.set_xlim(Xrng[1],Xrng[end])
		
		ax.text2D(0.25,0.83,title,fontsize=fs,transform=ax.transAxes)
		ax.set_xticks([-15,-10,-5,0,5,10,15])
		ax.set_yticks([-15,-10,-5,0,5,10,15])
		ax.set_ylim(Xrng[1],Xrng[end])
		ax.set_xlim(Xrng[1],Xrng[end])
		ax.set_xlabel(L"$X_1$",labelpad=10)
		ax.set_ylabel(L"$X_2$",labelpad=10)
	end
end
ax1.zaxis.set_rotate_label(false)
#ax1.set_zlabel(L"$\varepsilon\,[\mathrm{meV}]$",rotation=90,labelpad=10)
ax3.set_zlabel(L"$\varepsilon\,[\mathrm{meV}]$",rotation=90,labelpad=10)
#ax2.set_zticklabels([])
#ax3.set_zticklabels([])

savefig("PES_$(num_laps)x$(num_ions)_SOC$(tso)_tPh$(tPh)_SOCPh$(tsoPh)_multiState.pdf")
