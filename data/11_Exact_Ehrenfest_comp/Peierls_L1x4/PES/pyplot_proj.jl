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
tso = T0/40
tPh = 0.1*T0*ones(2)
tsoPh = 0.1*tso
model = FranssonHamiltonian
modelStr = "Fransson"
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


fig = plt.figure(figsize=(13,12))

px = -0.01
py = 0.65
sx = 0.3
sy = 0.37
dx = 0.02
dy = -0.03
phi1 = 17#12
th1 = 110#-80
phi2 = 0.0
th2 = 90.0
phi3 = 0
th3 = 180

ax11 = plt.axes([px,py,sx,sy],projection="3d")
ax11.view_init(phi1,th1)
ax12 = plt.axes([px,py-sy-dy,sx,sy],projection="3d")
ax12.view_init(phi2,th2)
ax13 = plt.axes([px,py-2*(sy+dy),sx,sy],projection="3d")
ax13.view_init(phi3,th3)

ax21 = plt.axes([px+sx+dx,py,sx,sy],projection="3d")
ax21.view_init(phi1,th1)
ax22 = plt.axes([px+sx+dx,py-sy-dy,sx,sy],projection="3d")
ax22.view_init(phi2,th2)
ax23 = plt.axes([px+sx+dx,py-2*(sy+dy),sx,sy],projection="3d")
ax23.view_init(phi3,th3)

ax31 = plt.axes([px+2*(sx+dx),py,sx,sy],projection="3d")
ax31.view_init(phi1,th1)
ax32 = plt.axes([px+2*(sx+dx),py-sy-dy,sx,sy],projection="3d")
ax32.view_init(phi2,th2)
ax33 = plt.axes([px+2*(sx+dx),py-2*(sy+dy),sx,sy],projection="3d")
ax33.view_init(phi3,th3)

axes = [[ax11,ax12,ax13],
	[ax21,ax22,ax23],
	[ax31,ax32,ax33]]
parameters = [[tPh,tsoPh],[0.0*tPh,tsoPh],[tPh,tsoPh]]
title_labels = [L"$t_1=4\, \mathrm{meV}, \ \lambda_1=0\, \mathrm{meV}$", L"$t_1=0\, \mathrm{meV}, \ \lambda_1=4\, \mathrm{meV}$", L"$t_1=4\, \mathrm{meV}, \ \lambda_1=4\, \mathrm{meV}$"]
for (axi,params,title) in zip(axes,parameters,title_labels)

	Xrng = -15:0.05:15 #-500:5:500
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
                                                            tPh,
                                                            tso,
                                                            tsoPh,
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

		for ax in axi
			if ax == ax11 
				ann = L"(a$_1$)"
			elseif ax == ax12 
				ann = L"(a$_2$)"
			elseif ax == ax13
				ann = L"(a$_3$)"
			elseif ax == ax21
                                ann = L"(b$_1$)"
                        elseif ax == ax22
                                ann = L"(b$_2$)"	
			elseif ax == ax23
                                ann = L"(b$_3$)"
                        elseif ax == ax31
                                ann = L"(c$_1$)"
			elseif ax == ax32
                                ann = L"(c$_2$)"
                        elseif ax == ax33
                                ann = L"(c$_3$)"
			end
			ax.text2D(0.05,0.9,ann,fontsize=fs,transform=ax.transAxes)
			if ax == ax11 || ax == ax21 || ax == ax31
				ax.text2D(0.25,0.87,title,fontsize=fs,transform=ax.transAxes)
				ax.set_xticks([-400,-200,0,200,400])
				ax.set_yticks([-400,-200,0,200,400])
				ax.set_xlabel(L"$X_1$",labelpad=7)
				ax.set_ylabel(L"$X_2$",labelpad=7)
			elseif ax == ax12 || ax == ax22 || ax == ax32
				ax.set_yticks([])
				ax.set_xticks([-400,-200,0,200,400])
				ax.set_xlabel(L"$X_1$",labelpad=10)
				ax.set_ylabel(L"$X_2$",labelpad=5)
			elseif ax == ax13 || ax == ax23 || ax == ax33
				ax.set_xticks([])
				ax.set_yticks([-400,-200,0,200,400])
				ax.set_xlabel(L"$X_1$",labelpad=5)
				ax.set_ylabel(L"$X_2$",labelpad=10)
			end
			ax.plot_surface(Xrng' .* ones(nx),ones(nx)' .* Xrng,PES,linewidth=0.0,color=col,antialiased=false,alpha = 0.5)
			ax.set_xlim(Xrng[1],Xrng[end])
			
			ax.set_ylim(Xrng[1],Xrng[end])
			ax.set_xlim(Xrng[1],Xrng[end])
			#ax.set_zlim(-0.3,0.05)
		end
	end
end

ax31.set_zlabel(L"$\varepsilon\,[\mathrm{eV}]$",rotation=90,labelpad=10)
ax32.set_zlabel(L"$\varepsilon\,[\mathrm{eV}]$",rotation=90,labelpad=10)
ax33.set_zlabel(L"$\varepsilon\,[\mathrm{eV}]$",rotation=90,labelpad=10)

savefig("PES_$(num_laps)x$(num_ions)_SOC$(tso)_tPh$(tPh)_SOCPh$(tsoPh)_multiState_multiAngle.pdf")
