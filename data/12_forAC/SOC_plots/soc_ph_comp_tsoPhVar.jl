src_path=expanduser("~/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v3_typed_PVar/")
include("$(src_path)hamiltonian.jl")
include("$(src_path)structures.jl")

using PyPlot
using DelimitedFiles
using LinearAlgebra
plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 15
lw = 2

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)


# Hamiltonian Parameters
num_laps = 2
num_ions = 4
N = num_laps * num_ions
ws  = 0.002*ones(N)
Hop = 0.04
T0  = Hop/ws[1]  # in eV
tso = T0/40
xs = readdlm("phonon.dat")[end,2:2:end] 
tsoPh = 0.1 * tso

function SOCPot(num_laps,num_ions,N,ws,Hop,T0,tso,tsoPh,a,c,w,x)

	Hsos = SOC_A_nn(N,tso,a,c,num_laps)
	HsoPh = Peierls_SOC(N,tsoPh,a,c,num_laps,w,x) 
	Hso = Hsos + HsoPh
	HsoUp = []
	HsoDn = []
	HsoUD = []
	HsoDU = []
	for i in 1:N-2
		j = 2*i-1
		jn = 2*(i+2)-1
		#push!(HsoUp,Hso[j,jn]) # 1,5
		#push!(HsoDn,Hso[j+1,jn+1]) # 2,6
		push!(HsoUD,Hso[j,jn+1]) # 1,6
		#push!(HsoDU,Hso[j+1,jn]) # 2,5
	end
	return HsoUD
end



sites = 3:N
tso_params = 0.0:0.05*tso:5*tso
c = 1.0
a = 1.0
HsoR = zeros(length(tso_params),length(sites))
HsoI = zeros(length(tso_params),length(sites)) 

px = 0.05
py = 0.1
sx = 0.38
sy = 0.88
dx = 0.05
dx1= 0.02

fig = plt.figure(figsize=(12,6))
ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px+sx+dx,py,sx,sy])
ax2c = plt.axes([px+2*sx+dx+dx1,py,0.01,sy])

for (i,tsoPh) in enumerate(tso_params)
	Hso = SOCPot(num_laps,num_ions,N,ws,Hop,T0,tso,tsoPh,a,c,ws,xs)
	HsoR[i,:]  = real.(Hso)
	HsoI[i,:]  = imag.(Hso)
end

cp = plt.cm.inferno
cmap = plt.cm.ScalarMappable(cmap=cp)
vm = minimum(HsoR)
vx = maximum(HsoR)

ax1.imshow(HsoR,aspect="auto",vmin=vm,vmax=vx,extent=[3,N,tso_params[end]*ws[1]*1e3,tso_params[1]*ws[1]*1e3],cmap="inferno")
image = ax2.imshow(HsoI,aspect="auto",vmin=vm,vmax=vx,extent=[3,N,tso_params[end]*ws[1]*1e3,tso_params[1]*ws[1]*1e3],cmap="inferno")
ax1.set_ylabel(L"$\lambda_{1}\,[\mathrm{meV}]$")
ax2.set_yticklabels([])
ax1.set_xlabel("site "*L"$j$")
ax2.set_xlabel("site "*L"$j$")
cbar = fig.colorbar(image,ax=ax2,cax = ax2c,aspect=90,pad=.001)
cbar.set_label(L"$H_{\mathrm{soc}}^{\uparrow\downarrow}\,[\mathrm{meV}]$")
plt.savefig("SOC_$(num_laps)x$(num_ions)_tsoPhVar_Ph.pdf")
