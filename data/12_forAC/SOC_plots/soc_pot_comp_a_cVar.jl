src_path=expanduser("~/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v3_typed_PVar/")
include("$(src_path)hamiltonian.jl")
include("$(src_path)structures.jl")

using PyPlot
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
T0  = 1e3*Hop #Hop/ws[1]  # in eV
tso = T0/40

function SOCPot(num_laps,num_ions,N,ws,Hop,T0,tso,a,c)

	Hso = SOC_A_nn(N,tso,a,c,num_laps)

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
#	return HsoUp,HsoDn,HsoUD,HsoDU
end

params = [[1.,1.],[T0,T0],[T0,0.1*T0],[0.1*T0,T0],[0.01*T0,T0]]
sites = 3:N

px = 0.1
py = 0.1
sx = 0.38
sy = 0.88
dx = 0.1
fig = plt.figure(figsize=(12,6))
ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px+sx+dx,py,sx,sy])

for param in params
	Hso  = SOCPot(num_laps,num_ions,N,ws,Hop,T0,tso,param[1],param[2]) 
	if param == params[2]
		ax1.plot(sites,real.(Hso),label="a=$(param[1]),c=$(param[2])",linestyle="dashed",linewidth=3)
		ax2.plot(sites,imag.(Hso),label="a=$(param[1]),c=$(param[2])",linestyle="dashed",linewidth=3)
	else
		ax1.plot(sites,real.(Hso),label="a=$(param[1]),c=$(param[2])",linewidth=3)
		ax2.plot(sites,imag.(Hso),label="a=$(param[1]),c=$(param[2])",linewidth=3)
	end
end
ax1.legend()
ax1.set_ylabel(L"$\mathrm{Re}\, H_{\mathrm{soc}}^{jj+2,\uparrow\downarrow}\,[\mathrm{meV}]$")
ax2.set_ylabel(L"$\mathrm{Im}\, H_{\mathrm{soc}}^{jj+2,\uparrow\downarrow}\,[\mathrm{meV}]$")
ax1.set_xlabel("site "*L"$j$")
ax2.set_xlabel("site "*L"$j$")
plt.savefig("SOC_$(num_laps)x$(num_ions)_aVar_cVar.pdf")
