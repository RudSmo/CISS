src_path=expanduser("~/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v3_typed_PVar/")
include("$(src_path)hamiltonian.jl")
include("$(src_path)structures.jl")

using PyPlot
using LinearAlgebra



# Hamiltonian Parameters
num_laps = 2
num_ions = 4
N = num_laps * num_ions
ws  = 0.002*ones(N)
Hop = 0.04
T0  = Hop/ws[1]  # in eV
tso = T0/40
a   = 1.0
c   = 1.0
Hso = SOC_A_nn(N,tso,a,c,num_laps)

HsoUp = []
HsoDn = []
HsoUD = []
HsoDU = []
for i in 1:N-2
	j = 2*i-1
 	jn = 2*(i+2)-1
	push!(HsoUp,Hso[j,jn]) # 1,5
	push!(HsoDn,Hso[j+1,jn+1]) # 2,6
	push!(HsoUD,Hso[j,jn+1]) # 1,6
	push!(HsoDU,Hso[j+1,jn]) # 2,5
end
sites = 3:N
plt.plot(sites,real.(HsoUp),color="red",label=L"$\uparrow\rightarrow \uparrow$")
plt.plot(sites,real.(HsoDn),color="blue",label=L"$\downarrow\rightarrow \downarrow$",linestyle="dashed")
plt.plot(sites,real.(HsoUD),color="purple",label=L"$\uparrow\rightarrow \downarrow$")
plt.plot(sites,real.(HsoDU),color="brown",label=L"$\downarrow\rightarrow \uparrow$",linestyle="dashed")
plt.ylabel(L"$H_{\mathrm{soc}}^{jj+2}$")
plt.xlabel("site "*L"$j$")
plt.legend()
plt.show()
