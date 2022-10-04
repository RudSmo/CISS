module InputFile

include("../../../hamiltonian.jl")

using .Hamiltonian
using LinearAlgebra
using PyPlot

# Geometry Parameters
E = [zeros(2),zeros(2),zeros(2)]; 
t = [1.0,1.0,1.0]; 
gSL = [[1.5 0; 0 1.5],[1.5 0; 0 1.5]]; 
a = ones(3); 
c = ones(3); 
num_laps = [1,4,1];
ions = [25,4,25]
N = Int.([n*i for (n,i) in zip(num_laps,ions)]); 
tph = zeros(3); 
tsoph = zeros(3); 
qs = [zeros(i) for i in N]; 
ws = qs; 
Tso = [1.0]# Array(0.0:0.001:4.0)
Nd = zeros(length(Tso),sum(N))
fil = open("extended_N1$(N[1])_NS$(N[2])_N2$(N[3])_eigvals_laps$(num_laps[2])_ions$(ions[2])_a$(a[2])_c$(c[2])_tsoVar.dat","w")
ev_inds = [Int(sum(N))]#Array(1:Int(sum(N)))#Int(sum(N))
sz = [1 0; 0 -1.]
id = [1 0; 0 1.]
op = id
for (ind,tso) in enumerate(Tso)
	tsoL = [0.0,tso,0.0]
	H = Hamiltonian.FranssonHamiltonianExtended(N,E,t,gSL,tsoL,a,c,num_laps,tph,tsoph,ws,qs)
	plt.imshow(real.(H),aspect="auto")
	plt.show()
	plt.imshow(imag.(H),aspect="auto")
	plt.show()
	ew,ev = eigen(H)
	ds = [sum([ev[i:i+1,ev_ind]'*op*ev[i:i+1,ev_ind] for ev_ind in ev_inds]) for i in 1:2:length(ew)]
	Nd[ind,:] = ds[:]
end
plt.imshow(Nd,aspect="auto")
plt.show()


end

