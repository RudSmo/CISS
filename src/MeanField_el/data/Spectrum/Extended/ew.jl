module InputFile

include("../../../hamiltonian.jl")

using .Hamiltonian
using LinearAlgebra

# Geometry Parameters
E = [zeros(2),ones(2),zeros(2)]; 
t = [10.0,1.0,10.0]; 
gSL = [[1.0 0; 0 1],[1.0 0; 0 1]]; 
a = ones(3); 
c = ones(3); 
num_laps = [1,8,1];
ions = [50,4,50]
N = Int.([n*i for (n,i) in zip(num_laps,ions)]); 
tph = zeros(3); 
tsoph = zeros(3); 
qs = [zeros(i) for i in N]; 
ws = qs; 

Tso = Array(0.0:0.01:20.0)
fil = open("extended_N1$(N[1])_NS$(N[2])_N2$(N[3])_eigvals_laps$(num_laps[2])_ions$(ions[2])_a$(a[2])_c$(c[2])_tsoVar.dat","w")
for tso in Tso
	tsoL = [0.0,tso,0.0]
	H = Hamiltonian.FranssonHamiltonianExtended(N,E,t,gSL,tsoL,a,c,num_laps,tph,tsoph,ws,qs)
	ew = eigvals(H)
	write(fil,"$(tso) ")
	for j in 1:length(ew)
		write(fil,"$(ew[j]) ")
	end
	write(fil,"\n")
end
close(fil)


end

