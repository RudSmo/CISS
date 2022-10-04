include("../../../hamiltonian.jl")
include("../../../electrons.jl")
using .Hamiltonian
using .Electrons
using PyPlot

E = [zeros(2),zeros(2),zeros(2)];
t = [10.0,1.0,10.0];
gSL = [[0.1 0.; 0. 0.1],[0.1 0; 0 0.1]];
a = ones(3);
c = ones(3);
num_laps = [1,10,1];
ions = [50,4,50]
N = Int.([n*i for (n,i) in zip(num_laps,ions)]);
tph = zeros(3);
tsoph = zeros(3);
qs = [zeros(i) for i in N];
ws = qs;
tsoL = [0.0,4.,0.0]
H = Hamiltonian.FranssonHamiltonianExtended(N,E,t,gSL,tsoL,a,c,num_laps,tph,tsoph,ws,qs)

l = sum(N)
rho = DensityFEx(l,1,H,l,0)

#plt.plot([rho[i,i]-rho[i+1,i+1] for i in 1:2:2*l-1])
plt.plot([rho[i,i] for i in 1:2:2*l],color="red")
plt.plot([rho[i,i] for i in 2:2:2*l],color="blue")
plt.show()


