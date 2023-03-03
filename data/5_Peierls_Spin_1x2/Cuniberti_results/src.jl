using PyPlot
using LinearAlgebra

a = 0.7 # in nm
b = 3.2 # in nm
Q = 2*pi/b # in 1/nm
dz = 0.32 # in nm 
alpha = 2 # in meV

function f(n,N,L0)
	g = 0
	for l in -L0/2:L0/2
		for j in 0:N-1
			g += exp(-1im * Q * j* dz) * 3 * (b*l-j*dz-n*dz)/(a^2 * (1 + (n*dz-l*b-j*dz)^2/a^2 )^(5/2))
		end
	end
	return g
end
function Psi(n,N,L0)
	psi = 0
	for l in -L0/2:L0/2
		for j in 0:N-1
			psi += exp(-1im * Q * j* dz)*1/( ( 1 + (n*dz-l*b-j*dz)^2/a^2 )^(3/2))
		end
	end
	return psi
end

ns = 4 # number of sites in a lap
L0 = Array(1:30) # number of turns
W00 = []
n = 1
for l0 in L0
	w00 = -alpha*f(n,ns*l0,l0)
	push!(W00,w00)
end
plt.plot(L0,abs.(W00))
plt.show()


