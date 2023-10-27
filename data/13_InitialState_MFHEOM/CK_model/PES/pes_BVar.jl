using LinearAlgebra
using PyPlot
using3D()

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3
rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)


function Hamiltonian(e1,e2,t,l1,l2,l12,B,w1,w2,X1,X2,k)
	Hel = [e1+l1*X1 t+l12*X1+1im*B*X2;
	       t+l12*X1-1im*B*X2 e2+l2*X1]
	return Hel 
end

fig = plt.figure(figsize=(10,6))
px = 0.09
py = 0.09
sx = 0.8
sy = 0.8
dx = 0.01
dy = 0.04

phi1 = 20
phi2 = 120
ax1 = plt.axes([px,py,sx,sy],projection="3d")
ax1.view_init(phi1,phi2)

Xrng = -0.3:0.005:0.3
nx = length(Array(Xrng))

e1 = 400
e2 = 300 
w1 = 100
w2 = 100
l1 = 80
l2 = -80
l12 = 80
B = -80
k = 40
t = 0

colors= ["blue","green"]
for (col,state) in zip(colors,1:2)
	PES = zeros(nx,nx)
	for (i,X1) in enumerate(Xrng)
		for (j,X2) in enumerate(Xrng)
			Hp = Hamiltonian(e1,e2,t,l1,l2,l12,B,w1,w2,X1,X2,k)
			Hm = Hamiltonian(e1,e2,t,l1,l2,l12,-B,w1,w2,X1,X2,k)
			ewP = eigvals(Hp)[state]
			ewM = eigvals(Hm)[state]
			PES[i,j] = ewP - ewM
		end
	end
	ax1.plot_surface(Xrng' .* ones(nx),ones(nx)' .* Xrng,PES,linewidth=0.0,color=col,antialiased=false,alpha = 0.5)
	ax1.set_xlabel(L"$X_1$",labelpad=10)
	ax1.set_ylabel(L"$X_2$",labelpad=10)
end

plt.savefig("CK_origin_params_pes_B_diff.pdf")
