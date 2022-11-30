using LinearAlgebra
using PyPlot
using PyCall

rc("text", usetex=true)
rc("font",family="serif",serif="Times")
#plt.style.use("classic")

fs = 15
rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

random=pyimport("random")

function Wigner(X,P,beta,w)
   W = zeros(length(P),length(X))
   for i in 1:length(P)
       for j in 1:length(X)
	   W[i,j] = 1/pi * tanh(0.5*beta*w) * exp(-tanh(0.5*beta*w) * (w* X[j]^2 + P[i]^2/w) )
       end
   end
   return W
end


P = Array(-1:0.002:1);
X = Array(-15:0.01:15);
beta = 1/0.025
w = 0.05
W = Wigner(X,P,beta,w)
m = 30; 
c = [(random.choices(P,W[:,j],k=m),random.choices(X,W[i,:],k=m)) for i in [random.randint(1,length(P)) for k in 1:m] for j in [random.randint(1,length(X)) for k in 1:m]];
print(size(c))

fig = plt.figure(figsize=(9,6))
px = 0.09
py = 0.55
sx = 0.8
sy = 0.43
dy = 0.05
dx = 0.02
ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px,py-dy-sy,sx,sy])
cax1 = plt.axes([px+sx+dx,py-dy-sy,0.02,2*sy+dy])

ax1.set_xticklabels([])
ax1.tick_params(direction="in")

cp = plt.cm.jet
cmap = plt.cm.ScalarMappable(cmap=cp)
vmn = minimum(W)
vmx = maximum(W)

ax1.imshow(W,aspect="auto",extent=[X[1],X[end],P[1],P[end]],vmin=vmn,vmax=vmx,cmap="jet")
ax1.annotate("Analytical",(0.05,0.9),xycoords="axes fraction",color="white",fontsize=fs)
ax1.set_ylabel(L"$P$")
cbar = fig.colorbar(cmap, ax=ax1, cax = cax1, aspect=90,pad=.001)
cbar.set_label(L"$p(X,P)$",labelpad=-3)
cbar.ax.tick_params(labelsize=fs)


ax2.hist2d([c[i][2][j] for i in 1:length(c) for j in 1:length(c[i][2])],[c[i][1][j] for i in 1:length(c) for j in 1:length(c[i][1])],bins=25,density=true,range=[[X[1],X[end]],[P[1],P[end]]],cmap="jet")
ax2.tick_params(direction="in")
ax2.annotate(L"Pseudo-random numbers $(n=900)$",(0.05,0.9),xycoords="axes fraction",color="white",fontsize=fs)
ax2.set_ylabel(L"$P$")
ax2.set_xlabel(L"$X$")
plt.savefig("Comparison_Pdist_analytical_numerical.pdf")
