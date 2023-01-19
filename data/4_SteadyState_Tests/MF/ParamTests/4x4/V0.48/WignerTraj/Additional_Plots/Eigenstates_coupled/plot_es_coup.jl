using LinearAlgebra
using PyPlot
src_path=expanduser("/Users/rudolfsmorka/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v2_weighting/")
include("$(src_path)hamiltonian.jl")
include("$(src_path)structures.jl")

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 15
lw = 3
rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

using .Hamiltonian
using .CustomStructures

num_laps = 4
num_ions = 4
N = num_laps * num_ions
t0  = 0.04
eu  = -6*t0
ed  = -6*t0
E   = [eu,ed]
ws  = 0.05*t0*ones(N)
gPh = [0.0,0.0]
a   = 1.0
c   = 1.0
tso = t0/40
tPh = 0.0
tsoPh = 0.0
#model = HolsteinHopping
#modelStr = "Holstein"
model = FranssonHamiltonian
modelStr = "Fransson"
xnIn = zeros(N)
pnIn = zeros(N)

P = 0.5
tL = t0/4 * (ones(2) .+ P * [1,-1])
tR = t0/4 * ones(2)
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
                            xnIn,
                            pnIn,
                            tL,
                            tR)


H = model(HamPar)
V = 0.48
Gamma = zeros(2*N,2*N)
Gamma[1,1] = tL[1]
Gamma[2,2] = tL[2]
Gamma[end-1,end-1] = tR[1]
Gamma[end,end] = tR[2]
Heff = H - 1im*Gamma

fig = plt.figure(figsize=(10,10))

px = 0.105
py = 0.827
sx = 0.87
sy = 0.168
dy = 0.051
dy2 = 0.018

ax0 = plt.axes([px,py,sx,sy])
ax1 = plt.axes([px,py - sy   - dy,sx,sy])
ax2 = plt.axes([px,py - 2*sy - dy - dy2,sx,sy])
ax3 = plt.axes([px,py - 3*sy - dy - 2*dy2,sx,sy])
ax4 = plt.axes([px,py - 4*sy - dy - 3*dy2,sx,sy])

ew,ev = eigen(Heff)
Sites = Array(1:N)
m = 2*N
colors = plt.cm.get_cmap("tab20c",m)
colors = [colors(i) for i in 0:1/(m-1):1]


matId = [1 0; 0 1.]
matSx = [0 1; 1 0.]
matSy = [0 -1im; 1im 0.]
matSz = [1 0; 0 -1.]
mats = [matId,matSx,matSy,matSz]
axes = [ax1,ax2,ax3,ax4]

X = Array(-1:0.01:1)

#Psi2 = zeros(m,N)

spc = [1,9,15,17,23,31]

for (i,col) in zip(Array(1:m),colors)
		
	ax0.plot([real(ew[i])],[imag(ew[i])],linewidth=0,marker="o",color=col,markersize=8)
end

for (i,col) in zip(spc,colors[spc])
	ax0.plot([real(ew[i])],[imag(ew[i])],linewidth=0,marker="x",color="k",markersize=7)
	for (axi,mat) in zip(axes,mats)
		psi2 = [ev[k:k+1,i]' * mat * ev[k:k+1,i] for k in 1:2:2*N-1]
		axi.plot(Sites,real.(psi2),color=col,linewidth=lw)
	end


	#Psi2[i,:] = psi2
	
end

ax0.set_ylim(-0.0032,0.0002)
ax0.set_xlabel(L"$\mathrm{Re}(\varepsilon_k)$")
ax0.set_ylabel(L"$\mathrm{Im}(\varepsilon_k)$")

labels = [L"$\vert \langle \psi_k\vert\psi_k\rangle\vert^2$", 
	  L"$\vert \langle \psi_k\vert\sigma^x\vert\psi_k\rangle\vert^2$",
	  L"$\vert \langle \psi_k\vert\sigma^y\vert\psi_k\rangle\vert^2$",
	  L"$\vert \langle \psi_k\vert\sigma^z\vert\psi_k\rangle\vert^2$"]

#ax1.imshow(Psi2,aspect="auto")
for (axi,lab) in zip(axes,labels)
	axi.set_ylabel(lab)
	axi.set_xlim(0.7,N+0.3)
	if axi != ax4
		axi.set_xticklabels([])
	end
end

annx = 0.01
anny = 0.88

axes = [ax0,ax1,ax2,ax3,ax4]
labels = ["(a)","(b)","(c)","(d)","(e)"]
for (ax,lab) in zip(axes,labels)
	ax.annotate(lab,(annx,anny),xycoords="axes fraction",fontsize=fs)
end

ax4.set_xlabel("site")
plt.savefig("Eigenstates_leads_$(N)_tso$(tso)_V$(V)_P$(P).pdf")



