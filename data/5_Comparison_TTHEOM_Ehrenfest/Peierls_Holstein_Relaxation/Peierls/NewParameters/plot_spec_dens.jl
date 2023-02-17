using LinearAlgebra
using PyPlot
using DelimitedFiles

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3
fac_c = 1e3 * 0.2433946
fac_c_SS = 1e6
fac_t = 1/(4.135668)

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.08
py = 0.12
sx = 0.89
sy = 0.85
dy = 0.05

fig = plt.figure(figsize=(10,5.5))
ax = plt.axes([px,py,sx,sy])

H(X,eps,t0,t1,w) = [eps -t0-t1*sqrt(2*w)*X; -t0-t1*sqrt(2*w)*X eps]
eps = 0.0
t0  = 2.0
t1  = 0.2
w0  = 0.1
v   = 7.0

dat = readdlm("L2_gPh0.0_tPh0.2_V$(v)_w00.1_dt0.0005_1_phonons.dat")
En = Array(-4:0.005:4)

t = dat[2:200:end,1]
Xph = dat[2:200:end,2]

el1 = []
el2 = []

SD = zeros(length(En),length(Xph))
Id = Matrix(I,2,2)

Gamma = 0.1*Id
for (i,x) in enumerate(Xph)
	Ham = H(x,eps,t0,t1,w0)
	for (j,E) in enumerate(En)
		Gr = inv(Id * E - Ham + 1im/2 * Gamma)
	        Ga = Gr'
        	Lambda = Gr * Gamma * Ga 
		SD[j,i] = real(tr(Lambda))
	end			
end

ax.imshow(SD,aspect="auto")

plt.savefig("Peierls_L2_NewParams_SpecDens_V$(v).pdf")
