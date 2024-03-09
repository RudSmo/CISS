using PyPlot
using DelimitedFiles


plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3
fac_c = 0.002*0.2*1e3 * 0.2433946
fac_c_SS = 1e6
fac_t = 1/(0.002*0.2*4.135668*1e3)

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.08
py = 0.12
sx = 0.89
sy = 0.85
dy = 0.05
dx = 0.01

fig = plt.figure(figsize=(10,10))
ax = plt.axes([px,py,sx,sy])

sp = "P"
dat = "electrons_avg_n101_$(sp).dat"

N = 21
SP = []
for i in 0:4
        datE = readdlm("../$(i)/$(dat)")
        push!(SP,datE[:,2:2:end] .- datE[:,3:2:end])
end


ls = [length(S[:,1]) for S in SP]
SPP = zeros(sum(ls),N)

SPP[1:ls[1],:] = SP[1]
SPP[ls[1]+1:ls[2]+ls[1],:] = SP[2]
SPP[ls[2]+ls[1]+1:ls[3]+ls[2]+ls[1],:] = SP[3]
SPP[ls[3]+ls[2]+ls[1]+1:ls[4]+ls[3]+ls[2]+ls[1],:] = SP[4]
SPP[ls[4]+ls[3]+ls[2]+ls[1]+1:ls[5]+ls[4]+ls[3]+ls[2]+ls[1],:] = SP[5]


sp = "M"
dat = "electrons_avg_n101_$(sp).dat"

N = 21
SP = []
for i in 0:4
	datE = readdlm("../$(i)/$(dat)")
	push!(SP,datE[:,2:2:end] .- datE[:,3:2:end])
end


ls = [length(S[:,1]) for S in SP]
SPm = zeros(sum(ls),N)

SPm[1:ls[1],:] = SP[1]
SPm[ls[1]+1:ls[2]+ls[1],:] = SP[2]
SPm[ls[2]+ls[1]+1:ls[3]+ls[2]+ls[1],:] = SP[3]
SPm[ls[3]+ls[2]+ls[1]+1:ls[4]+ls[3]+ls[2]+ls[1],:] = SP[4]
SPm[ls[4]+ls[3]+ls[2]+ls[1]+1:ls[5]+ls[4]+ls[3]+ls[2]+ls[1],:] = SP[5]


Dif = SPP + SPm

plt.imshow(Dif,aspect="auto")
ax.set_xlabel(L"site")
ax.set_ylabel(L"$t\,[\mathrm{ps}]$")
plt.savefig("L3x7_V60_elec_sz_PM_Dif.pdf")
