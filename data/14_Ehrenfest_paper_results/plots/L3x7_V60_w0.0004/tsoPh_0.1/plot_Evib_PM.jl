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

fig = plt.figure(figsize=(10,5.5))
ax = plt.axes([px,py,sx,sy])

nameP = "phonons_avg_n162_P.dat"
nameM = "phonons_avg_n162_M.dat"

function Evib(data)
	N = length(data[:,1])
	Ener = zeros(N)
	nv = length(data[1,:])
	for k in 1:N
		for j in 2:nv
			Ener[k] += data[k,j] ^ 2
		end
	end
	return Ener
end

for i in 0:1
	data = readdlm("../$(i)/$(nameP)")
	datb = readdlm("../$(i)/$(nameM)")
	Eviba = Evib(data)
	Evibb = Evib(datb)
	if i == 0

		plt.plot(data[:,1]*fac_t,Eviba,color="red",linewidth=lw,label=L"$A=+\frac{1}{2}$")
		plt.plot(datb[:,1]*fac_t,Evibb,color="blue",linewidth=lw,label=L"$A=-\frac{1}{2}$")
	else
                plt.plot(data[:,1]*fac_t,Eviba,color="red",linewidth=lw)
		plt.plot(datb[:,1]*fac_t,Evibb,color="blue",linewidth=lw)
	end
end
plt.legend()
plt.xlim(0,120)
ax.set_xlabel(L"$t\,[\mathrm{ps}]$")
ax.set_ylabel(L"$E_\mathrm{vib}\,[\hbar \omega]$")
plt.savefig("L3x7_V60_Evib_PM.pdf")
