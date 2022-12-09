using LinearAlgebra
using PyPlot
using DelimitedFiles

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3
rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)




fig = plt.figure(figsize = (10,4.5))
px = 0.08
py = 0.125
sx = 0.89
dx = 0.1
sy = 0.85
ax = plt.axes([px,py,sx,sy])


Gamma0 = 1/4

I50 = readdlm("currents_avg_n50.dat")
I25 = readdlm("currents_avg_n25.dat")
I5 = readdlm("currents_avg_n5.dat")


ax.plot(I50[2:end,1],(I50[2:end,2]-I50[2:end,3])/Gamma0,color="red", linewidth = lw, label = L"$n=1000$")
ax.plot(I25[2:end,1],(I25[2:end,2]-I25[2:end,3])/Gamma0,color="blue", linewidth = lw, label = L"$n=500$")
ax.plot(I5[2:end,1],(I5[2:end,2]-I5[2:end,3])/Gamma0,color="k", linewidth = lw, label = L"$n=100$")

ax.set_ylim(0.25,0.55)
ax.set_xlabel(L"$t\, [t_0^{-1}]$",labelpad=-2)
ax.set_ylabel(L"I_L-I_R \, [\Gamma_0]")
ax.legend(loc="upper right")
plt.savefig("Currents_1x3_Ehrenfest_V12.pdf")




