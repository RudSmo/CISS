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

#ax2 = plt.axes([px+sx+dx,py,sx,sy])

Vs = [0,50,80,100,120,140,160,180,200]
IV = []
Gamma0 = 1/4
avg = 5

function CurAnalyse(I)
	a = 1000 
	I = I[a:end]
	max = maximum(I)
	min = minimum(I)
	mean = 1/(length(I)) * sum(I)
	return [max,min,mean]
end


for V in Vs
	data = readdlm("V$(V)/currents_avg_n$(avg).dat")
	push!(IV,CurAnalyse((data[2:end,2]-data[2:end,3]) ./ Gamma0))
end

Vs = Vs / 10



ax.plot(Vs,[IV[i][1] for i in 1:length(IV)],color="red", linewidth = lw, label = "max")
ax.plot(Vs,[IV[i][2] for i in 1:length(IV)],color="blue", linewidth = lw, label = "min")
ax.plot(Vs,[IV[i][3] for i in 1:length(IV)],color="k", linewidth = lw, label = "mean")
ax.plot(Vs,14/15 * ones(length(Vs)),color = "purple", linewidth = lw, linestyle = "dashed", label =L"$\sum_\sigma \frac{\Gamma_{L\sigma}\Gamma_{R\sigma}}{\Gamma_{L\sigma}+ \Gamma_{R\sigma}}$")

datLB = readdlm("../../../../noPh/LandauerBuettiker/Fransson/1x3/fS3_nL1_tso0.025_P0.5_currents_coarse.dat")
ax.plot(datLB[2:end,1],datLB[2:end,2]/Gamma0,color="orange",linewidth=lw, label = "Landauer-Büttiker")

ax.set_xlabel(L"$V\, [t_0]$")
ax.set_ylabel(L"I_L-I_R \, [\Gamma_0]")
ax.set_ylim(-0.1,1.0)
ax.set_xlim(0,20)
ax.legend(loc="lower right")
plt.savefig("Currents_1x3_Ehrenfest_P0.5_IV_min_max_mean.pdf")





