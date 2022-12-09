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

numDats = 5
numFils = 20

dat = readdlm("data1/L3_tPh0.1_SOC0.025_SOCPh0.0025000000000000005_V18.0_dt0.002_1_currents.dat")
t = dat[2:end,1]
Gamma0 = 0.25 


let IT = zeros(length(dat[2:end,2]))
	let WT = 0.0

		for i in 1:numDats
			wT = 0.0
			for j in 1:numFils
				data = readdlm("data$(i)/L3_tPh0.1_SOC0.025_SOCPh0.0025000000000000005_V18.0_dt0.002_$(j)_currents.dat")
				I1,I2 = data[2:end,2],data[2:end,3]
				w = readdlm("data$(i)/L3_tPh0.1_SOC0.025_SOCPh0.0025000000000000005_V18.0_dt0.002_$(j)_weight.dat")[1]
				if i == 1 && j == 1
					ax.plot(t,(I1-I2)/Gamma0,color="k",alpha=0.3,label="single trajectory")
				else
					ax.plot(t,(I1-I2)/Gamma0,color="k",alpha=0.3)
				end
				wT += w
				IT .+= w * (I1-I2)/Gamma0
			end
			WT += wT
		end
		IT = IT ./ WT

		ax.plot(t,IT,color="red", linewidth = lw, label = L"weighted avg. $n=100$")
		I50 = readdlm("../currents_avg_n50.dat")
		ax.plot(I50[2:end,1],(I50[2:end,2]-I50[2:end,3]) / Gamma0,color="orange",linestyle="dashed", linewidth = lw, label = L"weighted avg. $n=1000$")
#		ax.plot(I50[2:end,1],14/15 * ones(length(I50[2:end,1])), color = "purple", linestyle = "dashed", linewidth = lw, label = L"$ \sum_{\sigma} \frac{\Gamma_{L\sigma}\Gamma_{R\sigma}}{\Gamma_{L\sigma}+\Gamma_{R\sigma}}$")
		
		ax.set_ylim(0.4,1.0)
		ax.set_xlabel(L"$t\, [t_0^{-1}]$",labelpad=-2)
		ax.set_ylabel(L"I_L-I_R \, [\Gamma_0]")
		ax.legend(loc="lower right")
		plt.savefig("Currents_1x3_Ehrenfest_noAvg_V18.pdf")
	end
end




