using LinearAlgebra
using DelimitedFiles
using PyPlot 

function DOS(En,Spec,sigma=1e-4)
        return 1/(sqrt(4*pi*sigma)*length(Spec))*sum([exp(-(En-i)^2/(4*sigma)) for i in Spec])
end


En = Array(-3.0:0.001:3.0)
dat = readdlm("eigvals_laps10_ions8_eu0.0_ed0.0_a1.0_c1.0_tsoVar.dat")
Tsos = Array(0.5:0.05:1.7)
for (ind,tso) in enumerate(Tsos)
	tso_ind = findall(x->x==tso,dat[:,1])
	ew_tso = dat[tso_ind,2:end]
	dos_tso = [DOS(E,ew_tso) for E in En]
	plt.plot(En,dos_tso .- 1.2*ind,label=L"$\lambda=$"*"$(tso)")
end
plt.legend()
plt.show()

