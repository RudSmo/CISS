using PyPlot
using LinearAlgebra
using DelimitedFiles


# In this file, the eigenvalues of input are read for each realization (one row)
# then the DOS is computed for each realization using a Gaussian
# and then the mean average of the DOS of all realizations is taken. 


input = readdlm("DOS_N_50_tso_0.0_tPhU_-0.004_tPhD_0.004_tsoPh_0.0_beta_40.0_tl_0.004_nS_1000_P_0.0_coupled.dat")

function DOS(N,eigv,Erng,sigma = 1e-7)
	return [1/(sqrt(2*pi)*sigma*N)*sum([exp(-(E-i)^2/(2*sigma^2)) for i in eigv]) for E in Erng]
end


Erng = -0.1:0.0001:0.1
N = length(input[1,2:end])
nS = length(input[1:end,1])
DOS_avg = DOS(N,input[1,2:end],Erng)
for i in 2:nS
	DOS_avg .+= DOS(N,input[i,2:end],Erng)
end
DOS_avg = DOS_avg ./ nS
plt.plot(Erng,DOS_avg)
plt.show()

