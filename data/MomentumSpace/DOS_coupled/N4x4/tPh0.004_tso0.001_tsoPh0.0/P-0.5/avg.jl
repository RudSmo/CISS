using LinearAlgebra
using DelimitedFiles

filname = "GF_DOS_N_16_eps_-0.24_tso_0.001_tPhU_-0.004_tPhD_-0.004_tsoPh_0.0_ws_0.002_beta_40.0_tl_0.004_nS_1000_P_-0.5_coupled"

numDat = 100
numAvg = 1000
totalAvg = numDat * numAvg

outfilname = "$(filname)_avg_$(totalAvg).dat"
dos0 = readdlm("data1/$(filname).dat")
outfil = zeros(length(dos0[:,1]),3)
outfil[:,1] = dos0[:,1]

# Collect all DOS

dos_all = zeros(length(dos0[:,1]),numDat)

for i in 1:numDat
	dos_n = readdlm("data$(i)/$(filname).dat")
	dos_all[:,i] = dos_n[:,2]
end

for i in 1:length(dos_all[:,1])
	outfil[i,2] = 1/numDat * sum(dos_all[i,:])
	outfil[i,3] = sqrt(1/(numDat-1) * sum((dos_all[i,:] .- outfil[i,2]) .^ 2))
end

writedlm(outfilname,outfil)
