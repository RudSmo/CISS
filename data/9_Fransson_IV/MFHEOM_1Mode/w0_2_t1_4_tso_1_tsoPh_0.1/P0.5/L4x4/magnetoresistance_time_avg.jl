using DelimitedFiles

mag_res = []
mag_res_std = []
n = 4
m = 4
P = 0.5
datP = readdlm("L$(n)x$(m)_currents_time_P0.5_avg.dat") 
datP_std = readdlm("L$(n)x$(m)_currents_time_P0.5_std.dat")
datM = readdlm("../../P-0.5/L$(n)x$(m)/L$(n)x$(m)_currents_time_P-0.5_avg.dat")
datM_std = readdlm("../../P-0.5/L$(n)x$(m)/L$(n)x$(m)_currents_time_P-0.5_std.dat")

Vs = datP[:,1]
curP = datP[:,2] - datP[:,3]
stdPL = datP_std[:,2]
stdPR = datP_std[:,3] 

curM = datM[:,2] - datM[:,3]
stdML = datM_std[:,2]
stdMR = datM_std[:,3]

	
curD = [curP[i] - curM[i] for i in 1:length(curP)]
curT = [curP[i] + curM[i] for i in 1:length(curP)]
mr_avg = [curD[i] / curT[i] * 100 for i in 1:length(curD)]
dPL = [(curT[i] - curD[i]) / curT[i]^2 for i in 1:length(curT)]
dPR = -dPL
dML = [-(curT[i] + curD[i]) / curT[i]^2 for i in 1:length(curT)]
dMR = -dML

mr_std = [100 * sqrt( (dPL[i] * stdPL[i]) ^ 2 + (dPR[i] * stdPR[i]) ^ 2 + (dML[i] * stdML[i]) ^ 2 + (dMR[i] * stdMR[i])^2 ) for i in 1:length(curD)]

fil = open("L$(n)x$(m)_magnetoresistance_avg_std.dat","w")
for i in 1:length(Vs)
	write(fil,"$(Vs[i]) $(mr_avg[i]) $(mr_std[i])\n")
end
close(fil)
