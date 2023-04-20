using DelimitedFiles

Vs = [48,55,65]
mag_res = []
mag_res_std = []
for v in Vs
	datP = readdlm("V$(v)/currents_time_avg.dat")
	datM = readdlm("../../P-0.5/L4x4/V$(v)/currents_time_avg.dat")
	
	curP = datP[1,1] - datP[1,2]
	stdPL,stdPR = datP[2,1:2] 
	curM = datM[1,1] - datM[1,2]
	stdML,stdMR = datM[2,1:2]
	
	curD = curP - curM
	curT = curP + curM
	mr_avg = curD / curT
	
	dPL = (curT - curD) / curT^2
	dPR = -dPL
	dML = -(curT+curD) / curT^2
	dMR = -dML
	mr_std = sqrt( (dPL * stdPL)^2 + (dPR * stdPR)^2 + (dML * stdML)^2 + (dMR * stdMR)^2 ) 	
	push!(mag_res,mr_avg)
	push!(mag_res_std,mr_std)
end


fil = open("Magnetoresistance_time_avg.dat","w")
for i in 1:length(Vs)
	write(fil,"$(Vs[i]) $(mag_res[i]) $(mag_res_std[i])\n")
end
close(fil)
