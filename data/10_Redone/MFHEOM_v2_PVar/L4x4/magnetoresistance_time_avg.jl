using DelimitedFiles

Vs = [50,60,75]
Ns = [11,11,11]
tA = 20000
mag_res = []
mag_res_std = []
dat0 = readdlm("V$(Vs[1])/currents_time_avg_nS$(Ns[1])_P.dat")
nt = length(dat0[tA:end,1])
for (v,nS) in zip(Vs,Ns)
	datP = readdlm("V$(v)/currents_time_avg_nS$(nS)_P_std.dat")
	datM = readdlm("V$(v)/currents_time_avg_nS$(nS)_M_std.dat")
	
	curP = datP[1,1] - datP[1,2]
	stdPL,stdPR = datP[3,1:2] 
	curM = datM[1,1] - datM[1,2]
	stdML,stdMR = datM[3,1:2]
	
	curD = curP - curM
	curT = curP + curM
	mr_avg = curD / curT * 100
	
	dPL = (curT - curD) / curT^2
	dPR = -dPL
	dML = -(curT+curD) / curT^2
	dMR = -dML
	mr_std = 100 * sqrt( (dPL * stdPL)^2 + (dPR * stdPR)^2 + (dML * stdML)^2 + (dMR * stdMR)^2 ) #/ sqrt(nt)	
	push!(mag_res,mr_avg)
	push!(mag_res_std,mr_std)
end


fil = open("Magnetoresistance_time_avg_stderr.dat","w")
for i in 1:length(Vs)
	write(fil,"$(Vs[i]) $(mag_res[i]) $(mag_res_std[i])\n")
end
close(fil)
