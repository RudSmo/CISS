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

fac_c = 0.2*0.002*1e3 * 0.2433946
fac_t = 1/(0.2*0.002*4.135668*1e3)

phonons = true
electrons = false
currents = false
energy = false

if phonons
	dat1 = readdlm("L3x6/L18_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_phonons.dat")
	dat2 = readdlm("L3x6/L18_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P-0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_phonons.dat")

	n1 = length(dat1[2:end,1])
	n2 = length(dat2[2:end,1])
	nf = min(n1,n2)
	X1 = dat1[2:nf,2:2:end]
	P1 = dat1[2:nf,3:2:end]
	X2 = dat2[2:nf,2:2:end]
	P2 = dat2[2:nf,3:2:end]
	N = length(X1[1,:])-1
	tim = dat1[2:nf,1]*fac_t
	
	lyapunov=[]
	for i in 1:N
		dist_i = sqrt.((X1[:,i]-X2[:,i]) .^ 2 .+ (P1[:,i]-P2[:,i]) .^ 2)
		dist_0 = sqrt.((X1[1,i]-X2[1,i]) .^ 2 .+ (P1[1,i]-P2[1,i]) .^ 2)
		log_dist = log.(dist_i/dist_0)
		push!(lyapunov,log_dist)
	end

elseif electrons

	dat1 = readdlm("unpert/L5_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_electrons.dat")
	dat2 = readdlm("pert/L5_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_electrons.dat")
	n1 = length(dat1[2:end,1])
        n2 = length(dat2[2:end,1])
        nf = min(n1,n2)	
	pop1 = dat1[2:nf,2:2:end] .+ dat1[2:nf,3:2:end]
	pop2 = dat2[2:nf,2:2:end] .+ dat2[2:nf,3:2:end]
	N = length(pop1[1,:])
	tim = dat1[2:nf,1]*fac_t

	for i in 1:N
		dist_i = sqrt.((pop1[:,i]-pop2[:,i]) .^ 2 .+ (pop1[:,i]-pop2[:,i]) .^ 2)
		dist_0 = sqrt.((pop1[1,i]-pop2[1,i]) .^ 2 .+ (pop1[1,i]-pop2[1,i]) .^ 2)
		log_dist = log.(dist_i/dist_0)
		push!(lyapunov,log_dist)
	end
elseif currents
	dat1 = readdlm("unpert/L5_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_currents.dat")
        dat2 = readdlm("pert/L5_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_currents.dat")
        n1 = length(dat1[2:end,1])
        n2 = length(dat2[2:end,1])
        nf = min(n1,n2)
        I1 = dat1[2:nf,2] .- dat1[2:nf,3]
        I2 = dat2[2:nf,2] .- dat2[2:nf,3]
        tim = dat1[2:nf,1]*fac_t
	lyapunov = []
        disti = sqrt.((I1-I2) .^ 2 .+ (I1-I2) .^ 2)
	dist0 = sqrt.((I1[1]-I2[1]) .^ 2 .+ (I1[1]-I2[1]) .^ 2)
	logDist = log.(disti/dist0)
	push!(lyapunov,logDist)
elseif energy
	dat1 = readdlm("unpert/L5_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_convergence_info.dat")
        dat2 = readdlm("pert/L5_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_convergence_info.dat")
        n1 = length(dat1[2:end,1])
        n2 = length(dat2[2:end,1])
        nf = min(n1,n2)
        I1 = dat1[2:nf,2] .- dat1[2:nf,3]
        I2 = dat2[2:nf,2] .- dat2[2:nf,3]
        tim = dat1[2:nf,1]*fac_t
        lyapunov = []
        disti = sqrt.((I1-I2) .^ 2 .+ (I1-I2) .^ 2)
        dist0 = sqrt.((I1[1]-I2[1]) .^ 2 .+ (I1[1]-I2[1]) .^ 2)
        logDist = log.(disti/dist0)
        push!(lyapunov,logDist)	

end

px = 0.16
py = 0.15
sx = 0.81
sy = 0.82
dy = 0.05
dx = 0.01

fig = plt.figure(figsize=(7,4))
ax = plt.axes([px,py,sx,sy])
colors = [get_cmap(:plasma)(i) for i in 0:1/length(lyapunov):1]
 
for i in 1:length(lyapunov)
	if i == 1
		lab = "mode "*L"$j=\ $"*"$(i)"
	else
		lab = L"$j=\ $"*"$(i)"
	end
	ax.plot(tim,lyapunov[i] ./ tim,linewidth=lw,label=lab,color=colors[i])
end
ax.set_ylabel("Lyapunov exponent "*L"$\lambda_j \, \left[\frac{1}{\mathrm{ps}}\right]$") 
ax.plot(tim,10 ./ tim, color = "red", linewidth=lw)
##=\frac{1}{t}\mathrm{log}\left(\frac{d_j(t)}{d_j(0)}\right)\, \left[\frac{1}{\mathrm{ps}}\right]$")

plt.legend(loc="upper right",ncol=2)
#ax.set_xlim(0,10)
ax.set_ylim(-0.05,10)
ax.set_xlabel("time "*L"$t\,[\mathrm{ps}]$")
plt.savefig("L3x6_V60_lyapunov_phonons_PM.pdf")




