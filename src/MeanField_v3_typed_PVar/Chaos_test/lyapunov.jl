using PyPlot
using DelimitedFiles

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 19
lw = 3
fac_c = 0.002*0.2*1e3 * 0.2433946
fac_t = 1/(0.002*0.2*4.135668*1e3)

phonons = true
electrons = false
currents = false
energy = false

if phonons
	dat1 = readdlm("unpert/L8_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_phonons.dat")
	dat2 = readdlm("pert/L8_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_phonons.dat")

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
		log_dist = log2.(dist_i/dist_0)
		push!(lyapunov,log_dist)
	end

elseif electrons

	dat1 = readdlm("unpert/L8_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_electrons.dat")
	dat2 = readdlm("pert/L8_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_electrons.dat")
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
		log_dist = log2.(dist_i/dist_0)
		push!(lyapunov,log_dist)
	end
elseif currents
	dat1 = readdlm("unpert/L8_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_currents.dat")
        dat2 = readdlm("pert/L8_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_currents.dat")
        n1 = length(dat1[2:end,1])
        n2 = length(dat2[2:end,1])
        nf = min(n1,n2)
        I1 = dat1[2:nf,2] .- dat1[2:nf,3]
        I2 = dat2[2:nf,2] .- dat2[2:nf,3]
        tim = dat1[2:nf,1]*fac_t
	lyapunov = []
        disti = sqrt.((I1-I2) .^ 2 .+ (I1-I2) .^ 2)
	dist0 = sqrt.((I1[1]-I2[1]) .^ 2 .+ (I1[1]-I2[1]) .^ 2)
	logDist = log2.(disti/dist0)
	push!(lyapunov,logDist)
elseif energy
	dat1 = readdlm("unpert/L8_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_convergence_info.dat")
        dat2 = readdlm("pert/L8_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_convergence_info.dat")
        n1 = length(dat1[2:end,1])
        n2 = length(dat2[2:end,1])
        nf = min(n1,n2)
        I1 = dat1[2:nf,2] .- dat1[2:nf,3]
        I2 = dat2[2:nf,2] .- dat2[2:nf,3]
        tim = dat1[2:nf,1]*fac_t
        lyapunov = []
        disti = sqrt.((I1-I2) .^ 2 .+ (I1-I2) .^ 2)
        dist0 = sqrt.((I1[1]-I2[1]) .^ 2 .+ (I1[1]-I2[1]) .^ 2)
        logDist = log2.(disti/dist0)
        push!(lyapunov,logDist)	

end

px = 0.12
py = 0.12
sx = 0.84
sy = 0.84
dy = 0.05
dx = 0.01

fig = plt.figure(figsize=(7,5))
ax = plt.axes([px,py,sx,sy])

for i in 1:length(lyapunov)
	ax.plot(tim,lyapunov[i] ./ tim,linewidth=lw,label="mode "*L"$j=$"*"$(i)")
end
ax.set_ylabel("Lyapunov exponent"*L"$\lambda_j$")
plt.legend(loc="lower right",ncol=2)
ax.set_xlim(0,120)
ax.set_xlabel("time "*L"$t\,[\mathrm{ps}]$")
plt.savefig("L2x4_V60_lyapunov_phonons.pdf")




