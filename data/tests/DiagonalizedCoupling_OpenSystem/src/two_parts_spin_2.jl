using LinearAlgebra
using PyPlot

function vdot(S)
	sx = [0 1; 1 0.]
	sy = [0 -1im; 1im 0]
	sz = [1 0; 0 -1.]
	return S[1] * sx + S[2] * sy + S[3] * sz
end

case = "normal"
Nl = 200
tl = 1.0
eps_l = 2*ones(2)
kL = 1.0*ones(2)
i0 = Nl+1
eps = eps_l
kR = kL
Nr = Nl
tr = tl
eps_r = eps_l

Nt = (Nl + Nr + 1)*2


# Left electrode

HL = zeros(Nt,Nt)

for i in 1:Nl
	k = 2*i-1
	if i < Nl
		kp = 2*(i+1)-1

		HL[k,kp] = -tl
		HL[kp,k] = -tl
		HL[k+1,kp+1] = -tl
		HL[kp+1,k+1] = -tl
	end

	if case == "gapped"
		HL[k,k] = (-1)^i * eps_l[1]
		HL[k+1,k+1] = (-1)^i * eps_l[2]
	else
		HL[k,k] = eps_l[1]
		HL[k+1,k+1] = eps_l[2]
	end
end

# Coupling Hamiltonian

VLC = zeros(Nt,Nt)
VCL = zeros(Nt,Nt)

ki0m = 2*(i0-1)-1
ki0 = 2*i0-1
VLC[ki0m,ki0] = kL[1]
VLC[ki0m+1,ki0+1] = kL[2]

VCL[ki0,ki0m] = kL[1]
VCL[ki0+1,ki0m+1] = kL[2]



# Central Hamiltonian

HC = zeros(Nt,Nt) * (1+0.0*1im)
S = [0.0,0.3,1.0]
Jsd = 4.5
Sv = Jsd*vdot(S)

if case == "gapped"
	HC[ki0,ki0] = Sv[1,1] + (-1)^(i0)*eps[1]
	HC[ki0+1,ki0+1] = Sv[2,2] + (-1)^(i0)*eps[2]
else
	HC[ki0,ki0] = Sv[1,1] + eps[1]
	HC[ki0+1,ki0+1] = Sv[2,2] + eps[2]
end
HC[ki0,ki0+1] = Sv[1,2]


HC[ki0+1,ki0] = Sv[2,1]




# Right Coupling

VRC = zeros(Nt,Nt)
VCR = zeros(Nt,Nt)
ki0p = 2*(i0+1)-1
VCR[ki0,ki0p] = kR[1]
VCR[ki0+1,ki0p+1] = kR[2]
VRC[ki0p,ki0] = kR[1]
VRC[ki0p+1,ki0+1] = kR[2]

# Right lead

HR = zeros(Nt,Nt)
for i in i0+1:(Nl + Nr + 1)
	k = 2*i-1
        if i < (Nl + Nr + 1)
                kp = 2*(i+1)-1

                HR[k,kp] = -tr
                HR[kp,k] = -tr
		HR[k+1,kp+1] = -tr
		HR[kp+1,k+1] = -tr
        end
	if case == "gapped"
	        HR[k,k] = (-1)^i * eps_r[1]
        	HR[k+1,k+1] = (-1)^i * eps_r[2]
	else
		HR[k,k] = eps_r[1]
		HR[k+1,k+1] = eps_r[2]
	end
end	


HT = HL + VLC + VCL + HC + VCR + VRC + HR


# Compute self energy of lead
# Sigma_l = VCl * Gl * VlC

function Gret(H,E,eps=1e-2)
	idm = Matrix(I,size(H))
	return inv((E+eps*1im)*idm-H)
end


# PLOT DOS
ew,ev = eigen(HL[1:2*Nl,1:2*Nl])
rng = -1:1
rng1 = ones(length(rng))
#index = 3
#psi = [ev[index,j:j+1]' * ev[index,j:j+1] for j in 1:2:length(ev[1,:])]
#plt.plot(psi)

for (i,e) in enumerate(ew)
	plt.plot(rng,e*rng1,color="k")
	#if e == 0.0	
	#	index = i
	#	psi = [ev[index,j:j+1]' * ev[index,j:j+1] for j in 1:2:length(ev[1,:])]
	#	plt.plot(psi)
	#end
end
plt.ylabel("Energy")
plt.savefig("DOS_HL_$(case).pdf")

# PLOT RETARDED GF AND SELFENERGY
Sig_E = []
GL_E = []
GL_E_00 = []
Erng = -5.0:0.1:7.0

for E in Erng
	GL = Gret(HL,E)
	push!(GL_E,sum([GL[l,l] for l in 1:size(GL)[1]]))
	push!(GL_E_00,sum(GL[l,l] for l in [2*(i0-1)-1,2*(i0-1)]))
	Sigma_L = VCL * GL * VLC
	GR = Gret(HR,E)
	Sigma_R = VCR * GR * VRC
	Sigma = Sigma_L + Sigma_R
	push!(Sig_E,sum([Sigma[l,l] for l in 1:size(Sigma)[1]]))
end
plt.clf()
plt.plot(Erng,real(Sig_E),color="k")
plt.plot(Erng,imag.(Sig_E),color="purple")
plt.ylabel("Self energy "*L"$\Sigma$")
plt.xlabel(L"$E$")
plt.savefig("SelfEnergy_HL_$(case).pdf")


plt.clf()
plt.plot(Erng,real.(GL_E)/100,color="k")
plt.plot(Erng,imag.(GL_E)/100,color="purple")
plt.ylabel("Retarded Greens function "*L"$G_L^R$")
plt.xlabel(L"$E$")
plt.savefig("GreenFunc_HL_$(case).pdf")


plt.clf()
plt.plot(Erng,real.(GL_E_00)/100,color="k")
plt.plot(Erng,imag.(GL_E_00)/100,color="purple")
plt.ylabel("Edge block retarded Greens function "*L"$G_{L,\mathrm{edge}}^R$")
plt.xlabel(L"$E$")
plt.savefig("GreenFunc00_HL_$(case).pdf")

#plt.show()




