using LinearAlgebra
using PyPlot


#eta = []
#gLRs = Array(0.0:0.5:100.0)
#for gLR in gLRs

function TB_SSH(L,l,gLR,tL,eL,t1,t2,eR)
	H = zeros(L,L);


	# tb chain
	for i in 1:l
		if i == l
			k = gLR
		else
			k = tL
		end

		H[i,i] = eL
		H[i,i+1] = k
		H[i+1,i] = k
	end	

	# ssh model
	for i in l+1:2:L-1
		H[i,i+1] = -t1
		H[i+1,i] = -t1
		if i < L-2
			H[i+1,i+2] = -t2
			H[i+2,i+1] = -t2
		end
	end
	return H
end

L = 1000;
l = 500
gLR = 1.0
tL = 10.0
t1 = 2.0
t2 = 1.0
eL = 0.0
eR = 0.0 
H = TB_SSH(L,l,gLR,tL,eL,t1,t2,eR)

#display(H)
ev = eigvecs(H)
#display(ev)

HL = H[1:l,1:l]
HR = H[l+1:L,l+1:L]

#print(size(HL))
#print(size(HR))

el,vl = eigen(HL)
er,vr = eigen(HR)


TH = zeros(L,L)
for i in 1:l
	TH[i,i] = el[i]
end
for i in l+1:L
	TH[i,i] = er[i-l]
end
#display(TH)

Cl = H[1:l,l+1:L]
Cr = H[l+1:L,1:l]

# Compute unitary for coupling
#Ue = ev[1:l,1:l]
#UCl = ev[1:l,l+1:L] 
#UCr = ev[l+1:L,1:l]
#Us  = ev[l+1:L,l+1:L]

#size(HL) = (l,l)
#size(HR) = (L-l-1,L-l-1)
#size(Cl) = (l,L-l-1)
#size(Cr) = (L-l-1,l)

# Plot block diagonalization elements
#A = Ue' * HL  * UCl  + UCl * Cl' * UCl + UCl * HR * Us  
#B = Ue * Cl * Us
#
#fig = plt.figure(figsize=(12,8))
#ax1 = plt.axes([0.1,0.1,0.35,0.8])
#ax2 = plt.axes([0.55,0.1,0.35,0.8])
#ax1.annotate("(a)",xy=(0.1,0.93),xycoords="axes fraction",fontsize=18)
#ax2.annotate("(b)",xy=(0.1,0.93),xycoords="axes fraction",fontsize=18)
#
#vm = 0.005
#ax1.imshow(A,aspect="auto",vmin=-vm,vmax=vm)
#im = ax2.imshow(B,aspect="auto",vmin=-vm,vmax=vm)
#fig.colorbar(im,ax=ax2)

#display(Cl)
#display(Cr)
#display(HR)
#print(eigvals(HR))

Tl = vl*Cl*vr'
Tr = vr*Cr*vl'

for i in 1:l
	for j in l+1:L
		TH[i,j] = Tl[i,j-l]
		TH[j,i] = Tr[j-l,i]
	end
end

#display(TH)
#print("\n")
# Test if spectra of blocked Hamiltonian and full one are the same

eH = eigvals(H)
#eTH = eigvals(TH)


#plt.plot(er)
#plt.plot(eH,label="exact")
#plt.plot(eTH,label="block")
#plt.legend()

#println(eH)
#println(eTH)

#print(eigvals(TH)≈eigvals(H))
#print(sum(eH-eTH)/length(eH))

#push!(eta,sum(abs.(abs.(eH).-abs.(eTH)))/length(eH))

# Plot block Hamiltonian

#display(TH)
#plt.imshow(TH)


# Plot spectrum of left and right parts

#plt.scatter(ones(length(el)),el,s=3)
#plt.scatter(1.1*ones(length(er)),er,s=3)
#plt.scatter(ones(length(eH)),eH,s=3)
plt.plot(eH)
plt.plot(-eH)

# Plot spectrum in k-space



#LR = length(er)
#plt.plot(er[1:Int(LR/2)],color="k")
#plt.plot(Array(0:-1:-Int(LR/2)+1),er[1:Int(LR/2)],color="k")
#
#plt.plot(er[end:-1:Int(LR/2)+1],color="k")
#plt.plot(Array(0:-1:-Int(LR/2)+1),er[end:-1:Int(LR/2)+1],color="k")
#
#LL = length(el)
#plt.plot(el[1:Int(LL/2)],color="red")
#plt.plot(Array(0:-1:-Int(LL/2)+1),el[1:Int(LL/2)],color="red")
#
#plt.plot(el[end:-1:Int(LL/2)+1],color="red")
#plt.plot(Array(0:-1:-Int(LL/2)+1),el[end:-1:Int(LL/2)+1],color="red")


# Plot coupling of left mode k to all right modes # Modes are stored in ascending format

k = Int(floor(l/2))
#plt.plot(er,TH[k,l+1:L],label=L"$\varepsilon_\mathrm{env} = 0.00$")
#plt.plot(er,TH[Int(floor(l/4)),l+1:L],label=L"$\varepsilon_\mathrm{env} = -14.23$")
#plt.plot(er,TH[1,l+1:L],label=L"$\varepsilon_\mathrm{env} = -20.00$")
#plt.xlabel(L"$\mathrm{System} \ \mathrm{mode} \ \mathrm{energy} \  \varepsilon_\mathrm{sys}$")
#plt.ylabel(L"$\mathrm{Coupling} \ A(\varepsilon_\mathrm{env},\varepsilon_\mathrm{sys})$ ")
#plt.legend()

# Plot coupling of right mode k to all left modes # Modes are stored in ascending format

#print(er[1])
#print(er[3])
#print(er[5])
#
#plt.plot(el,TH[L-9+1,1:l],label=L"$\varepsilon_\mathrm{sys} = -1.9$")
#plt.plot(el,TH[L-9+3,1:l],label=L"$\varepsilon_\mathrm{sys} = -1.17$")
#plt.plot(el,TH[L-9+5,1:l],label=L"$\varepsilon_\mathrm{sys} = 0.00$")
#plt.xlabel(L"$\mathrm{Environment} \ \mathrm{mode} \ \mathrm{energy} \  \varepsilon_\mathrm{sys}$")
#plt.ylabel(L"$\mathrm{Coupling} \ A(\varepsilon_\mathrm{sys},\varepsilon_\mathrm{env})$ ")
#plt.legend()

#print(er)

# Plot coupling of all left modes to all right modes
#print(er)
#fig = plt.figure(figsize=(12,8))
#ax = plt.axes([0.1,0.1,0.85,0.85])
#im = ax.imshow(TH[1:l,l+1:L],aspect="auto",extent=[er[1],er[end],el[1],el[end]],interpolation="none")
#ax.set_xlabel(L"$\mathrm{System} \ \mathrm{mode} \ \mathrm{energy} \  \varepsilon_\mathrm{sys}$")
#ax.set_ylabel(L"$\mathrm{Environment} \ \mathrm{mode} \ \mathrm{energy} \  \varepsilon_\mathrm{env}$")
#fig.colorbar(im,ax=ax,label=L"A(\varepsilon_\mathrm{env},\varepsilon_\mathrm{sys})")
#plt.matshow(TH[l+1:L,1:l],aspect="auto")

#end
#
#
#plt.plot(gLRs,eta)
#plt.xlabel(L"$\gamma$")
#plt.ylabel(L"$\eta$")


# Coupling as defined in NEGF

#Gamma = zeros(L-l,L-l)
#Gamma[1,1] = 2*pi*gLR
#vr = eigvecs(HR)
#GT = vr' * Gamma * vr
##plt.imshow(GT,aspect="auto")
#plt.plot(GT[5,:])
plt.show()


