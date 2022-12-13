using LinearAlgebra
using PyPlot
plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3
rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)


function VDot(a)
    sx = [0 1; 1 0.]
    sy = [0 -1im; 1im 0.]
    sz = [1 0; 0 -1.] 
    return a[1]*sx + a[2]*sy + a[3]*sz
end
function absv(a)
    return sqrt(dot(a,a))
end
function UnitVector(a,c,j,num_laps,N)
    phij = 2*pi*(j-1)*num_laps/(N-1)
    cj = c*(j-1)/(N-1)
    return [a*cos(phij),a*sin(phij),cj]
end


function HeffSpinless(N,t0,tph,w0,X,gL,gR)
	H = zeros(N,N) * (1+0.0*1im)
	H[1,1] = 1im/2 * gL
	H[end,end] = 1im/2 * gR
	for i in 1:N-1
		H[i,i+1] = -t0-tph * sqrt(2*w0) * X[i]
		H[i+1,i] = -t0-tph * sqrt(2*w0) * X[i]
	end
	return H
end

function Heff(N,t0,tso,tph,tsoph,gL,gR,w0,X)
	num_laps = 1
	num_ions = N
	a = 1.0
	c = 1.0
	Hs = zeros(2*N,2*N) * (1+ 0.0*1im)
	for i in 1:2*N-2
		Hs[i,i+2] = -t0
	end
	k = 1
	for j in 1:2:2*N-4
		rj = UnitVector(a,c,k,num_laps,N)
		rj1 = UnitVector(a,c,k+1,num_laps,N)
		rj2 = UnitVector(a,c,k+2,num_laps,N) 
		rjrj1 = rj-rj1
		rjrj2 = rj-rj2

		dj1 = rjrj1/absv(rjrj1)
		dj2 = rjrj2/absv(rjrj2)
		
		vj = cross(dj1,dj2)
		 
		Vj = 1im*tso*VDot(vj)

		Hs[j,j+4] = Vj[1,1]
		Hs[j+1,j+4] = Vj[2,1]
		Hs[j,j+5] = Vj[1,2]
		Hs[j+1,j+5] = Vj[2,2]  
		k += 1
	end
	k = 1
	for i in 1:2*N-2
       		Hs[i,i+2] = -tph * sqrt(2*w0)*X[k]  
        	if i % 2 == 0
	            k += 1
        	end
	end
	k = 1
	for j in 1:2:2*N-4
		n = Int((j+1)/2)
		rj = UnitVector(a,c,k,num_laps,N)
		rj1 = UnitVector(a,c,k+1,num_laps,N)
		rj2 = UnitVector(a,c,k+2,num_laps,N) 
		rjrj1 = rj-rj1
		rjrj2 = rj-rj2

		dj1 = rjrj1/absv(rjrj1)
		dj2 = rjrj2/absv(rjrj2)
		
		vj = cross(dj1,dj2) 
		Vj = 1im*tsoph*VDot(vj) * sqrt(2*w0) * X[k]

		Hs[j,j+4] = Vj[1,1]
		Hs[j+1,j+4] = Vj[2,1]
		Hs[j,j+5] = Vj[1,2]
		Hs[j+1,j+5] = Vj[2,2]  
		k += 1
	end
	HS = Hs + Hs'
	HS[1,1] += -1im/2*gL[1]
	HS[2,2] += -1im/2*gL[2]
	HS[end-1,end-1] += -1im/2*gR[1]
        HS[end,end] += -1im/2 * gR[2]
	return HS
end

function Liouvillian(H,rho)
	return -1im * (H * rho - rho * H)
end


function LiouvillianEigsXVar(N,P,Xs,X2)
	gL = 0.25
	gR = 0.25
	GL = gL * (ones(2) + P * [1,-1])
	GR = gR * ones(2)
	w0 = 0.05
	t0 = 1.0
	tso = t0/40
	tph = 0.1*t0
	tsoph = 0.1*tso
	#Id = Matrix(I,N,N)
	Id = Matrix(I,2*N,2*N)
	ew = []

	for X in Xs
		H = Heff(N,t0,tso,tph,tsoph,GL,GR,w0,[X,X2])
	#	H = HeffSpinless(N,t0,tph,w0,[X],gL,gR)
		L = -1im*(kron(Id,H) - kron(transpose(H),Id))
	#	Lr = L[1:end-1,1:end-1]
		#if X == Xs[1]
	#		println(eigvals(Lr))
	#	end
		push!(ew,eigvals(L))
	end
	return ew
end


fig = plt.figure(figsize = (10,6))
px = 0.08
py = 0.56
sx = 0.445
dx = 0.02
sy = 0.43
dy = 0.04

ax11 = plt.axes([px,py,sx,sy])
ax12 = plt.axes([px,py-sy-dy,sx,sy])

ax21 = plt.axes([px+dx+sx,py,sx,sy])
ax22 = plt.axes([px+dx+sx,py-sy-dy,sx,sy])

#ax31 = plt.axes([px+2*(dx+sx),py,sx,sy])
#ax32 = plt.axes([px+2*(dx+sx),py-sy-dy,sx,sy])

ax = [[ax11,ax12],[ax21,ax22]]
Ps = [0.0,0.5]
Pstr = [L"$0.0$",L"$0.5$"]
astr = [["(a)","(b)"],["(c)","(d)"]]
lw = 2
Xs = Array(0.0:0.0025:5.0)
N = 3
col = get_cmap(:brg)
X2 = 1.0

for i in 1:length(ax)
	ew = LiouvillianEigsXVar(N,Ps[i],Xs,X2)
	if i == 1
		println(ew[1])
	end
	inds_r = [sortperm(real.(e)) for e in ew]
	inds_i = [sortperm(imag.(e)) for e in ew]
	ew_s_re = [real.(e[i]) for (e,i) in zip(ew,inds_r)]
	ew_s_im = [imag.(e[i]) for (e,i) in zip(ew,inds_i)]
	for k in 1:length(ew_s_re[1])
		ax[i][1].plot(Xs,[ew_s_re[j][k]/0.25 for j in 1:length(ew_s_re)],linewidth=lw,color=col(k/(4*N^2)))
		ax[i][2].plot(Xs,[ew_s_im[j][k]/0.25 for j in 1:length(ew_s_re)],linewidth=lw,color=col(k/(4*N^2)))
	end
	ax[i][1].set_ylim(-0.9,0.9)
	ax[i][1].set_xticklabels([])
	ax[i][1].annotate(L"$P=$"*Pstr[i],(0.4,0.87),xycoords="axes fraction",fontsize=fs-2)
	ax[i][2].set_ylim(-0.8,0.8)
	ax[i][2].set_xlabel(L"$X_1$")
	if i > 1
		for j in 1:2
			ax[i][j].set_yticklabels([])
		end
	end
	for j in 1:2
                ax[i][j].annotate(astr[i][j],(0.02,0.92),xycoords="axes fraction",fontsize=fs-2)
        end	
end
ax[1][1].set_ylabel(L"$\mathrm{Re}(\lambda_n) \, [\Gamma_0]$")
ax[1][2].set_ylabel(L"$\mathrm{Im}(\lambda_n) \, [\Gamma_0]$")



plt.savefig("Liouvillian_Spec_X2_$(X2)_X1Var.pdf")
