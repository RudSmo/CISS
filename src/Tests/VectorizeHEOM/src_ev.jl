using LinearAlgebra
using PyPlot


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
	HS[1,1] += -1im/2*gL
	HS[2,2] += -1im/2*gL
	HS[end-1,end-1] += -1im/2*gR
        HS[end,end] += -1im/2 * gR
	return HS
end

function Liouvillian(H,rho)
	return -1im * (H * rho - rho * H)
end

N = 5
gL = 0.25
gR = 0.25
w0 = 0.05
t0 = 1.0
tso = t0/40
tph = 0.1*t0
tsoph = 0.1*tso
Id = Matrix(I,2*N,2*N)

Xs = [5.0] 
#Array(0.0:0.01:6)
ev = []


for X in Xs
	H = Heff(N,t0,tso,tph,tsoph,gL,gR,w0,[X,1.0,-0.2,1.0])
	display(H)
	L = -1im*(kron(transpose(H),Id) - kron(Id,H))
	EV = eigvecs(L)
	EV_ol = zeros((2*N)^2,(2*N)^2)
	for i in 1:length(EV[1,:])
		for j in 1:length(EV[1,:])
			EV_ol[i,j] = real(sqrt((EV[:,i]' * EV[:,j])) / (sqrt((EV[:,i]' * EV[:,i]) * (EV[:,j]' * EV[:,j]))))
		end
	end
	
	push!(ev,EV_ol)
end

plt.imshow(ev[1],aspect="auto")

#fig,ax = plt.subplots(2,1)

#ax[1].plot(Xs,real.(ew))
#ax[2].plot(Xs,imag.(ew))
plt.show()
