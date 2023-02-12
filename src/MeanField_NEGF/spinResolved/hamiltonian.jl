module Hamiltonian
using LinearAlgebra

export VDot,absv,UnitVector,FranssonHamiltonian,SingleLevel,HatanoNelson,HatanoNelson2D,NNN,FranssonHamiltonianExtended


function VDot(a)
    sx = [0 1; 1 0.]
    sy = [0 -1im; 1im 0.]
    sz = [1 0; 0 -1.] 
    return a[1]*sx + a[2]*sy + a[3]*sz
end
function absv(a)
    return sqrt(dot(a,a))
end
#-------------------------- Electronic Hamiltonian-------------------------------------
function UnitVector(a,c,j,num_laps,N)
    phij = 2*pi*(j-1)*num_laps/(N-1)
    cj = c*(j-1)/(N-1)
    return [a*cos(phij),a*sin(phij),cj]
end
function SOC_A_nn(N,tso,a,c,num_laps)
    Hs = zeros(2*N,2*N)*(1+0.0*1im)
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
    return Hs+Hs'
end

function Hopping_A_nn(N,t)
    Hh = zeros(2*N,2*N)
    for i in 1:2*N-2
        Hh[i,i+2] = -t
        Hh[i+2,i] = -t
    end
    return Hh
end

function Hopping_A_nn_spin(N,td,to)
    tu,td = td
    tud,tdu = to
    Hh = zeros(2*N,2*N)
    for i in 1:2:2*N-2
        Hh[i,i+2] = -tu
        Hh[i+2,i] = -tu
	Hh[i+1,i+3] = -td
	Hh[i+3,i+1] = -td
	Hh[i,i+3] = -tud
	Hh[i+3,i] = -tud
	Hh[i+1,i+2] = -tdu
	Hh[i+2,i+1] = -tdu
    end
    return Hh
end

function Onsite_A(N,E)
    eu,ed = E 
    H = zeros(2*N,2*N)
    for i in 1:2:2*N-1
        H[i,i] = eu 
    end
    for i in 2:2:2*N
        H[i,i] = ed 
    end
    return H
end


#-------------------------- Electron-Phonon Hamiltonian-------------------------------------
function Holstein_NN(N,tph,ws,qs)
    H = zeros(2*N,2*N)
    k = 1
    for i in 1:2*N-2
        H[i,i+2] = -tph * sqrt(2*ws[k])*qs[k]  
        H[i+2,i] = -tph * sqrt(2*ws[k])*qs[k]  
        if i % 2 == 0
            k += 1
        end
    end
    return H
end

function Holstein_SOC(N,tsoph,a,c,num_laps,ws,qs)
    Hs = zeros(2*N,2*N)*(1+0.0*1im)
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
        Vj = 1im*tsoph*VDot(vj) * sqrt(2*ws[n]) * qs[n]

        Hs[j,j+4] = Vj[1,1]
        Hs[j+1,j+4] = Vj[2,1]
        Hs[j,j+5] = Vj[1,2]
        Hs[j+1,j+5] = Vj[2,2]  
        k += 1
    end
    return Hs+Hs'
end



function FranssonHamiltonian(N,E,t,tso,a,c,num_laps,tph,tsoph,ws,qs)
	Ho = Onsite_A(N,E)
	Hh = Hopping_A_nn(N,t)
	HSo = SOC_A_nn(N,tso,a,c,num_laps)
    Heph = Holstein_NN(N,tph,ws,qs)
    Hsoph = Holstein_SOC(N,tsoph,a,c,num_laps,ws,qs)
	return Ho+Hh+HSo+Heph+Hsoph
end


function FranssonHamiltonianExtended(N,E,t,gSL,tso,a,c,num_laps,tph,tsoph,ws,qs)
	N1,NS,N2 = N
	Ntot = N1+NS+N2
	E1,ES,E2 = E
	t1,tS,t2 = t
	gS1,gS2 = gSL
	tso1,tsoS,tso2 = tso
	a1,aS,a2 = a
	c1,cS,c2 = c
	nL1,nLS,nL2 = num_laps
	tph1,tphS,tph2 = tph
	tsoph1,tsophS,tsoph2 = tsoph
	ws1,wsS,ws2 = ws
	qs1,qsS,qs2 = qs
	H1 = FranssonHamiltonian(N1,E1,t1,tso1,a1,c1,nL1,tph1,tsoph1,ws1,qs1)
	HS = FranssonHamiltonian(NS,ES,tS,tsoS,aS,cS,nLS,tphS,tsophS,wsS,qsS)
	H2 = FranssonHamiltonian(N2,E2,t2,tso2,a2,c2,nL2,tph2,tsoph2,ws2,qs2) 
	H = zeros(2*Ntot,2*Ntot)*(1+0*1im)
	H[1:2*N1,1:2*N1] += H1

	H[2*N1-1,2*N1+1] += gS1[1,1]
	H[2*N1+1,2*N1-1] += gS1[1,1]

	H[2*N1-1,2*N1+2] += gS1[1,2]
	H[2*N1+2,2*N1-1] += gS1[1,2]

	H[2*N1,2*N1+2] += gS1[2,2]	
	H[2*N1+2,2*N1] += gS1[2,2]

	H[2*N1,2*N1+1] += gS1[2,1]
	H[2*N1+1,2*N1] += gS1[2,1]

	H[2*N1+1:2*N1+2*NS,2*N1+1:2*N1+2*NS] += HS

	H[2*N1+2*NS-1,2*N1+2*NS+1] += gS2[1,1]
	H[2*N1+2*NS+1,2*N1+2*NS-1] += gS2[1,1]

	H[2*N1+2*NS,2*N1+2*NS+2] += gS2[2,2]
	H[2*N1+2*NS+2,2*N1+2*NS] += gS2[2,2]
	
	H[2*N1+2*NS,2*N1+2*NS+1] += gS2[2,1]
	H[2*N1+2*NS+1,2*N1+2*NS] += gS2[2,1]

	H[2*N1+2*NS-1,2*N1+2*NS+2] += gS2[1,2]
	H[2*N1+2*NS+2,2*N1+2*NS-1] += gS2[1,2]

	H[2*(N1+NS)+1:end,2*(N1+NS)+1:end] += H2
	return H
end

function HatanoNelson(N,tA,tB)
	H = zeros(N,N)
	for i in 1:N-1
		H[i,i+1] = -tA
		H[i+1,i] = -tB
	end
	return H
end
function HatanoNelson2D(Nx,Ny,tA,tB)
    l = Nx*Ny
    T = zeros(l,l)
    for i in 1:l
        if i==1
            T[i,i+1] = -tA
        end
        if i==l
            T[i,i-1] = -tB
        end
        if i != 1 && i != l && i%x !=0 && i%x != 1
            T[i,i-1] = -tB
            T[i,i+1] = -tA
        end
        if i%x == 1
            T[i,i+1] = -tA
        end
        if i%x == 0
            T[i,i-1] = -tB
        end
        if i > x
            T[i,i-x] = -tB
        end
        if i <= l-x
            T[i,i+x] = -tA
        end
    end
    return T
end

function SingleLevel(E,ws,gph,q)
    return [E[1] + sqrt(2*ws)*gph*q 0.0; 0.0 E[2] + sqrt(2*ws)*gph*q]
end

function NNN(N,tns)
	H = zeros((2*N,2*N))
	for j in 1:length(tns)	
		for i in 1:2:2*N-2*j
			H[i,i+2*j] = -tns[j]
			H[i+1,i+2*j+1] = -tns[j]
			H[i+2*j,i] = -tns[j]
			H[i+2*j+1,i+1] = -tns[j]
		end
	end
	return H
end

function SpinHopping(N,eps,td,to)
        Ho = Onsite_A(N,eps)
        Hh = Hopping_A_nn_spin(N,td,to)
        return Ho .+ Hh
end


end
