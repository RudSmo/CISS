module Hamiltonian
using LinearAlgebra

export VDot,absv,UnitVector,FranssonHamiltonian,SingleLevel


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
function Holstein_NNOld(N,tph,ws,qs)
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
    H[2*N,2*N-2] = -tph * sqrt(2*ws[N])*qs[N]
    H[2*N-2,2*N] = -tph * sqrt(2*ws[N])*qs[N]
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



function FranssonHamiltonian(N,E,t,tso,a,c,num_laps,tph,tsoph,ws,qs,ver)
	Ho = Onsite_A(N,E)
	Hh = Hopping_A_nn(N,t)
	HSo = SOC_A_nn(N,tso,a,c,num_laps)
	if ver == "0"
    		Heph = Holstein_NNOld(N,tph,ws,qs)
	else
		Heph = Holstein_NN(N,tph,ws,qs)
	end
    	Hsoph = Holstein_SOC(N,tsoph,a,c,num_laps,ws,qs)
	return Ho+Hh+HSo+Heph+Hsoph
end


function SingleLevel(E,ws,gph,q)
    return [E[1] + sqrt(2*ws)*gph*q 0.0; 0.0 E[2] + sqrt(2*ws)*gph*q]
end

end
