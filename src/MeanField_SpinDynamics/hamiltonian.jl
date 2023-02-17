module Hamiltonian
using LinearAlgebra

export VDot,Hamiltonian_Kondo,Hamiltonian_KondoPeierls


function VDot(a)
    sx = [0 1; 1 0.]
    sy = [0 -1im; 1im 0.]
    sz = [1 0; 0 -1.] 
    return a[1]*sx + a[2]*sy + a[3]*sz
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


#-------------------------- Electron-Spin Hamiltonian-------------------------------------
function Kondo(N,Jsd,M,pos)
	H = zeros(2*N,2*N) * (1+0.0*1im)
	for i in 1:length(pos)
		a = 0.5 * VDot(M[i])	
		j = 2*pos[i]-1
		H[j,j] = a[1,1]
		H[j,j+1] = a[1,2]
		H[j+1,j] = a[2,1]
		H[j+1,j+1] = a[2,2]
	end
	return Jsd * H
end

function Hamiltonian_Kondo(HamPar)
    Ho = Onsite_A(HamPar.N,HamPar.OnsiteEnergies) 
    Hh = Hopping_A_nn(HamPar.N,HamPar.Hopping)
    Hsp = Kondo(HamPar.N,HamPar.Jsd,HamPar.mag,HamPar.pos)
    return Ho .+ Hh .+ Hsp
end

function KondoPeierls(N,Jsd,M,pos)
	H = zeros(2*N,2*N)*(1+0.0*1im)
	for i in 1:length(pos)
		if 2*pos[i]-1 <= 2*N-3
			a = 0.5 * VDot(M[i])
        	        j = 2*pos[i]-1
                	H[j,j+2] = a[1,1]
	                H[j,j+3] = a[1,2]	
        	        H[j+1,j+2] = a[2,1]
                	H[j+1,j+3] = a[2,2]
		else
			print("Position not valid. This is the Kondo-Peierls model. Reconsider the position of the spin such 
			       that the next-nearest coherence is coupled.")
		end
        end
        return Jsd * (H + H')
end
function Hamiltonian_KondoPeierls(HamPar)
    Ho = Onsite_A(HamPar.N,HamPar.OnsiteEnergies)
    Hh = Hopping_A_nn(HamPar.N,HamPar.Hopping)
    Hsp = KondoPeierls(HamPar.N,HamPar.Jsd,HamPar.mag,HamPar.pos)
    return Ho .+ Hh .+ Hsp
end

end
