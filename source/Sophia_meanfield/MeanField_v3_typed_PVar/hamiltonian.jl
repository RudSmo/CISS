using LinearAlgebra

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
# -------------------------------Sophias auxilary functions-----------------------------------
function sigma_ort(
	phi::Float64, 
	theta::Float64
	)
    sig_x = Matrix([0 1; 1 0])
    sig_y = Matrix([0 -1im; 1im 0 ])
    sig_z = Matrix([1 0; 0 -1])
    return sig_x * sin(phi)*sin(theta) - sig_y* cos(phi)*sin(theta) + sig_z* cos(theta)
end
function sites_to_angle(sites_per_lap::Integer)# number of sites per helical turn
    Radius = 0.7* 10^(-9)
    delta_phi = pi* 2/sites_per_lap
    l_a = Radius/sites_per_lap
    return [delta_phi l_a]
end  
#-------------------------- Electronic Hamiltonian-------------------------------------


#--------------- if onsight e was ordered as x = Array[j = 1, j = 2]-----------
function H_Guo_dna_sort2(
	onsight_e::Array{Float64},
        thopping::Array{Float64}, 
        lambda::Float64,
        N::Integer
        )
   
    Mtrx = zeros(Complex{Float64}, 2*N,2*N)
    for i in 1:4:2*N # i = position of site n bloc in grand matrix
        i_steps = Integer((i+3)/4) # number of a site. Is one number per i
       # print(i_steps)
        Mtrx[i+2, i] = lambda # lambda is interstrand (inside one strand) hopping
        Mtrx[i, i+2] = lambda
        Mtrx[i+3, i+ 1] = lambda # j = 2
        Mtrx[i+ 1, i+3] = lambda # j = 2

        for spin in 1:1:2 # spin = 1 = up and spin = 2 = down 

        #diagonaleinträge
        Mtrx[(i+ spin -1), (i+ spin - 1)] = onsight_e[1] #  j = 1
        Mtrx[(i+ spin +1), (i+ spin +1)] = onsight_e[2] # j = 2
        # ------------------------------------------------use N/2 since N is the length of aray spin up or down 
        # for [j = 1, j=2] after half length of that array we switch to the other strand !
            
            if i<2*N-4# t is intrastrand hopping ( between strands)
            #rechte hopping matrix
            Mtrx[(i+3+spin), (spin+i-1)] = thopping[1] # j = 1
            Mtrx[(i+5+spin), (spin+i+1)] = thopping[2] # j = 2

            #untere hoping matrix
            Mtrx[ (spin+i-1), (i+3+spin)] = thopping[1] # j = 1
            Mtrx[ (spin+i+1), (i+5+spin)] = thopping[2] # j = 2

            end   
        end
       
    end
    return Matrix(Mtrx)
end

function H_Guo_so_sort2(
	t_so::Float64,
        sites_per_lap::Integer, 
        theta::Float64, 
        N::Int64 # total number of sites in both chains
        )
    Mtrx = zeros(Complex{Float64}, 2*N,2*N)# since we have two states: up and down for each site  
   
    alpha = 1
    l_a = sites_to_angle(sites_per_lap)[2]
    
    delta_phi = sites_to_angle(sites_per_lap)[1]
    
    # i = position of site n in full matrix, i_steps =  number of site n
    for i in 1:4:2*N 
        i_steps = Integer((i+3)/4)
        if i<2*N-4
            for z in 0:2:2 
                # für z  = 0 git j = 1 und z = 2 ist j = 2 
                if z ==0
                    #rechte hopping matrix
                    #* ist dies die h.c. seite oder die unteren matrizen?
                    sig_n_und_npl1 = 1im* (sigma_ort( i_steps*delta_phi, theta) .+ sigma_ort(( i_steps -1 ) *delta_phi, theta) )
                   #adresses a 2x2 block matrix
                    Mtrx[(i+z):(i+z+1), (i+4+z):(i+5+z)] = t_so*sig_n_und_npl1'
                    
                    #untere hopping matrix -> ohne compl conj?! 
                   Mtrx[(i+z+4):(i+z+5), (i+z):(i+1+z)] = t_so*(sig_n_und_npl1 )
                   
                end
                if z == 2
                    #plus pi in both sigma functions. As pi is half a circle adding this to phi results being on the other helical chain
                    sig_n_und_npl1 = 1im* (sigma_ort( i_steps*delta_phi + pi, theta) .+ sigma_ort(( i_steps -1 ) *delta_phi + pi, theta))
                    Mtrx[(i+z):(i+z+1), (i+4+z):(i+5+z)] = t_so*sig_n_und_npl1'
                    Mtrx[(i+z+4):(i+z+5), (i+z):(i+1+z)] = t_so*sig_n_und_npl1 
                end
            end
        end   
    end
    return Matrix(Mtrx)
end 


function H_Guo_elvib_sort2(
        g_coupling::Float64, 
        X::Array{Float64}, # array with enumerating over n first half of array for j = 1, second half for j = 2 starting with n=0 again
        N::Int64 # number of total sites in both chains
    )
    Mtrx = zeros(Complex{Float64},2*N, 2*N)
    
    for i in 1:4:2*N # i = number of site n bloc * 4 -> position in grand matrix
        i_steps = Integer((i+3)/4)
        for k in 1:1:4 # k = indexing  inside the bloc
            #diagonaleinträge
            if k == 1 || k ==2 #|| heißt osder, da spin up und spin down hier keine unterschiede aufweisen
            #erster und zweiter Blockelement der Untermatriyx j = 1 --> x[i_ssteps]
                Mtrx[(i+k-1), (i+k-1)] = X[i_steps] *g_coupling/ sqrt(2)
            else # für j = 2 -> ab hälfte der x liste ---> x[N + i_steps]
                Mtrx[(i+k-1), (i+k-1)] = X[Integer(N/2) + i_steps]*g_coupling/ sqrt(2)
            end
        end
    end
    return Matrix(Mtrx)
end
# is an imaginary matrix
function H_total_qm(onsight_e::Array{Float64},
        thopping::Array{Float64}, 
        lambda::Array{Float64},
        t_so::Float64,
	sites_per_turn::Integer, 
        theta::Float64,
        g_coupling::Float64, 
        X::Array{Float64},
        N::Integer
	)
	H_dna = H_Guo_dna_sort2(onsight_e, thopping, lambda, N )
        H_so = H_Guo_so_sort2(t_so, sites_per_turn, theta , N)
        H_elvib = H_Guo_elvib_sort2( g_coupling, X,  N )
        return H_dna .+ H_so .+ H_elvib
end

function H_total_hampar(HamPar::Any
	)
	H_dna = H_Guo_dna_sort2(HamPar.onsight_e, HamPar.thopping, HamPar.lambda, HamPar.N )
        H_so = H_Guo_so_sort2(HamPar.t_so, HamPar.sites_per_turn, HamPar.theta, HamPar.N)
        H_elvib = H_Guo_elvib_sort2(HamPar.g_coupling, HamPar.X,  HamPar.N )
        return H_dna .+ H_so .+ H_elvib
end

function Guo_Mode_v2(HamPar::Any,k::Int64)::Array 
        Xn = zeros(length(HamPar.X))
        Xn[k] = 1.0 # ableitung nach x = 1
        Heph = H_Guo_elvib_sort2(HamPar.g_coupling, Xn,  HamPar.N)
        # H_so not dependend on SO in this model?
        return Heph
end
#---------------------------------Rudis code --------------------------------------------

function SOC_A_nn(N,tso,a,c,num_laps)
    Hs = zeros(2*N,2*N)*(1+0.0*1im)
    for k in 1:N-2
	j = 2*k-1
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
    end
    return Hs+Hs'
end
function Hopping_A_nn(N,t)
    Hh = zeros(2*N,2*N)
    for i in 1:2*N-2
        Hh[i,i+2] = t
        Hh[i+2,i] = t
    end
    return Hh
end
function Hopping_A_nn_spin(N,t)
    tu,td = t
    Hh = zeros(2*N,2*N)
    for i in 1:2*N-2
        Hh[i,i+2] = tu
        Hh[i+2,i] = td
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
function Holstein(N,gPh,ws,xs)
    H = zeros(2*N,2*N)
    for k in 1:N
    	i = 2*k-1    
	H[i,i] = gPh[1] *  xs[k]/sqrt(2)
        H[i+1,i+1] = gPh[2] * xs[k]/sqrt(2)
    end
    return H
end


function Peierls(N,tph,ws,xs)
    H = zeros(2*N,2*N)
    for k in 1:N-1
	i = 2*k-1
        H[i,i+2] = tph[1] * xs[k]/sqrt(2) 
	H[i+1,i+3] = tph[2] * xs[k]/sqrt(2)
        H[i+2,i] = tph[1] * xs[k]/sqrt(2)
	H[i+3,i+1] = tph[2] * xs[k]/sqrt(2)
    end
    return H
end


function Peierls_SOC(N,tsoph,a,c,num_laps,ws,xs)
    Hs = zeros(2*N,2*N)*(1+0.0*1im)
    for k in 1:N-2
	j = 2*k-1
        rj = UnitVector(a,c,k,num_laps,N)
        rj1 = UnitVector(a,c,k+1,num_laps,N)
        rj2 = UnitVector(a,c,k+2,num_laps,N) 
        rjrj1 = rj-rj1
        rjrj2 = rj-rj2

        dj1 = rjrj1/absv(rjrj1)
        dj2 = rjrj2/absv(rjrj2)
        
        vj = cross(dj1,dj2) 
        Vj = 1im*tsoph*VDot(vj) * xs[k]/sqrt(2)

        Hs[j,j+4] = Vj[1,1]
        Hs[j+1,j+4] = Vj[2,1]
        Hs[j,j+5] = Vj[1,2]
        Hs[j+1,j+5] = Vj[2,2]  
    end
    return Hs+Hs'
end



function FranssonHamiltonian(HamPar::Any)::Array
	Ho = Onsite_A(HamPar.N,HamPar.OnsiteEnergies)
	Hh = Hopping_A_nn(HamPar.N,HamPar.Hopping)
	HSo = SOC_A_nn(HamPar.N,HamPar.SOC,HamPar.a_hel,HamPar.c_hel,HamPar.num_laps)
    	Heph = Peierls(HamPar.N,HamPar.HoppingPh,HamPar.ws,HamPar.X)
    	Hsoph = Peierls_SOC(HamPar.N,HamPar.SOCPh,HamPar.a_hel,HamPar.c_hel,HamPar.num_laps,HamPar.ws,HamPar.X)
	return Ho .+ Hh .+ HSo .+ Heph .+ Hsoph
end
function FranssonHamiltonian_Mode(HamPar,k)::Array
	Xn = zeros(length(HamPar.X))
	Xn[k] = HamPar.X[k]	
	Heph = Peierls(HamPar.N,HamPar.HoppingPh,HamPar.ws,Xn)
	Hsoph = Peierls_SOC(HamPar.N,HamPar.SOCPh,HamPar.a_hel,HamPar.c_hel,HamPar.num_laps,HamPar.ws,Xn)
	return Heph .+ Hsoph
end
# =#
function FranssonHamiltonian_Mode_v2(HamPar::Any,k::Int64)::Array
        Xn = zeros(length(HamPar.X))
        Xn[k] = 1.0
        Heph = Peierls(HamPar.N,HamPar.HoppingPh,HamPar.ws,Xn)
        Hsoph = Peierls_SOC(HamPar.N,HamPar.SOCPh,HamPar.a_hel,HamPar.c_hel,HamPar.num_laps,HamPar.ws,Xn)
        return Heph .+ Hsoph
end
##=
function FranssonHamiltonian_Mode_v3(HamPar::Any,k::Int64,X::Float64)::Array
        Xn = zeros(length(HamPar.X))
        Xn[k] = X
        Heph = Peierls(HamPar.N,HamPar.HoppingPh,HamPar.ws,Xn)
        Hsoph = Peierls_SOC(HamPar.N,HamPar.SOCPh,HamPar.a_hel,HamPar.c_hel,HamPar.num_laps,HamPar.ws,Xn)
        return Heph .+ Hsoph
end
	
function FranssonHamiltonian_H0(HamPar)
        Ho = Onsite_A(HamPar.N,HamPar.OnsiteEnergies)
        Hh = Hopping_A_nn(HamPar.N,HamPar.Hopping)
        HSo = SOC_A_nn(HamPar.N,HamPar.SOC,HamPar.a_hel,HamPar.c_hel,HamPar.num_laps)
        return Ho+Hh+HSo
end

function FranssonHamiltonian_HI(HamPar)
    	Heph = Peierls_Old(HamPar.N,HamPar.HoppingPh,HamPar.ws,HamPar.X)
    	Hsoph = Peierls_SOC(HamPar.N,HamPar.SOCPh,HamPar.a_hel,HamPar.c_hel,HamPar.num_laps,HamPar.ws,HamPar.X)
        return Heph+Hsoph
end

function FranssonHamiltonianNH(HamPar)
    HF = FranssonHamiltonian(HamPar)
    HF[1,1] -= 1im/2* HamPar.GL[1]
    HF[2,2] -= 1im/2* HamPar.GL[2]
    HF[end-1,end-1] -= 1im/2* HamPar.GR[1]
    HF[end,end] -= 1im/2 * HamPar.GR[2]
    return HF
end

function HolsteinHopping(HamPar)
    Ho = Onsite_A(HamPar.N,HamPar.OnsiteEnergies)
    Hh = Hopping_A_nn(HamPar.N,HamPar.Hopping)
    Hph = Holstein(HamPar.N,HamPar.gPh,HamPar.ws,HamPar.X)
    return Ho .+ Hh .+ Hph
end
# =#
function Holstein_Mode_v2(HamPar::Any,k::Int64)::Array
    Xn = zeros(length(HamPar.X))
    Xn[k] = 1.0
    Heph = Holstein(HamPar.N,HamPar.gPh,HamPar.ws,Xn)
    return Heph
end
##=

function SpinHopping(HamPar,t0)
	Ho = Onsite_A(HamPar.N,HamPar.OnsiteEnergies)
	Hh = Hopping_A_nn_spin(HamPar.N,t0)
	return Ho .+ Hh	
end

# =#

