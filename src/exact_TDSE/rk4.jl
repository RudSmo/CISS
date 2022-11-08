module RK4

include("helper.jl")
include("hamiltonian.jl")
using LinearAlgebra
using .Helper
using .Hamiltonian
using QuantumInformation
using HDF5

export Propagate, PropagateImag, PropagateKondo, PropagateKondoNI, InitialState, ElectronicStateNI, PhPosition, PhMomentum

function SchEq(psi0::Array,H::Array)::Array
    """
    Good ol' Schrödinger equation.
    Input:
        psi0 ::Array, the wavefunction
        H ::Array, the Hamiltonian
    Output:
        RHS of Schrödinger equation.
    """
    return -1im * H * psi0
end


function RK4Coeffs()::Array
    """
    Coefficients for RK4.
    """
    return [(1/3),
            (-1/3,1),
            (1,-1,1),
            (1/8,3/8,3/8,1/8)]
end

function RK4Step(H::Array,psi::Array,dt::Number)::Array
    """
    Makes a RK4 step. 
    Input:
        H ::Array, the Hamiltonian
        RK4 ::Array, RK4 Coefficients
        psi ::Array, wave-function at time t
        dt ::Float64, time-step
    Output:
        psi(t+dt)
    """
    RK4 = RK4Coeffs()
    rk1 = SchEq(psi,H)
    rk2 = SchEq(psi+dt* rk1*RK4[1][1],H)
    rk3 = SchEq(psi+dt*(rk1*RK4[2][1]+rk2*RK4[2][2]),H)
    rk4 = SchEq(psi+dt*(rk1*RK4[3][1]+rk2*RK4[3][2]+rk3*RK4[3][3]),H)
    rk = [rk1,rk2,rk3,rk4]
    psi_c = psi + dt * sum([rk[i]*RK4[4][i] for i in 1:4])
    return psi_c
end


function Measure(operator::Array, state::Array)::Float64
    """
    Expectation value of operator under wave-function state.

    Input:
        operator    ::Array     Operator
        state       ::Array     Wave-function

    Output:
        <state|operator|state>      ::Float64
    """
    R = real(state' * operator * state)
    return R[1]
end



function ElectronicState(s_state::Array,L::Int64)::Array
    """
    Creates an initial spin state |σ1,σ2,...,σL>, provided bit-array. 1 means up, 0 means down.
    Input:
        s_state ::Array, bit array (e.g. [0,1,1] for |↓↑↑>)
        L ::Array, lattice sites        
    """
    up = [1,0]
    dw = [0,1]
    st = 1
    for i in 1:2*L
        if s_state[i] == 1
            st = kron(st,up)
        else
            st = kron(st,dw)
        end
    end
    return st
end

function ElectronicStateNI(s_state::Array)::Array
    """
    Creates an initial spin state |σ1,σ2,...,σL>, provided bit-array. 1 means up, 0 means down.
    Input:
        s_state ::Array, bit array (e.g. [0,1,1] for |↓↑↑>)
        L ::Array, lattice sites        
    """
    psi_k = s_state[1:2]
    psi_e = s_state[3:end]
    up = [1,0]
    dw = [0,1]
    st = 1
    for i in 1:2
        if psi_k[i] == 1
            st = kron(st,up)
        else
            st = kron(st,dw)
        end
    end
    st = kron(st,psi_e)
    return st
end

function PhononState(ph_state::Array,M::Int64)::Array
    """
    Creates initial phonon state |n1,n2,n3,...,nL>, provided occupation-number array. Condition ni<M. 
    Input:
        ph_state::Array, phonon state (e.g. [0,5] for |0,5>, meaning HO on site 1 is in ground state, and on site 2 the 
        5th level is occupied only.)
        M ::Int64, total number of levels
    """
    L = length(ph_state)
    st = 1
    for i in 1:L
        ph = HOBasis(ph_state[i],M)
        st = kron(st,ph)
    end
    return st
end

function InitialState(el_state,ph_state,M,L)
    """
    Creates initial product state |σ1,σ2,...,σL> ⊗ |n1,n2,...,nL>.
    Input: 
        s_state ::Array, |σ1,σ2,...,σL>, as bit array
        ph_state::Array, |n1,n2,...,nL>, as occ. number array. 
    """
    return nkron([ElectronicState(el_state,L),PhononState(ph_state,M)])
end

function ElOccLattice(state::Array,L::Int64,M::Int64)::Array
    """
    Computes local spin-resolved fermionic occupation number.
    """
    IdPh = Matrix(I,M^L,M^L)
    return [[Measure(kron(ElNopUp(i,L),IdPh),state),Measure(kron(ElNopDw(i,L),IdPh),state)] for i in 1:L]
end

function ElOccLatticeKondo(state::Array,L::Int64)::Array
    """
    Computes local spin-resolved fermionic occupation number.
    """
    IdK = Matrix(I,4,4)
    return [[Measure(kron(IdK,ElNopUp(i,L)),state),Measure(kron(IdK,ElNopDw(i,L)),state)] for i in 1:L]
end

function ElOccLatticeKondoNI(state::Array,L::Int64)::Array
    rho = state * state' .+ 1e-15
    inds = [2,2,2*L]
    trinds = [1,2]
    rho_e = ptrace(rho,inds,trinds)
    return [real.([rho_e[i,i],rho_e[i+1,i+1]]) for i in 1:2:2*L-1]
end

function PhOccLattice(state::Array,L::Int64,M::Int64)::Array
    """
    Computes local total bosonic occupation number.
    """
    IdEl = Matrix(I,2^(2*L),2^(2*L))
    return [Measure(kron(IdEl,PhNop(M,L,i)),state) for i in 1:L]
end

function PhOccLevels(state::Array,L::Int64,M::Int64)::Array
    """
    Computes local level occupation on each site. Requires Pkg "QuantumInformation", because partial trace is employed.
    Input:
        state ::Array, wave-function
        L ::Array, lattice sites
        M ::Array, bosonic levels
    Output:
        occ ::Array = [occ1,occ2,...,occL], where occi = [occi_1,occi_2,...,occi_M] denotes the 
        occupations of each bosonic level on site i. 
    The occupations on site i are the diagonal entries of the reduced density matrix of the bosonic subspace on site i. 
    Note the convention used throughout this script.
    """
    Rho = state*state'
    occ = []
    idims = collect(Iterators.flatten([[2 for j in 1:2*L],[M for j in 1:L]]))
    sp_sites = Array(1:2*L)
    for i in 2*L+1:3*L
        trout = collect(Iterators.flatten([sp_sites,[j for j in 2*L+1:3*L if j!=i]])) 
        rho = QuantumInformation.ptrace(Rho,idims,trout)
        occ_l = diag(rho)
        push!(occ,real.(occ_l))
    end
    return occ
end

function PhPosition(psi0::Array,L::Int64,M::Int64)::Array
    ql = []
    for i in 1:L
#        bopL = nkron([Idkron(2*L,2),Idkron(i-1,M),bop(M),Idkron(L-i,M)])
#    	qop = 1/sqrt(2) * (bopL + bopL')
	bopL = 1/sqrt(2)*nkron([Idkron(2*L,2),Idkron(i-1,M),bop(M)+bop(M)',Idkron(L-i,M)])
#	push!(ql,Measure(qop,psi0))
	push!(ql,Measure(bopL,psi0))
    end
    return real.(ql)
end
function PhMomentum(psi0::Array,L::Int64,M::Int64)
    pl = []
    for i in 1:L
#        bopL = nkron([Idkron(2*L,2),Idkron(i-1,M),bop(M),Idkron(L-i,M)])
	bopL = -1/(sqrt(2)*1im) * nkron([Idkron(2*L,2),Idkron(i-1,M),bop(M)-bop(M)',Idkron(L-i,M)])

#        Pop = -1/(sqrt(2)*1im) * (bopL - bopL')
#	push!(pl,Measure(Pop,psi0))
	push!(pl,Measure(bopL,psi0))

    end
    return real.(pl)
end
function Loschmidt(psi0,psit)
	return psi0' * psit
end
function VonNeumannEntropy(psi::Array,L::Int64,M::Int64)
        rho = psi * psi' .+ 1e-15
        inds =[4 for i in 1:L]
        for i in 1:L
                push!(inds,M)
        end
	EE = []
	for bond in 1:2*L
		trinds = []
		for k in 1:2*L
			if k != bond
				push!(trinds,k)
			end
		end
        	rhoAB = ptrace(rho,inds,Int.(trinds))
        	push!(EE, real(-tr(rhoAB * log.(rhoAB))))
	end
	return EE
end


function Propagate(psi0,L,M,H,dt,filname,nsteps,nsave,nflush,pathIn = "")
    """
    Main RK4 propagation. Solves Schrödinger equation, computes expectation values, saves them to a file.
    """
    if pathIn != ""
        fid = h5open(pathIn,"r")
        psi0 = read(fid["psi"])
        close(fid)
    end
    nd = [ElOccLattice(psi0,L,M)]
    bl = [PhOccLevels(psi0,L,M)]
    ql = [PhPosition(psi0,L,M)]
    pl = [PhMomentum(psi0,L,M)]
    le = [Loschmidt(psi0,psi0)*(1+0.0*1im)]
    ee = [VonNeumannEntropy(psi0,L,M)]
    ti = [0.0]
    psi_p = psi0
    
    fil = open("$(filname).dat","w")
    write(fil,"#(1)t ")
    k = 1
    for i in 1:L
        write(fil,"($(k+1))ElN$(i)_up ($(k+2))ElN$(i)_dw ")
        k += 2
    end
    k = 1
    for i in 1:L
        for j in 0:M-1
            write(fil,"($(L+4+k))PhN$(i)_$(j) ")
            k += 1
        end
    end
    for i in 1:L
	write(fil,"($(L+4+k))PhNQ$(i) ($(L+4+k+1))PhNP$(i) ")
	k += 2
    end
    write(fil,"\n")
    close(fil)
    modeN = "w"
    normPsi = [norm(psi_p)]
    for i in 1:nsteps
        psi_n = RK4Step(H,psi_p,dt)
        psi_p = psi_n/norm(psi_n)
        if i % nsave == 0
            push!(nd,ElOccLattice(psi_p,L,M))
            #push!(bo,OccLattice(psi_p,L,M))
            push!(bl,PhOccLevels(psi_p,L,M))
	    push!(ql,PhPosition(psi_p,L,M))
	    push!(pl,PhMomentum(psi_p,L,M))
	    push!(le,Loschmidt(psi0,psi_p))
            push!(ee,VonNeumannEntropy(psi_p,L,M))
	    push!(normPsi,norm(psi_p))
            push!(ti,i*dt)
        end
        if i % nflush == 0
	    filN = open("$(filname)_norm_le.dat",modeN)
	    filEE = open("$(filname)_vonNeumann.dat",modeN)
	    modeN = "a"
	    fil = open("$(filname).dat","a")
            for j in 1:length(ti)
		write(filN,"$(ti[j]) $(normPsi[j]) $(real(le[j])) $(imag(le[j]))")
		write(filEE,"$(ti[j]) ")
		for k in 1:length(ee[j])
			write(filEE,"$(ee[j][k]) ")
		end
		write(filEE,"\n")

                write(fil,"$(ti[j]) ")
                for k in 1:length(nd[j])
                    write(fil,"$(nd[j][k][1]) $(nd[j][k][2]) ")
                end
               # for k in 1:length(bo[j])
               #     write(fil,"$(bo[j][k]) ")
               # end 
                for ki in 1:length(bl[j])
                    for li in 1:length(bl[j][ki])
                        write(fil,"$(bl[j][ki][li]) ")
                    end
                end
		for ki in 1:length(ql[j])
		    write(fil,"$(ql[j][ki]) $(pl[j][ki]) ")
		end
                write(fil,"\n")
		write(filN,"\n")
            end
            close(fil)
            close(filN)
	    close(filEE)
            nd = []
            bl = []
	    ql = []
	    pl = []
	    le = []
	    ee = []
	    normPsi = []
            ti = []
        end
        if i == nsteps
            fil = h5open("$(filname)_final_Quantities","w")
            fil["psi"] = psi_p
        end
    end
end



function KondoSpin(psi,L)
    rho = psi * psi' .+ 1e-15
    inds = Int.(2 * ones(2*(L+1)))
    trinds = Int.(Array(3:2*(L+1)))
    rho_spin = ptrace(rho,inds,Int.(trinds)) 
    c = fop()
    cd = c'
    id = Matrix(I,2,2)
    sx = 0.5 * [0 1; 1 0]
    sy = 0.5 * [0 -1im; 1im 0]
    sz = 0.5 * [1 0; 0 -1.]
    F = 2*sz

    Sx = tr(rho_spin * (sx[1,1] * kron(cd*c,id) + sx[1,2]*kron(cd*F,c) + sx[2,1]*kron(c*F,cd) + sx[2,2]*kron(id,cd*c)))
    Sy = tr(rho_spin * (sy[1,1] * kron(cd*c,id) + sy[1,2]*kron(cd*F,c) + sy[2,1]*kron(c*F,cd) + sy[2,2]*kron(id,cd*c)))
    Sz = tr(rho_spin * (sz[1,1] * kron(cd*c,id) + sz[1,2]*kron(cd*F,c) + sz[2,1]*kron(c*F,cd) + sz[2,2]*kron(id,cd*c)))  
   # Sx = tr(rho_spin * sx)
   # Sy = tr(rho_spin * sy)
   # Sz = tr(rho_spin * sz)
    return real.([Sx,Sy,Sz])
end

function KondoSpinNI(psi,L)
    rho = psi * psi' .+ 1e-15
    inds = [2,2,2*L]
    trinds = [3]
    rho_spin = ptrace(rho,inds,Int.(trinds)) 
    sx = 0.5 * [0 1; 1 0]
    sy = 0.5 * [0 -1im; 1im 0]
    sz = 0.5 * [1 0; 0 -1.]
    c = [0 1; 0 0.]
    cd = c' 
    id = Matrix(I,2,2)
    F = 2*sz
    Sx = tr(rho_spin * (sx[1,1] * kron(cd*c,id) + sx[1,2]*kron(cd*F,c) + sx[2,1]*kron(c*F,cd) + sx[2,2]*kron(id,cd*c)))
    Sy = tr(rho_spin * (sy[1,1] * kron(cd*c,id) + sy[1,2]*kron(cd*F,c) + sy[2,1]*kron(c*F,cd) + sy[2,2]*kron(id,cd*c)))
    Sz = tr(rho_spin * (sz[1,1] * kron(cd*c,id) + sz[1,2]*kron(cd*F,c) + sz[2,1]*kron(c*F,cd) + sz[2,2]*kron(id,cd*c)))   
    #Sx = tr(rho_spin * sx)
    #Sy = tr(rho_spin * sy)
    #Sz = tr(rho_spin * sz)
    return real.([Sx,Sy,Sz])
end

function PropagateKondo(psi0,L,H,dt,filname,nsteps,nsave,nflush,pathIn = "")
    """
    Main RK4 propagation. Solves Schrödinger equation, computes expectation values, saves them to a file.
    """
    if pathIn != ""
        fid = h5open(pathIn,"r")
        psi0 = read(fid["psi"])
        close(fid)
    end
    nd = [ElOccLatticeKondo(psi0,L)]
    sk = [KondoSpin(psi0,L)]
    ti = [0.0]
    psi_p = psi0
    
    fil = open("$(filname).dat","w")
    write(fil,"#(1)t ")
    k = 1
    for i in 1:L
        write(fil,"($(k+1))ElN$(i)_up ($(k+2))ElN$(i)_dw ")
        k += 2
    end
    write(fil,"($(k+1))Sx ($(k+2))Sy ($(k+3))Sz \n")
    close(fil)
    modeN = "w"
    for i in 1:nsteps
        psi_n = RK4Step(H,psi_p,dt)
        psi_p = psi_n # /norm(psi_n)
        if i % nsave == 0
            push!(nd,ElOccLatticeKondo(psi_p,L))
            push!(sk,KondoSpin(psi_p,L))
            push!(ti,i*dt)
        end
        if i % nflush == 0
	        fil = open("$(filname).dat","a")
            for j in 1:length(ti)
		        write(fil,"$(ti[j]) ")
                for k in 1:length(nd[j])
                    write(fil,"$(nd[j][k][1]) $(nd[j][k][2]) ")
                end
                for k in 1:3
                    write(fil,"$(sk[j][k]) ")
                end

                write(fil,"\n")
            end
            close(fil)
            nd = []
            sk = []
	        ti = []
        end
        if i == nsteps
            fil = h5open("$(filname)_final_Quantities","w")
            fil["psi"] = psi_p
        end
    end
end
function PropagateKondoNI(psi0,L,H,dt,filname,nsteps,nsave,nflush,pathIn = "")
    """
    Main RK4 propagation. Solves Schrödinger equation, computes expectation values, saves them to a file.
    """
    if pathIn != ""
        fid = h5open(pathIn,"r")
        psi0 = read(fid["psi"])
        close(fid)
    end
    nd = [ElOccLatticeKondoNI(psi0,L)]*1.0
    sk = [KondoSpinNI(psi0,L)]*1.0
    ti = [0.0]
    psi_p = psi0
    
    fil = open("$(filname).dat","w")
    write(fil,"#(1)t ")
    k = 1
    for i in 1:L
        write(fil,"($(k+1))ElN$(i)_up ($(k+2))ElN$(i)_dw ")
        k += 2
    end
    write(fil,"($(k+1))Sx ($(k+2))Sy ($(k+3))Sz \n")
    close(fil)
    modeN = "w"
    for i in 1:nsteps
        psi_n = RK4Step(H,psi_p,dt)
        psi_p = psi_n # /norm(psi_n)
        if i % nsave == 0
            push!(nd,ElOccLatticeKondoNI(psi_p,L))
            push!(sk,KondoSpinNI(psi_p,L))
            push!(ti,i*dt)
        end
        if i % nflush == 0
	        fil = open("$(filname).dat","a")
            for j in 1:length(ti)
		        write(fil,"$(ti[j]) ")
                for k in 1:length(nd[j])
                    write(fil,"$(nd[j][k][1]) $(nd[j][k][2]) ")
                end
                for k in 1:3
                    write(fil,"$(sk[j][k]) ")
                end

                write(fil,"\n")
            end
            close(fil)
            nd = []
            sk = []
	        ti = []
        end
        if i == nsteps
            fil = h5open("$(filname)_final_Quantities","w")
            fil["psi"] = psi_p
        end
    end
end


function PropagateImag(L,M,H,dt,filname,imagfilname,nsteps,imagsteps,nsave,nflush)
    """
    Main RK4 propagation with initial imaginary-time propagation. 
    Solves Schrödinger equation, computes expectation values, saves them to a file.
    """
    psi0 = rand(2^(2*L)*M^L)*(1+1im*0)
    psi0 = psi0/norm(psi0)
    nd = [ElOccLattice(psi0,L,M)]
    bl = [PhOccLevels(psi0,L,M)]
    ti = [0.0]
    psi_p = psi0
    
    fil = open(filname,"w")
    ifil = open(imagfilname,"w")
     
    write(fil,"#(1)t ")
    write(ifil,"#(1)t ")
    k = 1
    for i in 1:L
        write(fil,"($(k+1))ElN$(i)_up ($(k+2))ElN$(i)_dw ")
        write(ifil,"($(k+1))ElN$(i)_up ($(k+2))ElN$(i)_dw ")
        k += 2
    end
    k = 1
    for i in 1:L
        for j in 0:M-1
            write(fil,"($(L+3+k))PhN$(i)_$(j) ")
            write(ifil,"($(L+3+k))PhN$(i)_$(j) ")
            k += 1
        end
    end
    write(fil,"\n")
    write(ifil,"\n")
     
    close(fil)
    close(ifil)
     
    for i in 1:imagsteps
        psi_n = RK4Step(H,psi_p,-1im*dt)
        psi_p = psi_n/norm(psi_n)
        if i % nsave == 0
            push!(nd,ElOccLattice(psi_p,L,M))
            #push!(bo,OccLattice(psi_p,L,M))
            push!(bl,PhOccLevels(psi_p,L,M))
            push!(ti,i*dt)
        end
        if i % nflush == 0
            fil = open(imagfilname,"a")
            for j in 1:length(ti)
                write(fil,"$(ti[j]) ")
                for k in 1:length(nd[j])
                    write(fil,"$(nd[j][k][1]) $(nd[j][k][2]) ")
                end
               # for k in 1:length(bo[j])
               #     write(fil,"$(bo[j][k]) ")
               # end 
                for ki in 1:length(bl[j])
                    for li in 1:length(bl[j][ki])
                        write(fil,"$(bl[j][ki][li]) ")
                    end
                end
                write(fil,"\n")
            end
            close(fil)
            
            nd = []
            bl = []
            ti = []
        end    
    end
    psi0 = psi_p
    ti = [0.0]
    nd = [ElOccLattice(psi0,L,M)]
    bl = [PhOccLevels(psi0,L,M)]
    for i in 1:nsteps
        psi_n = RK4Step(H,psi_p,dt)
        psi_p = psi_n#/norm(psi_n)
        if i % nsave == 0
            push!(nd,ElOccLattice(psi_p,L,M))
            push!(bl,PhOccLevels(psi_p,L,M))
            push!(ti,i*dt)
        end
        if i % nflush == 0
            fil = open(filname,"a")
            for j in 1:length(ti)
                write(fil,"$(ti[j]) ")
                for k in 1:length(nd[j])
                    write(fil,"$(nd[j][k][1]) $(nd[j][k][2]) ")
                end
                for ki in 1:length(bl[j])
                    for li in 1:length(bl[j][ki])
                        write(fil,"$(bl[j][ki][li]) ")
                    end
                end
                write(fil,"\n")
            end
            close(fil)
            
            nd = []
            bl = []
            ti = []
        end
    end
end



end
