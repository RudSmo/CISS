module RK4

include("helper.jl")
include("hamiltonian.jl")
using LinearAlgebra
using .Helper
using .Hamiltonian
using QuantumInformation
export Propagate, PropagateImag, InitialState

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

function Propagate(psi0,L,M,H,dt,filname,nsteps,nsave,nflush)
    """
    Main RK4 propagation. Solves Schrödinger equation, computes expectation values, saves them to a file.
    """
    nd = [ElOccLattice(psi0,L,M)]
    bl = [PhOccLevels(psi0,L,M)]
    ti = [0.0]
    psi_p = psi0
    
    fil = open(filname,"w")
    write(fil,"#(1)t ")
    k = 1
    for i in 1:L
        write(fil,"($(k+1))ElN$(i)_up ($(k+2))ElN$(i)_dw ")
        k += 2
    end
    k = 1
    for i in 1:L
        for j in 0:M-1
            write(fil,"($(L+3+k))PhN$(i)_$(j) ")
            k += 1
        end
    end
    write(fil,"\n")
    close(fil)

    for i in 1:nsteps
        psi_n = RK4Step(H,psi_p,dt)
        psi_p = psi_n#/norm(psi_n)
        if i % nsave == 0
            push!(nd,ElOccLattice(psi_p,L,M))
            #push!(bo,OccLattice(psi_p,L,M))
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