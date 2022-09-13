module JaynesCumming

using LinearAlgebra
using QuantumInformation

function TLSBasis()
    e = [1, 0]
    g = [0, 1]
    return e,g
end

function PhotonOperators(N)
    a = zeros(N,N)
    for i in 1:N-1
        a[i,i+1] = sqrt(i)
    end 
    return a,a'
end

function HOBasis(k::Int64,M::Int64)::Array
    """
    Generates harmonic oscillator basis state for M levels, occupied at k.
    Input:
        k ::Int64, occupied level
        M ::Int64, total number of levels
    """
    st = zeros(M)
    st[k+1] = 1.0
    return st
end

function JCHamiltonian(wa,Ee,Eg,gep,N)
    IdP = Matrix(I,N,N)
    IdE = Matrix(I,2,2)

    e,g = TLSBasis()
    a,ad = PhotonOperators(N)
    H = wa * kron(ad*a,IdE) + Ee*kron(IdP,e*e') 
    H += Eg*kron(IdP,g*g') 
    H += gep*(kron(a,g*e')+kron(ad,e*g'))

    return H
end

function SchEq(psi,H)
    return -1im*H*psi
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

function ElectronicState(s_state::Array)::Array
    """
    Creates an initial spin state |σ1,σ2,...,σL>, provided bit-array. 1 means up, 0 means down.
    Input:
        s_state ::Array, bit array (e.g. [0,1,1] for |↓↑↑>)
        L ::Array, lattice sites        
    """
   # e,g = TLSBasis()
   # if s_state == "1"
   #     st = e#kron(e,e)
   # elseif s_state == "0"
   #     st = g#kron(g,g)
   # elseif s_state == "01"
   #     st = #kron(g,e) + kron(e,g)
   #     st = st/norm(st)
   # elseif s_state == "0+1"
   #     st = kron(g,g) + kron(e,e)
   #     st = st/norm(st)
   # end
    return s_state
end

function PhononState(ph_state::Number,M::Int64)::Array
    """
    Creates initial phonon state |n1,n2,n3,...,nL>, provided occupation-number array. Condition ni<M. 
    Input:
        ph_state::Array, phonon state (e.g. [0,5] for |0,5>, meaning HO on site 1 is in ground state, and on site 2 the 
        5th level is occupied only.)
        M ::Int64, total number of levels
    """
    ph = HOBasis(ph_state,M)
    return ph
end

function nkron(A::Array)::Any
    """
    Tensor product of all arrays, stored in A.
    Input:
        A ::Array = [A1,A2,...,An]
    Output:
        C ::Array = A1⊗A2⊗...⊗An
    """
    if length(A) == 1
        return A[1]
    else
        return kron(A[1],nkron(A[2:end]))
    end
end
function InitialState(el_state,ph_state,M)
    """
    Creates initial product state |σ1,σ2,...,σL> ⊗ |n1,n2,...,nL>.
    Input: 
        s_state ::Array, |σ1,σ2,...,σL>, as bit array
        ph_state::Array, |n1,n2,...,nL>, as occ. number array. 
    """
    return nkron([PhononState(ph_state,M),ElectronicState(el_state)])
end

function ElOccLattice(state::Array,M::Int64)::Array
    """
    Computes local spin-resolved fermionic occupation number.
    """
    IdPh = Matrix(I,M,M)
    e,g = TLSBasis()
    ee = e*e'
    gg = g*g'
    return [Measure(kron(IdPh,ee),state),Measure(kron(IdPh,gg),state)]
end


function PhOccLevels(state::Array,M::Int64)::Array
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
    idims = collect(Iterators.flatten([[M],[2]]))
    trout = collect(Iterators.flatten([2])) 
    rho = QuantumInformation.ptrace(Rho,idims,trout)
    occ = real.(diag(rho))
    return occ
end

function Propagate(psi_p,M,H,dt,filname,nsteps,nsave,nflush)
    """
    Main RK4 propagation. Solves Schrödinger equation, computes expectation values, saves them to a file.
    """
    
    nd = [ElOccLattice(psi_p,M)]
    bl = [PhOccLevels(psi_p,M)]
    ti = [0.0]
    
    fil = open(filname,"w")
    write(fil,"#(1)t ")
    write(fil,"(2)E (3)G ")
    k = 1
    for j in 0:M-1
        write(fil,"($(3+j+1))PhM$(j) ")
    end
    write(fil,"\n")
    close(fil)

    for i in 1:nsteps
        psi_n = RK4Step(H,psi_p,dt)
        psi_p = psi_n#/norm(psi_n)
        if i % nsave == 0
            push!(nd,ElOccLattice(psi_p,M))
            push!(bl,PhOccLevels(psi_p,M))
            push!(ti,i*dt)
        end
        if i % nflush == 0
            fil = open(filname,"a")
            for j in 1:length(ti)
                write(fil,"$(ti[j]) ")
                for k in 1:length(nd[j])
                    write(fil,"$(nd[j][k]) ")
                end
                for ki in 1:length(bl[j])
                    write(fil,"$(bl[j][ki][1]) ")
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

wa = 1.0
Ee = 1.0
Eg = -1.0
gep = 1.0
N = 8
H = JCHamiltonian(wa,Ee,Eg,gep,N)

el_state = [0,1.]
ph_state = 7
psi0 = InitialState(el_state,ph_state,N)
print(PhOccLevels(psi0,N))
dt = 0.001
filname = "wa$(wa)_Ee$(Ee)_Eg$(Eg)_gep$(gep)_N$(N).dat"
tend = 10
nsteps = Int(tend/dt)
nsave = 10
nflush = 100

Propagate(psi0,N,H,dt,filname,nsteps,nsave,nflush)

end