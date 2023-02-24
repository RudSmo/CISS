module RK4

using LinearAlgebra

export Propagate,bop,HOBasis


function Peierls(w0,nvib,dif,t1)
	a = bopD(nvib,t1,w0,dif)
        ad = a'
        H  = w0 * ad * a - t1 * dif * (ad + a)
        return H
end
function DisOp(alpha,nvib)
	a = bop(nvib)
	ad = a'
	return exp(alpha*a' - alpha' * a)
end
function bop(M)::Array
    """
    Bosonic annihilation operator with M levels.
    Input:
        M ::Int64, number of levels
    """
    b = zeros(M,M)
    for i in 1:M-1
                b[i,i+1] = sqrt(i)
    end
    return b
end
function bopD(M,t1,w0,dif)::Array
    """
    Bosonic annihilation operator with M levels.
    Input:
        M ::Int64, number of levels
    """
    alpha = sqrt(2*w0)/2 * t1 * sqrt(2/w0^3) * dif
    D = DisOp(alpha,M)	
    b = D' * bop(M) * D
    return b
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
    occ = [Rho[i,i] for i in 1:M]
    return real.(occ)
end

function PhXP(psi0::Array,M::Int64,w0::Float64,t1::Float64,dif::Float64)::Tuple
    a = bopD(M,t1,w0,dif)  
    X = 1/sqrt(2*w0) * (a' + a)
    P = -1im * sqrt(w0/2) * (a'-a)
    return Measure(X,psi0),Measure(P,psi0)
end
function Propagate(psi0,w0,t1,nvib,dif,dt,filname,nsteps,nsave,nflush)
    """
    Main RK4 propagation. Solves Schrödinger equation, computes expectation values, saves them to a file.
    """
    X,P = PhXP(psi0,nvib,w0,t1,dif[1])
    Xl = [X]
    Pl = [P]
    PhO = [PhOccLevels(psi0,nvib)]
    ti = [0.0]
    psi_p = psi0
    fil = open("$(filname)","w")
    write(fil,"(1)t (2)X (3)P ")
    for i in 1:nvib
	write(fil,"($(3+i))Occ$(i) ")
    end
    write(fil,"\n")

    modeN = "w"
    normPsi = [norm(psi_p)]
    for i in 1:nsteps
    	H  = Peierls(w0,nvib,dif[i],t1) 
        psi_n = RK4Step(H,psi_p,dt)
        psi_p = psi_n/norm(psi_n)
        if i % nsave == 0
	    X,P = PhXP(psi_p,nvib,w0,t1,dif[i])
	    push!(Xl,X)
	    push!(Pl,P)
	    push!(PhO,PhOccLevels(psi_p,nvib))
            push!(ti,i*dt)
        end
        if i % nflush == 0
	    fil = open("$(filname)","a")
            for j in 1:length(ti)
		write(fil,"$(ti[j]) $(Xl[j]) $(Pl[j]) ")
		for k in 1:length(PhO[j])
			write(fil,"$(PhO[j][k]) ")
		end
		write(fil,"\n")
            end
            close(fil)
	    Xl = []
	    Pl = []
	    PhO = []
            ti = []
        end
    end
    display(abs.(psi_p))
end

end
