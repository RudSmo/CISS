module RK4

using LinearAlgebra

export Propagate,bop,HOBasis


function Peierls(dx,Xrng,w0,dif,t1)
	N = length(Xrng)
	H = zeros(N,N)
	for i in 1:N
		H[i,i] = 1/2 * w0^2 * Xrng[i]^2 - 1/(2*dx^2) - sqrt(2*w0)*t1*Xrng[i]*dif
		if i < N
			H[i,i+1] = -1/(dx^2)
			H[i+1,i] = -1/(dx^2)
		end
	end
        return H
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

function Propagate(psi0,w0,t1,dif,dt,dx,Xrng,filname,nsteps,nsave,nflush)
    """
    Main RK4 propagation. Solves Schrödinger equation, computes expectation values, saves them to a file.
    """
    Xl = [abs.(psi0)]
    ti = [0.0]
    psi_p = psi0
    fil = open("$(filname)","w")
    write(fil,"(1)t ")
    N = length(Xrng)
    for i in 1:N
	write(fil,"($(i+1))X$(i) ")
    end
    write(fil,"\n")
    modeN = "w"
    for i in 1:nsteps
    	H  = Peierls(dx,Xrng,w0,dif[i],t1) 
        psi_n = RK4Step(H,psi_p,dt)
        psi_p = psi_n /norm(psi_n)
        if i % nsave == 0
	    X = abs.(psi_p)
	    push!(Xl,X)
            push!(ti,i*dt)
        end
        if i % nflush == 0
	    fil = open("$(filname)","a")
            for j in 1:length(ti)
		write(fil,"$(ti[j]) ")
		for k in 1:length(Xl[j])
			write(fil,"$(Xl[j][k]) ")
		end
		write(fil,"\n")
            end
            close(fil)
	    Xl = []
            ti = []
        end
    end
end

end
