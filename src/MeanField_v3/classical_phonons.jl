module ClassicalPhonons
include("hamiltonian.jl")
using LinearAlgebra
using .Hamiltonian
export SampleWigner_Multi,SymplecticEulerStep,SymplecticVerletStep,SymplecticRuth3Step,SymplecticRuth4Step,RK4_intertwined,RK4_sequential

function SampleWigner(beta,w)
	a = tanh(beta*w/2)
	u1,u2 = rand(2)
	z1 = 1/sqrt(a)*sqrt(-2*log(u1))*cos(2*pi*u2)
	z2 = 1/sqrt(a)*sqrt(-2*log(u1))*sin(2*pi*u2)
	return z1,z2
end

function SampleWigner_Multi(beta,w,m)
	W = [SampleWigner(beta,w[i]) for i in 1:m]
	return [W[i][1] for i in 1:m],[W[i][2] for i in 1:m]
end


function HigherOrderSymplecticStep(HamPar,rho,dt,coeffs)
	Pn = []
        Xn = []
        for i in 1:HamPar.N
            H_red = Hamiltonian.FranssonHamiltonian_Mode_v2(HamPar,i)
            derv  = tr(rho*H_red)
            Pn_c  = HamPar.P[i]
            Xn_c  = HamPar.X[i]
            for (ai,bi) in zip(coeffs[1],coeffs[2])
                    Pn_c += -bi*dt*(derv + Xn_c)
                    Xn_c += ai*dt*Pn_c
            end
            push!(Pn,real(Pn_c))
            push!(Xn,real(Xn_c))
        end
        return Xn,Pn
end

function SymplecticEulerStep(HamPar,rho,dt)
        coeffs = [[1.0],[1.0]]
        HigherOrderSymplecticStep(HamPar,rho,dt,coeffs)
end
function SymplecticVerletStep(HamPar,rho,dt)
	coeffs = [[0,1.],[0.5,0.5]]
	HigherOrderSymplecticStep(HamPar,rho,dt,coeffs)
end
function SymplecticRuth3Step(HamPar,rho,dt)
	coeffs = [[1,-2/3,2/3],[-1/24,3/4,7/24]]
	HigherOrderSymplecticStep(HamPar,rho,dt,coeffs)
end
function SymplecticRuth4Step(HamPar,rho,dt)
	a = 1/(2-2^(1/3))
        coeffs = a*[[1/2,(1-2^(1/3))/2,(1-2^(1/3))/2,1/2],
		    [0.,1.,-(2)^(1/3),1.]]
        HigherOrderSymplecticStep(HamPar,rho,dt,coeffs)
end

function P_RHS(derv,X)
	return -(derv+X)	
end
function X_RHS(P)
	return P
end
function RK4_intertwined(HamPar,rho,dt)	
	Pn = []
	Xn = []
	for i in 1:HamPar.N
		H_red = Hamiltonian.FranssonHamiltonian_Mode_v2(HamPar,i)
		derv  = tr(rho*H_red)
		p1 = dt*P_RHS(derv,HamPar.X[i])	
		x1 = dt*p1
		p2 = dt*P_RHS(derv,x1 + 1/3 * p1)
                x2 = dt*(p2 + 1/3 * x1)
		p3 = dt*P_RHS(derv,x2 - 1/3 * p1 + p2)
		x3 = dt*(p3 - 1/3 * x1 + x2)
		p4 = dt*P_RHS(derv,x3 + 1.0 * p1 - p2 + p3)	
		Pn_c = HamPar.P[i] + 1/8 * (p1 + 3*p2 + 3*p3 + p4)
                x4 = dt*(p4 + 1.0 * x1 - x2 + x3)
                Xn_c = HamPar.X[i] + 1/8 * (x1 + 3*x2 + 3*x3 + x4)
		push!(Pn,real(Pn_c))
		push!(Xn,real(Xn_c))
	end
	return Xn,Pn		
end
function RK4_sequential(HamPar,rho,dt)
        Pn = []
        Xn = []
        for i in 1:HamPar.N
                H_red = Hamiltonian.FranssonHamiltonian_Mode_v2(HamPar,i)
                derv  = tr(rho*H_red)
                p1 = dt*P_RHS(derv,HamPar.X[i])
                p2 = dt*P_RHS(derv,HamPar.X[i] + 1/3 * p1)
                p3 = dt*P_RHS(derv,HamPar.X[i] - 1/3 * p1 + p2)
                p4 = dt*P_RHS(derv,HamPar.X[i] + 1.0 * p1 - p2 + p3)
                Pn_c = HamPar.P[i] + 1/8 * (p1 + 3*p2 + 3*p3 + p4)

                x1 = dt*HamPar.P[i]
                x2 = dt*(HamPar.P[i] + 1/3 * x1)
                x3 = dt*(HamPar.P[i] - 1/3 * x1 + x2)
                x4 = dt*(HamPar.P[i] + 1.0 * x1 - x2 + x3)
                Xn_c = HamPar.X[i] + 1/8 * (x1 + 3*x2 + 3*x3 + x4)
                push!(Pn,real(Pn_c))
                push!(Xn,real(Xn_c))
        end
        return Xn,Pn
end




end
