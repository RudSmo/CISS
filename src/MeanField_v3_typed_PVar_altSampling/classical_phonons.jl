include("hamiltonian.jl")

using LinearAlgebra



function SampleWigner(beta::Float64,w::Float64)::Tuple{Float64,Float64}
	a = tanh((beta/w)*w/2)
	u1,u2 = rand(2)
	z1 = 1/sqrt(a)*sqrt(-2*log(u1))*cos(2*pi*u2)
	z2 = 1/sqrt(a)*sqrt(-2*log(u1))*sin(2*pi*u2)
	return z1,z2
end

function SampleGS()::Tuple{Float64,Float64}
        u1,u2 = rand(2)
        z1 = sqrt(-2*log(u1))*cos(2*pi*u2)
        z2 = sqrt(-2*log(u1))*sin(2*pi*u2)
        return z1,z2
end



function SampleWigner_Multi(beta::Float64,w::Array{Float64},m::Int64)::Tuple{Array{Float64},Array{Float64}}
	Xn = zeros(m)
	Pn = zeros(m)
	for i in 1:m
		W = SampleWigner(beta,w[i])
		Xn[i] = W[1]
		Pn[i] = W[2]
	end
	return Xn, Pn	
end
function SampleWigner_Multi_ntraj(beta::Float64,w::Array{Float64},m::Int64,num_avg::Int64)::Tuple{Array{Float64},Array{Float64}}
        Xn = zeros(m,num_avg)
        Pn = zeros(m,num_avg)
        for j in 1:num_avg
                for i in 1:m
                        W = SampleWigner(beta,w[i])
                        Xn[i,j] = W[1]
                        Pn[i,j] = W[2]
                end
        end
        return Xn, Pn
end


function SampleGS_Multi(m::Int64)::Tuple{Array{Float64},Array{Float64}}
        Xn = zeros(m)
        Pn = zeros(m)
        for i in 1:m
                W = SampleGS()
                Xn[i] = W[1]
                Pn[i] = W[2]
        end
        return Xn, Pn
end


function HigherOrderSymplecticStep(HamPar::Any,rho::Array{ComplexF64},dt::Float64,coeffs::Array)::Tuple{Array{Float64}, Array{Float64}}
	Pn = similar(HamPar.P)
        Xn = similar(HamPar.X)
        for i in 1:HamPar.N
            H_red = FranssonHamiltonian_Mode_v2(HamPar,i)
            derv  = tr(@views rho*H_red)
	    Xn_c = HamPar.X[i]
	    Pn_c = HamPar.P[i]
            for j in 1:length(coeffs[1])
                    Pn_c += real(-coeffs[2][j]*dt*(derv + Xn_c))
                    Xn_c += real(coeffs[1][j]*dt*Pn_c)
            end
	    Pn[i] = Pn_c
	    Xn[i] = Xn_c
        end
	return Xn,Pn
end

function SymplecticEulerStep(HamPar::Any,rho::Array,dt::Float64)
        coeffs = [[1.0],[1.0]]
        HigherOrderSymplecticStep(HamPar,rho,dt,coeffs)
end
function SymplecticVerletStep(HamPar::Any,rho::Array,dt::Float64)
	coeffs = [[0,1.],[0.5,0.5]]
	HigherOrderSymplecticStep(HamPar,rho,dt,coeffs)
end
function SymplecticRuth3Step(HamPar::Any,rho::Array,dt::Float64)
	coeffs = [[1,-2/3,2/3],[-1/24,3/4,7/24]]
	HigherOrderSymplecticStep(HamPar,rho,dt,coeffs)
end
function SymplecticRuth4Step(HamPar::Any,rho::Array{ComplexF64},dt::Float64)
	a = 1/(2-2^(1/3))
        coeffs = a*[[1/2,(1-2^(1/3))/2,(1-2^(1/3))/2,1/2],
		    [0.,1.,-(2)^(1/3),1.]]
        HigherOrderSymplecticStep(HamPar,rho,dt,coeffs)
end

function P_RHS(derv::ComplexF64,X::Float64)::ComplexF64
        return -(derv+X)
end

function X_RHS(P::Float64)::Float64
        return P
end
function P_RHS_RK(derv::Number,X::Number)::Number
	return -(derv+X)	
end
function X_RHS_RK(P::Number)::Number
	return P
end
function RK4_intertwined(HamPar,rho,dt)	
	Pn = []
	Xn = []
	for i in 1:HamPar.N
		H_red = Hamiltonian.FranssonHamiltonian_Mode_v2(HamPar,i)
		derv  = tr(rho*H_red)
		p1 = dt*P_RHS_RK(derv,HamPar.X[i])	
		x1 = dt*p1
		p2 = dt*P_RHS_RK(derv,x1 + 1/3 * p1)
                x2 = dt*(p2 + 1/3 * x1)
		p3 = dt*P_RHS_RK(derv,x2 - 1/3 * p1 + p2)
		x3 = dt*(p3 - 1/3 * x1 + x2)
		p4 = dt*P_RHS_RK(derv,x3 + 1.0 * p1 - p2 + p3)	
		Pn_c = HamPar.P[i] + 1/8 * (p1 + 3*p2 + 3*p3 + p4)
                x4 = dt*(p4 + 1.0 * x1 - x2 + x3)
                Xn_c = HamPar.X[i] + 1/8 * (x1 + 3*x2 + 3*x3 + x4)
		push!(Pn,real(Pn_c))
		push!(Xn,real(Xn_c))
	end
	return Xn,Pn		
end
function RK4_sequential(HamPar::Any,rho::Array,dt::Float64)::Tuple
        Pn = []
        Xn = []
        for i in 1:HamPar.N
                H_red = FranssonHamiltonian_Mode_v2(HamPar,i)
                derv  = tr(rho*H_red)
                p1 = dt*P_RHS_RK(derv,HamPar.X[i])
                p2 = dt*P_RHS_RK(derv,HamPar.X[i] + 1/3 * p1)
                p3 = dt*P_RHS_RK(derv,HamPar.X[i] - 1/3 * p1 + p2)
                p4 = dt*P_RHS_RK(derv,HamPar.X[i] + 1.0 * p1 - p2 + p3)
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




