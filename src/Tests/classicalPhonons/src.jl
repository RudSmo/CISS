using LinearAlgebra


function SymplecticEulerStep(X0,P0,w0,dt)
	Pn = P0 - dt * w0^2 * X0
        Xn = X0 + dt * Pn
        return Xn,Pn
end

function StoermerVerletStep(X0,P0,w0,dt)
	Xn = X0 + dt * (P0 + dt/2 * w0^2 * X0)
	Pn = P0 - dt/2 * (w0^2 * X0 + w0^2 * Xn) 
	return Xn,Pn
end

function Propagate(X0,P0,w0,dt,numStep)
	Xi = X0
	Pi = P0
	for i in 1:numStep	
		Xn,Pn = SymplecticEulerStep(Xi,Pi,w0,dt)
#		Xn,Pn = StoermerVerletStep(Xi,Pi,w0,dt)
		Xi = Xn
		Pi = Pn
		Etot = 0.5 * w0^2 * Xi^2 + 0.5 * Pi^2
		println(Etot)
	end
end

X0 = 0.4
P0 = -0.01
w0 = 1.0
dt = 0.0001
numStep = 10000
Propagate(X0,P0,w0,dt,numStep)
