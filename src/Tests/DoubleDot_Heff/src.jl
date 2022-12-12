using LinearAlgebra
using PyPlot

function Hel(t0,tph,w0,X,gL,gR)
	return [-1im*gL/2 -t0-tph*sqrt(2*w0)*X; -t0-tph*sqrt(2*w0)*X -1im*gR/2]
end

function VNeumann(rho,H)
	return -1im * (H*rho-rho*H)
end

function Xfunc(w0,P)
	return w0 * P
end

function Pfunc(w0,X,tph,rho)
	return - (-tph*sqrt(2*w0)*(rho[1,2] + rho[2,1]) + w0^2*X)
end

function Simulate(rho0,x0,p0,gL,gR,w0,t0,tph,dt,numsteps)
	p1l = [rho0[1,1]]
	p2l = [rho0[2,2]]
	xl = [x0]
	pl = [p0]
	ti = [0.0]
	H = Hel(t0,tph,w0,x0,gL,gR)
	for i in 1:numsteps
		vn1	= dt * VNeumann(rho0,H)
		vn2	= dt * VNeumann(rho0+1/3*vn1,H)
		vn3	= dt * VNeumann(rho0-1/3*vn1+vn2,H)
		vn4	= dt * VNeumann(rho0+vn1 - vn2 + vn3,H)	
		rho_c	= rho0 + 1/8 * (vn1 + 3*vn2 + 3*vn3 + vn4)	
		
		
		p_c = p0 + dt * Pfunc(w0,x0,tph,rho_c)
		x_c = x0 + dt * Xfunc(w0,p_c)
		
		push!(p1l,real(rho_c[1,1]))
		push!(p2l,real(rho_c[2,2]))
		push!(xl,real(x_c))
		push!(pl,real(p_c))
		push!(ti,i*dt)

		rho0 = rho_c
		x0 = x_c
		p0 = p_c
		H = Hel(t0,tph,w0,x_c,gL,gR)
	end
	return ti, p1l,p2l,xl,pl
end

rho0 = [1.0 0; 0 0.0]
t0 = 1.0
X0 = 1.0
P0 = 0.0
gL = 0.25
gR = 0.25
w0 = 0.1
tph = 0.0
dt = 0.0001
numsteps = Int(500/dt)


ti, p1l,p2l,xl,pl = Simulate(rho0,X0,P0,gL,gR,w0,t0,tph,dt,numsteps)

obs_str = ["electrons","phonons"]
obs_val = [[p1l,p2l],[xl,pl]]
folder = "data"
if !isdir(folder)
	mkdir(folder)
end
filname = "L2_t0$(t0)_gL$(gL)_gR$(gR)_w0_$(w0)_tph$(tph)"


for (ostr,oval) in zip(obs_str,obs_val)
	fil = open("$(folder)/$(filname)_$(ostr).dat","w")
	for k in 1:100:length(ti)
		write(fil,"$(ti[k]) ")
		for j in 1:length(oval)
			write(fil, "$(oval[j][k]) ")
		end
		write(fil, "\n")
	end
	close(fil)
end
	
