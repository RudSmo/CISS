using LinearAlgebra
using PyPlot

function Xfunc(w0,P)
	return w0 * P
end

function Pfunc(w0,X,tph,C)
	return - (-tph*sqrt(2*w0)*C + w0^2*X)
end

function Simulate(C,x0,p0,w0,tph,dt,numsteps)
	xl = [x0]
	pl = [p0]
	ti = [0.0]
	for i in 1:numsteps
		p_c = p0 + dt * Pfunc(w0,x0,tph,C)
		x_c = x0 + dt * Xfunc(w0,p_c)
		
		push!(xl,real(x_c))
		push!(pl,real(p_c))
		push!(ti,i*dt)

		x0 = x_c
		p0 = p_c
	end
	return ti, xl,pl
end

C = 0.3

w0 = 0.1
tph = 0.05
dt = 0.001
numsteps = Int(500/dt)

for X0 in -1:0.1:1
	for P0 in -1:0.1:1
		ti, xl,pl = Simulate(C,X0,P0,w0,tph,dt,numsteps)

		obs_str = ["phonons"]
		obs_val = [[xl,pl]]
		folder = "data"
		if !isdir(folder)
			mkdir(folder)
		end
		filname = "x$(X0)_p$(P0)_w0_$(w0)_tph$(tph)"


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
	end
end
	
