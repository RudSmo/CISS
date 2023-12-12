using LinearAlgebra

function SampleElectrons(pops,gamma)
	# Sample electronic mapping variables from unit circle
	# pop == population of the state
	# gamma == ZPE parameter

	xl = []
	pl = []
	for pop in pops
		r = sqrt(2*(pop+gamma))
		phi = rand()*2*pi
		x = r*cos(phi)
		p = r*sin(phi)
		push!(xl,x)
		push!(pl,p)
	end

	return xl,pl
end
function SampleWigner(beta::Float64,w::Float64)::Tuple{Float64,Float64}
        a = tanh((beta/w)*w/2)
        u1,u2 = rand(2)
        z1 = 1/sqrt(a)*sqrt(-2*log(u1))*cos(2*pi*u2)
        z2 = 1/sqrt(a)*sqrt(-2*log(u1))*sin(2*pi*u2)
        return z1,z2
end
function SampleVibrations(beta,w,m)
	# Sample from thermal Wigner distribution, given w0, beta and number of vib modes
	Xn = zeros(m)
        Pn = zeros(m)
	a = tanh(beta*w/2)
        for i in 1:m
        	u1,u2 = rand(2)
	        z1 = 1/sqrt(a)*sqrt(-2*log(u1))*cos(2*pi*u2)
	        z2 = 1/sqrt(a)*sqrt(-2*log(u1))*sin(2*pi*u2)		
                Xn[i] = z1
                Pn[i] = z2
        end
        return Xn, Pn		
end

function VibrationalEOM(X,P,x,p,t1)
	# EOM for the vibrational degrees of freedoms
	N = length(X)
	X_l = zeros(N)
	P_l = zeros(N)
	for i in 1:N
		X_l[i] = P[i]
		P_l[i] = -X[i] - t1/sqrt(2) * (x[i]*x[i+1] + p[i]*p[i+1])
	end
	return X_l,P_l
end

function IntegrateVib(X,P,x,p,t1,dt)
	a = 1/(2-2^(1/3))
        coeffs = a*[[1/2,(1-2^(1/3))/2,(1-2^(1/3))/2,1/2],
                    [0.,1.,-(2)^(1/3),1.]]
	N = length(X)
        Pn = similar(P)
        Xn = similar(X)
	XEOM,PEOM = VibrationalEOM(X,P,x,p,t1)
        for i in 1:N
            Xn_c = X[i]
            Pn_c = P[i]
            for j in 1:length(coeffs[1])
                    Pn_c += -coeffs[2][j]*dt*(Xn_c + t1/(2*sqrt(2)) * (x[i]*x[i+1] + p[i]*p[i+1]))
                    Xn_c += coeffs[1][j]*dt*Pn_c
            end
            Pn[i] = Pn_c
            Xn[i] = Xn_c
        end
        return Xn,Pn
end


function ElectronicEOM(x,p,X, e0, t0, t1)
	# EOM for the electronic mapping variables
	N = length(x)
	x_l = zeros(N)	
	p_l = zeros(N)
	for i in 1:N
		if i == N
			x_l[i] = 2 * e0 * p[i] 
			p_l[i] = -2 * e0 * x[i]  
		else
			x_l[i] = 2 * e0 * p[i] + p[i+1] * (-t0 + t1*X[i]/sqrt(2))
			p_l[i] = -2 * e0 * x[i] - x[i+1] * (-t0 + t1*X[i]/sqrt(2))	
		end
	end
	return x_l,p_l	
end

function IntegrateEl(x,p,X,e0,t0,t1,dt)
        a = 1/(2-2^(1/3))
        coeffs = a*[[1/2,(1-2^(1/3))/2,(1-2^(1/3))/2,1/2],
                    [0.,1.,-(2)^(1/3),1.]]
        N = length(x)
        pn = similar(p)
        xn = similar(x)
        for i in 1:N-1
            xn_c = x[i]
            pn_c = p[i]
            for j in 1:length(coeffs[1])
                    pn_c += -coeffs[2][j]* dt * (e0 * xn_c + x[i+1]/2 * (-t0 + t1*X[i]/sqrt(2))) 
                    xn_c += coeffs[1][j] * dt * (e0 * pn_c + p[i+1]/2 * (-t0 + t1*X[i]/sqrt(2)))  
            end
            pn[i] = pn_c
            xn[i] = xn_c
        end
	xn_c = x[N]
        pn_c = p[N]
        for j in 1:length(coeffs[1])
                pn_c += -coeffs[2][j]* dt *  e0 * xn_c
                xn_c += coeffs[1][j] * dt *  e0 * pn_c
        end
        pn[N] = pn_c
        xn[N] = xn_c	
        return xn,pn
end

function Population_el(x,p,gamma)
	return real.([1/2*(xi^2 + pi^2 - 2*gamma) for (xi,pi) in zip(x,p)])
end

function Propagate(pop0,w0,beta,e0,t0,t1,gamma,dt,numstep)
	# Propagate and save stuff to file
	x0,p0 = SampleElectrons(pop0,gamma)
	X0,P0 = SampleVibrations(beta,w0,length(pop0)-1)
	for iN in 1:numstep

		Xc,Pc = IntegrateVib(X0,P0,x0,p0,t1,dt/4)
		xc,pc = IntegrateEl(x0,p0,Xc,e0,t0,t1,dt/2)
		Xc,Pc = IntegrateVib(Xc,Pc,xc,pc,t1,dt/4)	

		Xc,Pc = IntegrateVib(Xc,Pc,xc,pc,t1,dt/4)
		xc,pc = IntegrateEl(xc,pc,Xc,e0,t0,t1,dt/2)
		Xc,Pc = IntegrateVib(Xc,Pc,xc,pc,t1,dt/4)
		X0 = Xc
		P0 = Pc
		x0 = xc
		p0 = pc
		if iN == 1
			fil_pop = open("electronic.dat","w")
			fil_vib = open("vibrational.dat","w")
			write(fil_vib,"# \n")
			write(fil_pop,"# \n")
			close(fil_pop)
			close(fil_vib)
		end
		if iN % 20 == 1
			
			file_pop = open("electronic.dat","a")
			file_vib = open("vibrational.dat","a")
			

			write(file_pop,"$(iN*dt) ")
			write(file_vib,"$(iN*dt) ")
			pop = Population_el(xc,pc,gamma)
			for k in 1:length(pop)
				write(file_pop,"$(pop[k]) ")
			end
			for k in 1:length(Xc)
				write(file_vib,"$(Xc[k]) $(Pc[k]) ")
			end
			write(file_pop,"\n")
			write(file_vib,"\n")
			
			
			close(file_pop)
			close(file_vib)

		end
	end		
end

w0 = 0.002
beta = w0/0.025
pop0 = [1,0.0001,0.0001,0.0001]
pop0 = pop0/norm(pop0)
T0 = 0.04/w0
e0 = -6*T0
t0 = T0
t1 = 0.0 * T0
N = length(pop0)
gamma = (sqrt(N+1)-1) / N
dt = 0.00001
tend = 16
numstep = Int(round(tend/abs(dt)))
Propagate(pop0,w0,beta,e0,t0,t1,gamma,dt,numstep)
