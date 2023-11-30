

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

function SampleWigner_Multi(beta::Float64,w::Array{Float64},m::Int64,num_avg::Int64)::Tuple{Array{Float64},Array{Float64}}
	Xn = zeros(num_avg,m)
	Pn = zeros(num_avg,m)
	for j in 1:num_avg
		for i in 1:m
			W = SampleWigner(beta,w[i])
			Xn[j,i] = W[1]
			Pn[j,i] = W[2]
		end
	end
	return Xn, Pn	
end

function SampleGS_Multi(m::Int64,num_avg::Int64)::Tuple{Array{Float64},Array{Float64}}
        Xn = zeros(num_avg,m)
        Pn = zeros(num_avg,m)
        for j in 1:num_avg
                for i in 1:m
                        W = SampleGS()
                        Xn[j,i] = W[1]
                        Pn[j,i] = W[2]
                end
        end
        return Xn, Pn
end






