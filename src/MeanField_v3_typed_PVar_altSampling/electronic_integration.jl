include("observables.jl")
using LinearAlgebra

function RK4_VNL!(shx	::Int64,
		 shy	::Int64,
		 nP	::Int64,
		 beta	::Float64,
		 hcen_i	::Array{ComplexF64},
		 Gamma1	::Array{Float64},
		 Gamma2	::Array{Float64},
		 GammaT	::Array{Float64},
		 kappa	::Array{Float64},
		 chi1	::Array{ComplexF64},
		 chi2	::Array{ComplexF64},
		 curp1_i::Array{ComplexF64},
		 curp2_i::Array{ComplexF64},
		 cur_i	::Array{ComplexF64},
		 rho_i	::Array{ComplexF64},
		 dt	::Float64)::Tuple{Array{ComplexF64}, Array{ComplexF64}, Array{ComplexF64}, Array{ComplexF64}, Array{ComplexF64}, Array{ComplexF64}}
#::Tuple{Array{ComplexF64}, Array{ComplexF64}}#, Array{ComplexF64}, Array{ComplexF64}, Array{ComplexF64}, Array{ComplexF64}}

	
	idt::ComplexF64 = -1im*dt		

	id  	= Matrix(I,size(rho_i))
	rhoA	= 1/4 .* (id .- 2 .* rho_i)
	vnl_b   = cur_i .+ conj(transpose(cur_i))		

	# SLOWEST and MOST MEMORY USAGE	
#	vn1	= dt .* VNeumannL_v2(rho_i,hcen_i,vnl_b)
#	vn2     = dt .* VNeumannL_v2(rho_i .+ 1/3 .* vn1,hcen_i,vnl_b)
#	vn3     = dt .* VNeumannL_v2(rho_i .- 1/3 .* vn1 .+ vn2,hcen_i,vnl_b)
#	vn4     = dt .* VNeumannL_v2(rho_i .+ vn1 .- vn2 .+ vn3,hcen_i,vnl_b)

	# FASTEST AND SMALLEST MEMORY USAGE
	vn1	= dt .* (-1im .* (hcen_i * rho_i .- rho_i * hcen_i) .+ vnl_b)
	hv1	= hcen_i * vn1 .- vn1 * hcen_i
	vn2	= vn1 .+ idt/3 .* hv1
	hv2 	= hcen_i * vn2 .- vn2 * hcen_i    
	vn3	= vn1 .- idt/3 .* hv1 .+ idt .* hv2
	vn4	= vn1 .+ idt .* hv1 .- idt .* hv2 .+ idt * (hcen_i * vn3 .- vn3 * hcen_i)	

	# FASTER THAN FIRST, BUT MOST MEMORY EXPENSIVE 
#	vn1     = dt .* (-1im .* (hcen_i * rho_i .- rho_i * hcen_i) .+ vnl_b)
#	vn2     = vn1 .+ idt/3 .* (hcen_i * vn1 .- vn1 * hcen_i)
#	vn3     = vn1 .- idt/3 .* (hcen_i * vn1 .- vn1 * hcen_i) .+ idt .* (hcen_i * vn2 .- vn2 * hcen_i)
#	vn4     = vn1 .+ idt .* (hcen_i * vn1 .- vn1 * hcen_i) .- idt .* (hcen_i * vn2 .- vn2 * hcen_i) .+ idt * (hcen_i * vn3 .- vn3 * hcen_i)
	
	rho_i	.+=  1/8 .* (vn1 .+ 3 .* vn2 .+ 3 .* vn3 .+ vn4)

#	curp1_c = similar(curp1_i, ComplexF64)
#	curp2_c = similar(curp2_i, ComplexF64)
			
	curp1_c_sum = zeros(ComplexF64,size(cur_i))
	curp2_c_sum = zeros(ComplexF64,size(cur_i))

	@inbounds @simd for j in 1:nP
		curp0 = @view curp1_i[:,:,j]
		curpA,curpB = CurrentEOM_v2(beta,Gamma1,hcen_i,GammaT,kappa[j],chi1[j]) 

		cu1     = idt .* (curpA .+ curpB * curp0)
		curpBcu1= curpB * cu1
		cu2     = cu1 .+ 1/3*idt .* curpBcu1
		curpBcu2= curpB * cu2
		cu3     = cu1 .- 1/3*idt .* curpBcu1 .+ idt .* curpBcu2
		cu4     = cu1 .+ idt .* curpBcu1 .- idt .* curpBcu2 .+ idt .* curpB * cu3
		
		@views curp1_i[:,:,j] = curp0 .+ 1/8 .* (cu1 .+ 3 .* cu2 .+ 3 .* cu3 .+ cu4)
		curp1_c_sum .+= curp1_i[:,:,j]
	
		curp0 = @view curp2_i[:,:,j]
		curpA,curpB = CurrentEOM_v2(beta,Gamma2,hcen_i,GammaT,kappa[j],chi2[j])	
		
		cu1     = idt .* (curpA .+ curpB * curp0)
                curpBcu1= curpB * cu1
                cu2     = cu1 .+ 1/3*idt .* curpBcu1
                curpBcu2= curpB * cu2
                cu3     = cu1 .- 1/3*idt .* curpBcu1 .+ idt .* curpBcu2
                cu4     = cu1 .+ idt .* curpBcu1 .- idt .* curpBcu2 .+ idt .* curpB * cu3
		
		@views curp2_i[:,:,j] = curp0 .+ 1/8 .* (cu1 .+ 3 .* cu2 .+ 3 .* cu3 .+ cu4)
		curp2_c_sum .+= curp2_i[:,:,j]
	end
	
        @views cur1_c   = rhoA*Gamma1 .+ curp1_c_sum
        @views cur2_c   = rhoA*Gamma2 .+ curp2_c_sum
	cur_i	.= cur1_c .+ cur2_c
#	return cur1_c,cur2_c
	return curp1_c_sum,curp2_c_sum,cur1_c,cur2_c,cur_i,rho_i
end


function RK45_VNL_RK4_Curs(beta,hcen_i,Gamma1,Gamma2,GammaT,kappa,chi1,chi2,curp1_i,curp2_i,cur_i,rho_i,dt)
        curp1_c = []
        curp2_c = []
        for j in 1:length(chi1)
                cu1     = dt * CurrentEOM(beta,Gamma1,hcen_i,GammaT,kappa[j],chi1[j],curp1_i[j])
                cu2     = dt * CurrentEOM(beta,Gamma1,hcen_i,GammaT,kappa[j],chi1[j],curp1_i[j]+1/3*cu1)
                cu3     = dt * CurrentEOM(beta,Gamma1,hcen_i,GammaT,kappa[j],chi1[j],curp1_i[j]-1/3*cu1+cu2)
                cu4     = dt * CurrentEOM(beta,Gamma1,hcen_i,GammaT,kappa[j],chi1[j],curp1_i[j]+cu1-cu2+cu3)
                push!(curp1_c,curp1_i[j] + 1/8 * (cu1 +3*cu2 + 3*cu3 + cu4))

                cu21    = dt * CurrentEOM(beta,Gamma2,hcen_i,GammaT,kappa[j],chi2[j],curp2_i[j])
                cu22    = dt * CurrentEOM(beta,Gamma2,hcen_i,GammaT,kappa[j],chi2[j],curp2_i[j]+1/3*cu21)
                cu23    = dt * CurrentEOM(beta,Gamma2,hcen_i,GammaT,kappa[j],chi2[j],curp2_i[j]-1/3*cu21+cu22)
                cu24    = dt * CurrentEOM(beta,Gamma2,hcen_i,GammaT,kappa[j],chi2[j],curp2_i[j]+cu21-cu22+cu23)

                push!(curp2_c,curp2_i[j] + 1/8 * (cu21 +3*cu22 + 3*cu23 + cu24))
        end

        cur1_c  = CurrentMatrix(rho_i,Gamma1,curp1_c)
        cur2_c  = CurrentMatrix(rho_i,Gamma2,curp2_c)
        cur_c   = cur1_c + cur2_c

	TOL = 2*1e-5
        c  = [[1/4],
              [3/32,9/32],
              [1932/2197,-7200/2197,7296/2197],
              [439/216,-8,3680/513,-845/4104],
              [-8/27,2,-3544/2565,1859/4104,-11/40]]
        c4  = [25/216,0,1408/2565,2197/4101,-1/5,0]
        c5  = [16/135,0,6656/12825,28561/56430,-9/50,2/55]

        vn1     = dt * VNeumannL(rho_i,hcen_i,cur_i)
        vn2     = dt * VNeumannL(rho_i+c[1][1]*vn1,hcen_i,cur_i)
        vn3     = dt * VNeumannL(rho_i+c[2][1]*vn1+c[2][2]*vn2,hcen_i,cur_i)
        vn4     = dt * VNeumannL(rho_i+c[3][1]*vn1+c[3][2]*vn2+c[3][3]*vn3,hcen_i,cur_i)
	vn5     = dt * VNeumannL(rho_i+c[4][1]*vn1+c[4][2]*vn2+c[4][3]*vn3+c[4][4]*vn4,hcen_i,cur_i)
	vn6     = dt * VNeumannL(rho_i+c[5][1]*vn1+c[5][2]*vn2+c[5][3]*vn3+c[5][4]*vn4+c[5][5]*vn5,hcen_i,cur_i)

	rk_stp = [vn1,vn2,vn3,vn4,vn5,vn6]
        rho_c_4   = rho_i + sum([c4[i]*rk_stp[i] for i in 1:length(rk_stp)])
	rho_c_5   = rho_i + sum([c5[i]*rk_stp[i] for i in 1:length(rk_stp)])
	
	if norm(rho_c_4-rho_c_5) < TOL
		dt *= 2
	else
		dt /= 2
	end
        return curp1_c,curp2_c,cur1_c,cur2_c,cur_c,rho_c_5,dt
end


function RK4_VNL_RK45_Curs(beta,hcen_i,Gamma1,Gamma2,GammaT,kappa,chi1,chi2,curp1_i,curp2_i,cur_i,rho_i,dt)
	TOL = 2*1e-5
	c  = [[1/4],
	      [3/32,9/32],	
	      [1932/2197,-7200/2197,7296/2197],
	      [439/216,-8,3680/513,-845/4104],
	      [-8/27,2,-3544/2565,1859/4104,-11/40]]
	c4  = [25/216,0,1408/2565,2197/4101,-1/5,0]
	c5  = [16/135,0,6656/12825,28561/56430,-9/50,2/55]
	
	nP  = length(chi1)
	curp1_c_4 = []
	curp1_c_5 = []
	for j in 1:nP
	
	        cu1     = dt * CurrentEOM(beta,Gamma1,hcen_i,GammaT,kappa[j],chi1[j],curp1_i[j]) 
        	cu2     = dt * CurrentEOM(beta,Gamma1,hcen_i,GammaT,kappa[j],chi1[j],curp1_i[j]+c[1][1]*cu1) 
	        cu3     = dt * CurrentEOM(beta,Gamma1,hcen_i,GammaT,kappa[j],chi1[j],curp1_i[j]+c[2][1]*cu1+c[2][2]*cu2) 
        	cu4     = dt * CurrentEOM(beta,Gamma1,hcen_i,GammaT,kappa[j],chi1[j],curp1_i[j]+c[3][1]*cu1+c[3][2]*cu2+c[3][3]*cu3) 
	        cu5     = dt * CurrentEOM(beta,Gamma1,hcen_i,GammaT,kappa[j],chi1[j],curp1_i[j]+c[4][1]*cu1+c[4][2]*cu2+c[4][3]*cu3+c[4][4]*cu4) 
		cu6     = dt * CurrentEOM(beta,Gamma1,hcen_i,GammaT,kappa[j],chi1[j],curp1_i[j]+c[5][1]*cu1+c[5][2]*cu2+c[5][3]*cu3+c[5][4]*cu4+c[5][5]*cu5) 
		
        	push!(curp1_c_4, curp1_i[j] + c4[1]*cu1 +c4[2]*cu2 + c4[3]*cu3 + c4[4]*cu4 + c4[5]*cu5 + c4[6]*cu6)
		push!(curp1_c_5, curp1_i[j] + c5[1]*cu1 +c5[2]*cu2 + c5[3]*cu3 + c5[4]*cu4 + c5[5]*cu5 + c5[6]*cu6)
	end
        cur1_c_4  = CurrentMatrix(rho_i,Gamma1,curp1_c_4)
	cur1_c_5  = CurrentMatrix(rho_i,Gamma1,curp1_c_5)
	
	curp2_c_4 = []
	curp2_c_5 = []
	for j in 1:nP

		cu21    = dt * CurrentEOM(beta,Gamma2,hcen_i,GammaT,kappa[j],chi2[j],curp2_i[j])
        	cu22    = dt * CurrentEOM(beta,Gamma2,hcen_i,GammaT,kappa[j],chi2[j],curp2_i[j]+c[1][1]*cu21)
	        cu23    = dt * CurrentEOM(beta,Gamma2,hcen_i,GammaT,kappa[j],chi2[j],curp2_i[j]+c[2][1]*cu21+c[2][2]*cu22)
        	cu24    = dt * CurrentEOM(beta,Gamma2,hcen_i,GammaT,kappa[j],chi2[j],curp2_i[j]+c[3][1]*cu21+c[3][2]*cu22+c[3][3]*cu23)
	        cu25    = dt * CurrentEOM(beta,Gamma2,hcen_i,GammaT,kappa[j],chi2[j],curp2_i[j]+c[4][1]*cu21+c[4][2]*cu22+c[4][3]*cu23+c[4][4]*cu24)
        	cu26    = dt * CurrentEOM(beta,Gamma2,hcen_i,GammaT,kappa[j],chi2[j],curp2_i[j]+c[5][1]*cu21+c[5][2]*cu22+c[5][3]*cu23+c[5][4]*cu24+c[5][5]*cu25)

		push!(curp2_c_4,curp2_i[j] + c4[1]*cu21 + c4[2]*cu22 + c4[3]*cu23 + c4[4]*cu24 + c4[5]*cu25 + c4[6]*cu26)
		push!(curp2_c_5,curp2_i[j] + c5[1]*cu21 + c5[2]*cu22 + c5[3]*cu23 + c5[4]*cu24 + c5[5]*cu25 + c5[6]*cu26)
	end
	
        cur2_c_4  = CurrentMatrix(rho_i,Gamma2,curp2_c_4)
	cur2_c_5  = CurrentMatrix(rho_i,Gamma2,curp2_c_5)
        cur_c_4   = cur1_c_4 + cur2_c_4
	cur_c_5   = cur1_c_5 + cur2_c_5

	dif_cur = norm(cur_c_4-cur_c_5)
	s_cur   = 0.84 * (TOL * dt/dif_cur)^(1/4)
	dt_cur  = s_cur * dt
	
        vn1     = dt * VNeumannL(rho_i,hcen_i,cur_i)
        vn2     = dt * VNeumannL(rho_i+1/3*vn1,hcen_i,cur_i)
        vn3     = dt * VNeumannL(rho_i-1/3*vn1+vn2,hcen_i,cur_i)
        vn4     = dt * VNeumannL(rho_i+vn1 - vn2 + vn3,hcen_i,cur_i)
        rho_c   = rho_i + 1/8 * (vn1 + 3*vn2 + 3*vn3 + vn4)
	
        return curp1_c_5,curp2_c_5,cur1_c_5,cur2_c_5,cur_c_5,rho_c,dt_cur
end


