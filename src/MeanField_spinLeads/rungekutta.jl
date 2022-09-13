module RungeKutta
include("electrons.jl")
using .Electrons 

export RK_Currents, RK_Vneumann, RK_VneumannL

function RK_Currents(beta,hcen_i,Gamma1u,Gamma1d,Gamma2u,Gamma2d,GammaT,kappa,chi1u,chi1d,chi2u,chi2d,curp1u_i,curp1d_i,curp2u_i,curp2d_i,rho_i,dt)
	betaLu,betaLd,betaRu,betaRd = beta
	cu1u	= [dt * CurrentEOM(betaLu,Gamma1u,hcen_i,GammaT,kappa[j],chi1u[j],curp1u_i[j]) for j in 1:length(chi1u)]
	cu2u 	= [dt * CurrentEOM(betaLu,Gamma1u,hcen_i,GammaT,kappa[j],chi1u[j],curp1u_i[j]+1/3*cu1u[j]) for j in 1:length(chi1u)]
	cu3u 	= [dt * CurrentEOM(betaLu,Gamma1u,hcen_i,GammaT,kappa[j],chi1u[j],curp1u_i[j]-1/3*cu1u[j]+cu2u[j]) for j in 1:length(chi1u)]
	cu4u 	= [dt * CurrentEOM(betaLu,Gamma1u,hcen_i,GammaT,kappa[j],chi1u[j],curp1u_i[j]+cu1u[j]-cu2u[j]+cu3u[j]) for j in 1:length(chi1u)]
	curp1u_c= [curp1u_i[j] + 1/8 * (cu1u[j] +3*cu2u[j] + 3*cu3u[j] + cu4u[j]) for j in 1:length(chi1u)]
	cur1u_c	= CurrentMatrix(rho_i,Gamma1u,curp1u_c)

	cu1d	= [dt * CurrentEOM(betaLd,Gamma1d,hcen_i,GammaT,kappa[j],chi1d[j],curp1d_i[j]) for j in 1:length(chi1d)]
	cu2d 	= [dt * CurrentEOM(betaLd,Gamma1d,hcen_i,GammaT,kappa[j],chi1d[j],curp1d_i[j]+1/3*cu1d[j]) for j in 1:length(chi1d)]
	cu3d 	= [dt * CurrentEOM(betaLd,Gamma1d,hcen_i,GammaT,kappa[j],chi1d[j],curp1d_i[j]-1/3*cu1d[j]+cu2d[j]) for j in 1:length(chi1d)]
	cu4d 	= [dt * CurrentEOM(betaLd,Gamma1d,hcen_i,GammaT,kappa[j],chi1d[j],curp1d_i[j]+cu1d[j]-cu2d[j]+cu3d[j]) for j in 1:length(chi1d)]
	curp1d_c= [curp1d_i[j] + 1/8 * (cu1d[j] +3*cu2d[j] + 3*cu3d[j] + cu4d[j]) for j in 1:length(chi1d)]
	cur1d_c	= CurrentMatrix(rho_i,Gamma1d,curp1d_c)
	

	cu21u	= [dt * CurrentEOM(betaRu,Gamma2u,hcen_i,GammaT,kappa[j],chi2u[j],curp2u_i[j]) for j in 1:length(chi2u)]
	cu22u 	= [dt * CurrentEOM(betaRu,Gamma2u,hcen_i,GammaT,kappa[j],chi2u[j],curp2u_i[j]+1/3*cu21u[j]) for j in 1:length(chi2u)]
	cu23u 	= [dt * CurrentEOM(betaRu,Gamma2u,hcen_i,GammaT,kappa[j],chi2u[j],curp2u_i[j]-1/3*cu21u[j]+cu22u[j]) for j in 1:length(chi2u)]
	cu24u 	= [dt * CurrentEOM(betaRu,Gamma2u,hcen_i,GammaT,kappa[j],chi2u[j],curp2u_i[j]+cu21u[j]-cu22u[j]+cu23u[j]) for j in 1:length(chi2u)]	
	curp2u_c= [curp2u_i[j] + 1/8 * (cu21u[j] +3*cu22u[j] + 3*cu23u[j] + cu24u[j]) for j in 1:length(chi2u)]
	cur2u_c	= CurrentMatrix(rho_i,Gamma2u,curp2u_c)

	cu21d	= [dt * CurrentEOM(betaRd,Gamma2d,hcen_i,GammaT,kappa[j],chi2d[j],curp2d_i[j]) for j in 1:length(chi2d)]
	cu22d 	= [dt * CurrentEOM(betaRd,Gamma2d,hcen_i,GammaT,kappa[j],chi2d[j],curp2d_i[j]+1/3*cu21d[j]) for j in 1:length(chi2d)]
	cu23d 	= [dt * CurrentEOM(betaRd,Gamma2d,hcen_i,GammaT,kappa[j],chi2d[j],curp2d_i[j]-1/3*cu21d[j]+cu22d[j]) for j in 1:length(chi2d)]
	cu24d 	= [dt * CurrentEOM(betaRd,Gamma2d,hcen_i,GammaT,kappa[j],chi2d[j],curp2d_i[j]+cu21d[j]-cu22d[j]+cu23d[j]) for j in 1:length(chi2d)]
	curp2d_c= [curp2d_i[j] + 1/8 * (cu21d[j] +3*cu22d[j] + 3*cu23d[j] + cu24d[j]) for j in 1:length(chi2d)]
	cur2d_c	= CurrentMatrix(rho_i,Gamma2d,curp2d_c)

	cur_c	= cur1u_c + cur1d_c + cur2u_c + cur2d_c
	return curp1u_c,curp1d_c,curp2u_c,curp2d_c,cur1u_c,cur1d_c,cur2u_c,cur2d_c,cur_c
end
function RK_Vneumann(rho_i,hcen_i,dt)
	vn1	= dt * VNeumann(rho_i,hcen_i)
	vn2	= dt * VNeumann(rho_i+1/3*vn1,hcen_i)
	vn3	= dt * VNeumann(rho_i-1/3*vn1+vn2,hcen_i)
	vn4	= dt * VNeumann(rho_i+vn1 - vn2 + vn3,hcen_i)	
	rho_c	= rho_i + 1/8 * (vn1 + 3*vn2 + 3*vn3 + vn4)
	return rho_c
end


function RK_VneumannL(rho_i,hcen_i,cur_c,dt)
	vn1	= dt * VNeumannL(rho_i,hcen_i,cur_c)
	vn2	= dt * VNeumannL(rho_i+1/3*vn1,hcen_i,cur_c)
	vn3	= dt * VNeumannL(rho_i-1/3*vn1+vn2,hcen_i,cur_c)
	vn4	= dt * VNeumannL(rho_i+vn1 - vn2 + vn3,hcen_i,cur_c)
	rho_c	= rho_i + 1/8 * (vn1 + 3*vn2 + 3*vn3 + vn4)
	return rho_c

end


end
