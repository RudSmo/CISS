module RungeKutta
include("electrons.jl")
include("classical_spins.jl")
using .Electrons 
using .ClassicalSpins
using LinearAlgebra
export RK_Currents, RK_Vneumann, RK_VneumannL, RK_Spins

function RK_Currents(beta,hcen_i,Gamma1,Gamma2,GammaT,kappa,chi1,chi2,curp1_i,curp2_i,rho_i,dt)
	cu1	= [dt * CurrentEOM(beta,Gamma1,hcen_i,GammaT,kappa[j],chi1[j],curp1_i[j]) for j in 1:length(chi1)]
	cu2 	= [dt * CurrentEOM(beta,Gamma1,hcen_i,GammaT,kappa[j],chi1[j],curp1_i[j]+1/3*cu1[j]) for j in 1:length(chi1)]
	cu3 	= [dt * CurrentEOM(beta,Gamma1,hcen_i,GammaT,kappa[j],chi1[j],curp1_i[j]-1/3*cu1[j]+cu2[j]) for j in 1:length(chi1)]
	cu4 	= [dt * CurrentEOM(beta,Gamma1,hcen_i,GammaT,kappa[j],chi1[j],curp1_i[j]+cu1[j]-cu2[j]+cu3[j]) for j in 1:length(chi1)]
	curp1_c	= [curp1_i[j] + 1/8 * (cu1[j] +3*cu2[j] + 3*cu3[j] + cu4[j]) for j in 1:length(chi1)]#nPoles]
	cur1_c	= CurrentMatrix(rho_i,Gamma1,curp1_c)

	cu21	= [dt * CurrentEOM(beta,Gamma2,hcen_i,GammaT,kappa[j],chi2[j],curp2_i[j]) for j in 1:length(chi2)]
	cu22 	= [dt * CurrentEOM(beta,Gamma2,hcen_i,GammaT,kappa[j],chi2[j],curp2_i[j]+1/3*cu21[j]) for j in 1:length(chi2)]
	cu23 	= [dt * CurrentEOM(beta,Gamma2,hcen_i,GammaT,kappa[j],chi2[j],curp2_i[j]-1/3*cu21[j]+cu22[j]) for j in 1:length(chi2)]
	cu24 	= [dt * CurrentEOM(beta,Gamma2,hcen_i,GammaT,kappa[j],chi2[j],curp2_i[j]+cu21[j]-cu22[j]+cu23[j]) for j in 1:length(chi2)]

	curp2_c	= [curp2_i[j] + 1/8 * (cu21[j] +3*cu22[j] + 3*cu23[j] + cu24[j]) for j in 1:length(chi2)]#nPoles]
	cur2_c	= CurrentMatrix(rho_i,Gamma2,curp2_c)
	cur_c	= cur1_c + cur2_c
	return curp1_c,curp2_c,cur1_c,cur2_c,cur_c
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

function normMag(mag,nfac)
        normed = []
        for mg in mag
                nm = norm(mg)
                if nm != 0.0
                        push!(normed,nfac*mg/nm)
                else
                        push!(normed,[0.,0.,0.])
                end
        end
        return normed
end
function RK_Spins(HamPar,scd_i,dt)
	Jsd = HamPar.Jsd
	mag_i = HamPar.mag
	pos = HamPar.pos
	hext = HamPar.Hext
        mk_1    = SpinEOM(Jsd,scd_i,mag_i,pos,hext)
        mk_2    = SpinEOM(Jsd,scd_i,[mag_i[j]+(1.0/3)*dt*mk_1[j] for j in 1:length(mag_i)],pos,hext)
        mk_3    = SpinEOM(Jsd,scd_i,[mag_i[j]-(1.0/3)*dt*mk_1[j]+dt*mk_2[j] for j in 1:length(mag_i)],pos,hext)
        mk_4    = SpinEOM(Jsd,scd_i,[mag_i[j]+dt*(mk_1[j]-mk_2[j]+mk_3[j]) for j in 1:length(mag_i)],pos,hext)
        mag_p   = [mag_i[j] + (1.0/8)*dt*(mk_1[j]+3*mk_2[j]+3*mk_3[j]+mk_4[j]) for j in 1:length(mag_i)]
        mag_c   = normMag(mag_p,1.0)
        return mag_c
end


end
