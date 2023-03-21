module RungeKutta
include("electrons.jl")
using .Electrons 

export RK_Auxiliary


function RK_Auxiliary(nAll,H,ChiPade,EtaPade,BetaL,DeltaL,Gamma,Axi,Omega0,Psi0,Rho0,dt)
	# RK4 (3/8)
	# 1 | 1 
	# 2 | 1 1/3
 	# 3 | 1 -1/3 1
	# 4 | 1 1 -1 1 
	# ------------
	#     1 3 3 1 * 1/8
	omega_1 = dt * OmegaEOM(nAll,ChiPade,DeltaL,Omega0)
	omega_2 = dt * OmegaEOM(nAll,ChiPade,DeltaL,Omega0 + 1/3 * omega_1)
	omega_3 = dt * OmegaEOM(nAll,ChiPade,DeltaL,Omega0 - 1/3 * omega_1 + omega_2)	
	omega_4 = dt * OmegaEOM(nAll,ChiPade,DeltaL,Omega0 + omega_1 - omega_2 + omega_3)
	omega_c = Omega0 + 1/8 * (omega_1 + 3*omega_2 + 3*omega_3 + omega_4)
	
	psi_1   = dt * PsiEOM(nAll,H,DeltaL,BetaL,ChiPade,EtaPade,Gamma,Axi,omega_c,Psi0)	
	psi_2   = dt * PsiEOM(nAll,H,DeltaL,BetaL,ChiPade,EtaPade,Gamma,Axi,omega_c,Psi0 + 1/3 * psi_1)
	psi_3   = dt * PsiEOM(nAll,H,DeltaL,BetaL,ChiPade,EtaPade,Gamma,Axi,omega_c,Psi0 - 1/3 * psi_1 + psi_2)
	psi_4   = dt * PsiEOM(nAll,H,DeltaL,BetaL,ChiPade,EtaPade,Gamma,Axi,omega_c,Psi0 + psi_1 - psi_2 + psi_3)
	psi_c 	= Psi0 + 1/8 * (psi_1 + 3*psi_2 + 3*psi_3 + psi_4)
	
	pi_c 	= CurrentMatrix(nAll,psi_c,Axi)	

	rho_1   = dt * RhoEOM(nAll,H,pi_c,Rho0) 
	rho_2   = dt * RhoEOM(nAll,H,pi_c,Rho0 + 1/3 * rho_1)
	rho_3   = dt * RhoEOM(nAll,H,pi_c,Rho0 - 1/3 * rho_1 + rho_2)
	rho_4   = dt * RhoEOM(nAll,H,pi_c,Rho0 + rho_1 - rho_2 + rho_3)
	rho_c	= Rho0 + 1/8 * (rho_1 + 3*rho_2 + 3*rho_3 + rho_4) 
	
	return omega_c, psi_c, rho_c
end


end
