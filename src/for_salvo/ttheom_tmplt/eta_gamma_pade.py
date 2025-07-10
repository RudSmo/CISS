# This module calculates the coefficients and exponents in the decomposition of the bath correlation functions. 
# In order to be consistent with Yaling's code, we have deleted factors of 2*pi from the coefficients. 
# Also note that the dimensions of the coeffs/exponents in the wide-band case is one less than in the non wide-band case 
# (since there is no pole from the lorentzian). This is modified from Sam's original code in which the dimensions are the same.

import numpy as np

# Generate bath coefficients in the non wide-band case
def genbathcoeff(kappa,zeta,specwidth,muvec,temperature,N_leads,N_fermipoles,N_sign):
    
                                                                   
    gamma_vec = np.zeros((N_leads,N_sign,N_fermipoles+1),dtype=complex) # Initialize gamma_vec and eta_vec as 0 with double pre1jsion
    eta_vec = np.zeros((N_leads,N_sign,N_fermipoles+1),dtype=complex) 

    for itri in range(0,N_leads):                                                              # Loop through the leads index
        gamma_vec[itri,0,0] = specwidth[itri] - 1j*muvec[itri]                          # Generate the gamma exponent for this lead, the '0th' Pade pole, and sigma = -
        gamma_vec[itri,1,0] = specwidth[itri] + 1j*muvec[itri]                          # Generate the gamma exponent for this lead, the '0th' Pade pole, and sigma = -
        for itrj in range(1,N_fermipoles+1):                                                              # Loop through the remaining Pade poles
            gamma_vec[itri,0,itrj] = np.sqrt(-zeta[itrj-1])*temperature - 1j*muvec[itri]          # Generate the gamma exponent for this lead, the itrj-th Pade pole, and sigma = -
            gamma_vec[itri,1,itrj] = np.sqrt(-zeta[itrj-1])*temperature + 1j*muvec[itri]          # Generate the gamma exponent for this lead, the itrj-th Pade pole, and sigma = +
    
        
        F_upper = genfermiapp(1j*specwidth[itri]/temperature,kappa,zeta,N_fermipoles)                   # Calculate f_{approx}(iW)
        F_lower = genfermiapp(-1j*specwidth[itri]/temperature,kappa,zeta,N_fermipoles)                   # Calculate f_{approx}(iW)
        # do itrj = 0,(N_sign-1)                                                             # Loop through the electron levels
        eta_vec[itri,0,0] = 0.5*specwidth[itri]*(0.5 + F_upper)
        eta_vec[itri,1,0] = 0.5*specwidth[itri]*(0.5 - F_lower)
                                                                                    # For each electron level, generate the eta coefficient of the '0th' Pade pole
        for itrk in range(1,N_fermipoles+1):                                                          # Loop through the Pade poles
            eta_vec[itri,0,itrk] = -1j*temperature*kappa[itrk-1] \
                                    *(specwidth[itri]**2)/(zeta[itrk-1]*(temperature**2) + specwidth[itri]**2)   # Generate eta_{l}*V_{K,m} (K = itri,m = itrj,l = itrk)

            eta_vec[itri,1,itrk] = -1j*temperature*kappa[itrk-1] \
                                    *(specwidth[itri]**2)/(zeta[itrk-1]*(temperature**2) + specwidth[itri]**2)
          

    return eta_vec, gamma_vec 

# ---------------------------------------------------------------------------
# 
#     GENERATE COEFFICIENTS AND EXPONENTS OF BATH-CORRELATION EXPANSION
#                   UNDER THE WIDE-BAND APPROXIMATION
#
# ---------------------------------------------------------------------------
#
# Same as genbathcoeff, except under the wide-band approximation

def genbathcoeff_wbl(kappa,zeta,muvec,temperature,N_leads,N_fermipoles,N_sign):
        
    gamma_vec = np.zeros((N_leads,N_sign,N_fermipoles),dtype=complex)                                                                  # Initialize gamma_vec and eta_vec as 0 with double pre1jsion
    eta_vec = np.zeros((N_leads,N_sign,N_fermipoles),dtype=complex)                                                                  # Initialize gamma_vec and eta_vec as 0 with double pre1jsion
    

    for itri in range(0,N_leads):                                                              # Loop through the leads index
        for itrj in range(0,N_fermipoles):                                                            # Loop through the remaining Pade poles
            gamma_vec[itri,0,itrj] = np.sqrt(-zeta[itrj])*temperature - 1j*muvec[itri]            # Generate the gamma exponent for this lead, the itrj-th Pade pole, and sigma = -
            gamma_vec[itri,1,itrj] = np.sqrt(-zeta[itrj])*temperature + 1j*muvec[itri]            # Generate the gamma exponent for this lead, the itrj-th Pade pole, and sigma = +
        
        
        for itrj in range(0,N_sign):                                                           # Loop through the electron levels
            for itrk in range(0,N_fermipoles):                                                        # Loop through the Pade poles)
                                                                                                   # Generate eta_{l}*V_{K,m} (K = itri,m = itrj,l = itrk)
                eta_vec[itri,itrj,itrk] = -1j*temperature*kappa[itrk] 
    
    return eta_vec, gamma_vec 


# ---------------------------------------------------------------------------
#
#     APPROXIMATE FERMI-DIRAC FUNCTION BASED ON PADE DECOMPOSITION
#
# ---------------------------------------------------------------------------
#
# This FORTRAN subroutine is meant to be converted to a Python wrapper using f2py, for use 
# in the main Python code. It takes the system parameters and results of the Pade approximation
# and outputs the exponents and coeffi1jents of the power series expansion of the bath-correlation functions,
# which are calculated from residue theory.
#
# USAGE (IN FORTRAN):
#       call genfermiapp(w,kappa,zeta,F,N_fermipoles)


def genfermiapp(w,kappa,zeta,N_fermipoles):



    print('zeta=')
    print(zeta)
    print('kapppa=')
    print(kappa)

    # x = w/temperature
    F = 0.
    for itri in range(0,N_fermipoles):                            # Implement Pade approximation
        F = F - (2*kappa[itri]*w)/(-zeta[itri] + w**2)
    
    return F 

