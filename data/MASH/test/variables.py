#!/usr/bin/python

from __future__ import print_function

from runelib import units
import numpy as np

S = 8
F = 28

with open('ortho-couplings-28-modes.inp') as f:
    dict_param = {}
    file = f.readlines()

    for row in file:
        row.splitlines()

        if row.startswith('#'):
            name = row.rstrip('\r\n')
            dict_param[name] = {}
        else:
            val   = row.split(',')[0]
            key   = val.split('=')[0]
            value = val.split('=')[1:]

            for i in value:
                num = float(i)
                dict_param[name][key] = num

"""Frequencies and Mass """
omega_list = np.zeros(F,dtype=np.float64)
mass_list  = np.zeros(F,dtype=np.float64)

for idx, value in dict_param['#frequencies'].items():
    ind = idx.strip(' ')
    for f in range(F):
        name = 'w' +str(f+1)
        if ind == name:
            omega_au = value*units.eV
            omega_list[f] = omega_au
            mass_list[f]  = 1./omega_au

#print(omega_list)

"""State Energies"""
stateE_list = np.zeros((S,S),dtype=np.float64)
stateE_sWAP = np.zeros((S,S),dtype=np.float64)
for s1 in range(S):
    for s2 in range(S):
        for idx, value in dict_param['#diagonal zeroth order'].items():
            ind = idx.strip(' ')
            name = 'e' +str(s1+1)+str(s1+1)
            if ind == name:
                stateE_list[s1,s1] = value*units.eV
        for idx,value in dict_param['#off_diagonal zeroth order'].items():
            ind = idx.strip(' ')
            name = 'e' +str(s1+1)+str(s2+1)
            if ind == name:
                stateE_list[s1,s2] = value*units.eV
                stateE_list[s2,s1] = stateE_list[s1,s2]

#print(stateE_list)
#
# stateE_sWAP[1,1] = stateE_list[3,3]
# stateE_sWAP[3,3] = stateE_list[1,1]
#
# stateE_sWAP[0,0] = stateE_list[0,0]
# stateE_sWAP[2,2] = stateE_list[2,2]
#
#
# stateE_sWAP[0,2] = stateE_list[0,2]
# stateE_sWAP[2,0] = stateE_sWAP[0,2]
#
#
# stateE_sWAP[1,2] = stateE_list[3,2]
# stateE_sWAP[3,2] = stateE_list[1,2]
# stateE_sWAP[2,1] = stateE_sWAP[1,2]
# stateE_sWAP[2,3] = stateE_sWAP[3,2]
#
# stateE_sWAP[1,0] = stateE_list[3,0]
# stateE_sWAP[0,1] = stateE_sWAP[1,0]
#
# stateE_sWAP[3,0] = stateE_list[0,1]
# stateE_sWAP[0,3] = stateE_sWAP[3,0]
#
# stateE_sWAP[3,1] = stateE_list[1,3]
# stateE_sWAP[1,3] = stateE_sWAP[3,1]
#
# #
# for i in range(4):
#  for j in range(4):
#
#      #stateE_sWAP[i,1] = stateE_list[i,3]
#      #stateE_sWAP[1,i] = stateE_list[i,3]
#
#      #stateE_sWAP[i,3] = stateE_list[i,1]
#      #stateE_sWAP[3,i] = stateE_list[i,1]
#      stateE_sWAP[i,j] = stateE_list[i,j]
#      stateE_sWAP[j,i] = stateE_sWAP[i,j]
#
#      stateE_sWAP[7,i] = stateE_list[7,i]
#      stateE_sWAP[i,7] = stateE_sWAP[7,i]
#      stateE_sWAP[6,i] = stateE_list[6,i]
#      stateE_sWAP[i,6] = stateE_sWAP[6,i]
#

##################################
#
#for i in range(4,8):
#    for j in range(4,8):
#        stateE_sWAP[i,3] = stateE_list[i,1]
#        stateE_sWAP[3,i] = stateE_sWAP[i,3]
#        stateE_sWAP[i,1] = stateE_list[i,3]
#        stateE_sWAP[1,i] = stateE_sWAP[i,1]
#        stateE_sWAP[i,j] = stateE_list[i,j]
#        stateE_sWAP[j,i] = stateE_list[i,j]
##
# #
# # stateE_sWAP[4,4] = stateE_list[5,5]
# #
# #stateE_sWAP[5,5] = stateE_list[4,4]
# #
# stateE_sWAP[4,0] = stateE_list[5,0]
# stateE_sWAP[0,4] = stateE_sWAP[4,0]
# stateE_sWAP[5,0] = stateE_list[0,4]
# stateE_sWAP[0,5] = stateE_sWAP[5,0]
# stateE_sWAP[4,1] = stateE_list[5,3]
# stateE_sWAP[1,4] = stateE_sWAP[4,1]
# stateE_sWAP[5,1] = stateE_list[4,3]
# stateE_sWAP[1,5] = stateE_sWAP[5,1]
# stateE_sWAP[2,5] = stateE_list[2,4]
# stateE_sWAP[5,2] = stateE_list[2,4]
# stateE_sWAP[2,4] = stateE_list[2,5]
# stateE_sWAP[4,2] = stateE_list[2,5]
# stateE_sWAP[4,3] = stateE_list[5,1]
# stateE_sWAP[3,4] = stateE_sWAP[4,3]
# stateE_sWAP[5,3] = stateE_list[4,1]
# stateE_sWAP[3,5] = stateE_sWAP[5,3]
# #
# #
# stateE_sWAP[5,4] = stateE_list[4,5]
# stateE_sWAP[4,5] = stateE_sWAP[5,4]
#
#
# # print(stateE_sWAP )

# for i in range(len(stateE_sWAP)):
#  for j in range(len(stateE_sWAP)):
#   if j != len(stateE_sWAP)-1:
#    print( "{} & ".format(stateE_sWAP[i, j]), end='')
#   else:
#    print( "{} \\\\".format(stateE_sWAP[i, j]) )
#

"""Linear Coupling """
linear_coupling = np.zeros((F,S,S),dtype=np.float64)

for l in range(F):
    for i in range(S):
        for idx, value in dict_param.items():
            if idx.startswith('#mode'):
                for ind, val in dict_param[idx].items():
                    key = ind.strip(' ')
                    name_ii = "H" +str(i+1)+str(i+1)+"_k"+str(l+1)

                    if key == name_ii :
                        linear_coupling[l,i,i] = val*units.eV
                    for j in range(S):
                        name_ij = "H" +str(i+1)+str(j+1)+"_k"+str(l+1)
                        if key == name_ij :
                            linear_coupling[l,i,j] = val*units.eV
                            linear_coupling[l,j,i] = linear_coupling[l,i,j]


"""Bilinear Coupling """
bilin_coupling = np.zeros((F,S,S),dtype=np.float64)

for l in range(F):
    for i in range(S):
        for idx, value in dict_param.items():
            if idx.startswith('#mode'):
                for ind, val in dict_param[idx].items():
                    key = ind.strip(' ')
                    name_ii = "H" +str(i+1)+str(i+1)+"_g"+str(l+1)
                    if key == name_ii :
                        bilin_coupling[l,i,i] = val*units.eV
                    for j in range(S):
                        name_ij = "H" +str(i+1)+str(j+1)+"_g"+str(l+1)
                        if key == name_ij :
                            bilin_coupling[l,i,j] = val*units.eV
                            bilin_coupling[l,j,i] = bilin_coupling[l,i,j]





"""Potential Energy Matrix (from paper)"""

# StateE_paper = np.zeros((S,S),dtype=np.float64)
# off_diag = np.zeros((S,S))

# diag = np.array([0, 1576, 1577,1595,1744,1744,1809,1821])


# off_diag[0,:] = [0,15,-15,-18,-108,107,-426, 416]
# off_diag[1,:] = [15,1576,23,7,-30,-43,52,-7]
# off_diag[2,:] = [-15,23,1577,-7, -44, -31, -9, 51]
# off_diag[3,:] = [-18,7,-7,1595,-273, 273, -26,24]
# off_diag[4,:] = [-108,-30,-44,-273,1744,26, 43, -36]
# off_diag[5,:] = [107,-43,-31,273,26,1744,-39, 42]
# off_diag[6,:] = [-426,52,-9,-26,43,-39,1809,176]
# off_diag[7,:] = [416, -7,51,24,-36,42,176,1821]

# for i in range(S):
#     for j in range(S):
#         StateE_paper[i,i] = float(diag[i]*10**(-3))*units.eV
#         a = off_diag[i,j]*10**(-3)

#         StateE_paper[i,j] = a*units.eV


#         StateE_paper[j,i] = StateE_paper[i,j]

#print(StateE_paper)
