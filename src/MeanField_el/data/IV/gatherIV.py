import numpy as np
import fnmatch,os

VsFolders = np.array(fnmatch.filter(os.listdir("."),"V*"))
Vs = np.array([float(V.split("V")[1]) for V in VsFolders])

order = np.argsort(Vs)
VsFolders = VsFolders[order]
Vs = Vs[order]

IV = []

fil = open("IV.dat","w")
for v,vf in zip(Vs,VsFolders):
    dat = np.genfromtxt(str(vf)+"/"+np.array(fnmatch.filter(os.listdir(str(vf)),"*current*"))[0])
    fil.write(f"{v} ")
    for k in dat[-1,1::]:
        fil.write(f"{k} ")
    fil.write("\n")

fil.close()




