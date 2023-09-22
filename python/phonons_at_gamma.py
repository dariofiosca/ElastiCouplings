import numpy as np
from MagInt.utils import print_arr
from math import *


f_meV=16.159 # convertion eigenvalues to meV for oxygen mass
f_cm1=130.34 # convertion eigenvalues to cm-1 for oxygen mass
#f_cm1=15.633302*15.633302

#q1_on_site=34.468798345807365
t2g_on_site=10.11186654199756 # on-site in eV/Ang^2
eg_on_site=20.306082230657907
N=5
NR=12

arr=np.loadtxt('Elastic_couplings.dat')
J0=np.zeros((N,N))
for iR in range(NR):
    J0+=arr[iR*N:iR*N+N,:]

print_arr(J0,log='Q = 0 elastic meanfield:')

on_site=np.zeros((N,N))

#on_site[0,0]=q1_on_site
for i in [0,1,3]: on_site[i,i]=t2g_on_site
for i in [2,4]: on_site[i,i]=eg_on_site

H=on_site + J0
#H=on_site

e,v=np.linalg.eigh(H)

print("\nPhonons at Gamma in cm-1:")
for i in range(N):
    print("%s %10.0f"%(i+1,np.sqrt(e[i])*f_cm1))


