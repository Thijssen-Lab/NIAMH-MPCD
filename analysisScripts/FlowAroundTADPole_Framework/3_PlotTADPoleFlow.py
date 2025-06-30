from scipy.special import eval_legendre
import json
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd

from scipy.interpolate import griddata
from numba import jit

from scipy import integrate
from tqdm import tqdm

from matplotlib.colors import LogNorm

plt.style.use('shendrukGroupStyle')
import shendrukGroupFormat as ed

import scipy

# Loading the in-silico experiment to be plotted

plt.figure(figsize=(10, 6), dpi=80)

Exp='../../CLUSTER/zeroNetMom/0.0001/1/700/1'

xyU=np.load(Exp+'/xyU.npy')

InputFile=Exp+'/../input.json'
with open(InputFile, 'r') as InputF:
    tion=json.load(InputF)
L=tion['domain']
f=tion['species'][0].get('damp',0)

for i in range(len(xyU[:,0,0])):
        for j in range(int(len(xyU[0,:,0]))):
            if (i-int(L[0]/2))*(i-int(L[0]/2))+j*j>=L[0]*L[0]/4:
                xyU[i,j,0]=0
                xyU[i,j,1]=0

xyU=xyU[:,0:int(L[1]/2),:]

xyU2=np.sqrt(xyU[:,:,0]*xyU[:,:,0]+xyU[:,:,1]*xyU[:,:,1]).T


u,v=np.meshgrid(np.arange(len(xyU[:,0,0])),np.arange(len(xyU[0,:,0])))
plt.xlabel('Distance along swimmer axis, $r$',fontsize=32)
plt.ylabel(r'Distance away\\ from swimmer axis',fontsize=32)

skip=1

w2=np.sqrt(xyU[:,:,0]*xyU[:,:,0]+xyU[:,:,1]*xyU[:,:,1])
w=[min(i[0],np.power(10.0,-15)) for i in w2]

a2=plt.imshow(xyU2,norm=LogNorm(),label='Vel')
plt.quiver(u[::skip,::skip],v[::skip,::skip],xyU[::skip,::skip,0].T,-xyU[::skip,::skip,1].T,linewidth=w)

plt.show()