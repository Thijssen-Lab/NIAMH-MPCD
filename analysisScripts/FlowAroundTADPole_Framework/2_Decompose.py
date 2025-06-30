from scipy.special import eval_legendre
from scipy import stats
import json
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
import math
from scipy.optimize import curve_fit
from scipy.interpolate import griddata
from scipy.interpolate import interpn
from numba import jit

from scipy import integrate

from tqdm import tqdm

from matplotlib.colors import LogNorm
from matplotlib import cm
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

plt.style.use('shendrukGroupStyle')
import shendrukGroupFormat as ed

import scipy
import os

# Loading the in-silico experiment, then optional plotting 

def find_files(root_dir, target_file):
    # List to store the file paths
    file_paths = []

    # Walk through the directory
    for dirpath, _, filenames in os.walk(root_dir):
        if target_file in filenames:
            # Construct full file path
            file_path = os.path.join(dirpath, target_file)
            # file_paths.append(file_path)
            file_paths.append(dirpath)

    return file_paths

# Specify the root directory, target file name, and output file name
root_directory = "../../CLUSTER/PhDiaFric/"  # Replace with the root directory path

target_filename = "swimmerflowfield.dat"
output_filename = "found_paths.txt"

# Call the function
filepath=find_files(root_directory, target_filename)

for files in tqdm(filepath):
    xyU=np.load(files+'/xyU.npy')

    InputFile=files+'/../input.json'
    with open(InputFile, 'r') as InputF:
        tion=json.load(InputF)
    L=tion['domain']

    if len(L)==3:
        dim=3
        Lx=L[0]
        Ly=int(L[1]/np.sqrt(2))

    if len(L)==2:
        dim=2
        Lx=L[0]
        Ly=int(L[1]/2)


    # Construct Ur (different from xyU2, dot product)
    # Also creating data structures for the interpolation
    Ur=[]
    grid=[]

    for i in range(Lx):
        for j in range(Ly):
            x=float(i-L[0]/2)+0.5
            y=1*(float(j)+0.5)

            Ur.append((xyU[i,j,0]*x+xyU[i,j,1]*y)/np.sqrt(x*x+y*y))
            grid.append([x,y])
        
    for i in range(Lx):
        for j in range(Ly):
            x=float(i-L[0]/2)+0.5
            y=-1*(float(j)+0.5)

            Ur.append((xyU[i,j,0]*x-xyU[i,j,1]*y)/np.sqrt(x*x+y*y))
            grid.append([x,y])

    a=1
    n=5

    # Creating the array of radii that will be considered (staying close to the swimmer for decent stats)
    rl=[0.5+i*a for i in range(int(np.ceil(Lx/2)/a)-1)]

    # The data structure for each U_{r,n}(r)
    url=np.zeros((len(rl),n))

    er=-1
    d=Ly

    @jit
    def Interp(x,y):
        return griddata(grid,Ur,(x,y),method='cubic')

    N=300

    dt=np.pi/N
    Th=np.zeros((N,n))

    for i in range(N):
        for k in range(n):
            if dim==2:
                Th[i,k]=np.cos(dt*k*i)
            if dim==3:
                Th[i,k]=np.sin(dt*i)*eval_legendre(k,np.cos(dt*i))

    # Loop over radii
    for r in rl:

        er+=1
        
        dt=np.pi/N
        th_list=np.zeros((N,n))
        tint=np.zeros(N)

        X=[]
        Y=[]

        # Loop over polar shells
        for ti in (range(N)):
            th=ti*dt
            
            tint[ti]=th
            # Finding the x,y position within the loop 
            c=np.cos(th)
            s=np.sin(th)
            x=r*c
            y=r*s
            X.append(x)
            Y.append(y)
            A=1

        for k in (range(n)):
            th_list[:,k]=griddata(grid,Ur,(X,Y),method='cubic')
            th_list[:,k]*=Th[:,k]
            url[er,k]=integrate.simpson(th_list[:,k],tint,dt)
            url[er,k]*=(0.5*(2.0*k+1.0))


    # Plotting the result
    for k in range(n):
        plt.plot(rl,abs(url[:,k]),label=k)
    plt.xscale('log',base=2)
    plt.yscale('log')
    plt.legend()
    plt.savefig(files+'/Moments.png')
    # plt.show()

    np.save(files+'/url',url)

    for k in range(n):
        np.save(files+'/url'+str(k),url[:,k])