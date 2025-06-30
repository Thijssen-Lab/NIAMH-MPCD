from scipy.special import eval_legendre
from scipy import stats
import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from tqdm import tqdm
import os

# Fetches file with a given name 
def find_files(root_dir, target_file):
    # List to store the file paths
    file_paths = []

    # Walk through the directory
    for dirpath, _, filenames in os.walk(root_dir):
        if target_file in filenames:
            # Construct full file path
            file_paths.append(dirpath)

    return file_paths

# Specify the root directory, target file name, and output file name
root_directory = "PhDiaFric/"  # Replace with the root directory path
target_filename = "swimmerflowfield.dat"

# Call the function
filepath=find_files(root_directory, target_filename)

# First loop reads data and averages
for files in tqdm(filepath):

    cnt=0

    # Extracts simulation input parameters using the json file 
    InputFile=files+'/../input.json'
    with open(InputFile, 'r') as InputF:
        tion=json.load(InputF)

    L=tion['domain']
    dt=int(tion.get('dt',10))
    t=int(tion['simSteps']/dt)
    swt=tion["swFlowOut"]/dt


    if len(L)==2:
        Lx=tion['domain'][0]
        Ly=tion['domain'][1]
        Lz=1

    if len(L)==3:
        Lx=tion['domain'][0]
        Ly=tion['domain'][1]
        Lz=tion['domain'][2]

    Size=Lx*Ly*Lz

    X=np.zeros(Lx)
    Y=np.zeros(Ly)
    Z=np.zeros(Lz)

    vRav=np.zeros((Lx,Ly,Lz,3))

    file=files+'/swimmerflowfield.dat'

    infile = open(file,"r")
    # print( '\tToss header ...' )
    for i in range(13):
        #toss header
        line = infile.readline()
        #print line

    # print( '\tRead data ...' )
    i=0
    j=0
    n=-1
    while infile:
        
        line = infile.readline()
        if not line:
            break
        else:
            split = line.split("\t")
        
        if len(split) == 7:
            time,Qx,Qy,Qz,Vx,Vy,Vz = split # to take into account old vs new format
        else:
            break # error handling

        j=int(i%(Size))
        i=i+1

        kx=int((j//Lz)//Ly)
        ky=int((j//Lz)%Ly)
        kz=int(j%Lz)
        
        Qx=float(Qx);Qy=float(Qy);Qz=float(Qz)
        if Vx!='' and Vy!='' and Vz!='':
            Vx=float(Vx);Vy=float(Vy);Vz=float(Vz)
        else:
            Vx=0.0;Vy=0.0;Vz=0.0
        

        dx=Qx-Lx/2+0.5
        dy=Qy-Ly/2+0.5
        dz=Qz-Lz/2+0.5
        cnt+=1

        X[kx]=dx
        Y[ky]=dy
        Z[kz]=dz

        vRav[kx,ky,kz,0]+=Vx/swt
        vRav[kx,ky,kz,1]+=Vy/swt
        vRav[kx,ky,kz,2]+=Vz/swt

    np.save(files+'/x',X)
    np.save(files+'/y',Y)
    np.save(files+'/z',Z)
    np.save(files+'/v',vRav/cnt)

# Second loop performs an azimuthal average (in 2-D, this is just a folding-over of the top and the bottom half)
for files in tqdm(filepath):

    Exp=files
    X=np.load(Exp+'/x.npy')
    Y=np.load(Exp+'/y.npy')
    Z=np.load(Exp+'/z.npy')

    V=np.load(Exp+'/v.npy')

    InputFile=Exp+'/../input.json'
    with open(InputFile, 'r') as InputF:
        tion=json.load(InputF)

    L=tion['domain']

    test=[]

    if len(L)==3:
        diags=int(np.ceil(len(Y)/np.sqrt(2)))

    if len(L)==2:
        diags=int(np.ceil(len(Y)/2))

    xyU=np.zeros((len(X),diags,3))
    cnt=np.zeros((len(X),diags))

    ix=-1
    for x in X:
        ix+=1
        iy=-1
        for y in Y:
            iy+=1
            iz=-1
            for z in Z:
                iz+=1

                phi=0-np.arctan2(z,y)
                test.append(float(phi))
                y2=y*np.cos(phi)-z*np.sin(phi)
                z2=y*np.sin(phi)+z*np.cos(phi)

                y2=int(y2+len(Y)/2)
                z2=int(z2+len(Z)/2)

                while y2>=len(Y):
                    y2-=len(Y)
                while y2<0:
                    y2+=len(Y)
                
                y2-=int(len(Y)/2)

                while y2<0:
                    y2+=len(Y)
            
                xyU[ix,y2,0]+=V[ix,iy,iz,0]
                xyU[ix,y2,1]+=V[ix,iy,iz,1]*np.cos(phi)-V[ix,iy,iz,2]*np.sin(phi)
                xyU[ix,y2,2]+=V[ix,iy,iz,1]*np.sin(phi)+V[ix,iy,iz,2]*np.cos(phi)
                cnt[ix,y2]+=1

    offset=-tion['sigSwim']*tion['dsSwim']
    offset=-tion['sigSwim']-tion['dsSwim']
    offset=int(tion['sigSwim']*(1-2*tion['dsSwim'])/4)

    a=xyU[:,:,0]/cnt
    b=xyU[:,:,1]/cnt
    c=xyU[:,:,2]/cnt

    for i in range(len(a[:,0])):
        j=i+offset
        if j<0:
            j+=L[0]
        xyU[j,:,0]=a[i,:]
        xyU[j,:,1]=b[i,:]
        xyU[j,:,2]=c[i,:]

    # Optional plotting as a sanity check 

    # from matplotlib.colors import LogNorm

    # # plt.imshow(xyU[:,:,0].T,norm=LogNorm())
    # # plt.colorbar()
    # # plt.show()

    # from matplotlib.patches import Circle

    # # plt.imshow(xyU[:,:,0])
    # # plt.colorbar()
    # # plt.show()
    # # plt.imshow(xyU[:,:,1])
    # # plt.colorbar()
    # # plt.show()
    np.save(Exp+'/xyU',xyU)

# Averaging every repeat 
xyUtot=xyU
cnt=0
for files in tqdm(filepath):
    Exp=files

    xyU=np.load(Exp+'/xyU.npy')
    xyUtot+=xyU
    cnt+=1

xyUtot/=cnt

np.save(Exp+'/../xyU',xyUtot)
