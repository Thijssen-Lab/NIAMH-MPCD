import matplotlib
from matplotlib import pyplot as plt
from matplotlib import ticker
from pylab import *
from numpy import ma
from subprocess import call
from numpy import linalg as LA
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import os
import csv
import sys
import numpy as np
import math
import gyrationTensor_2D
###########################################################
### Calculates and plots end-to-end for 2D
###########################################################
# Use our custom style
plt.style.use('shendrukGroupStyle')
# Use our custom colours
import shendrukGroupFormat as ed
# Colour map to use
myMap=ed.plasma
FS=20
test = False

def myfunc():
    # os.chdir(beginning path if not already navigated here)
    time = []
    reelist = []
    cwd = os.getcwd()

    for i in range(1,21):
        temptime=[]
        tempree=[]
        path = str(cwd+"/result"+str(i))
        os.chdir(path)
        # create gyration.dat for that file
        gyrationTensor_2D.func(cwd)
        # access gyration.dat
        name = path+'/gyration.dat'
        file=open(name, "r")
        line=file.readline() #Gyration tensor measures and end to end distance with time 
        line=file.readline() #time	rg	rg_0	rg_1	k^2	ree
        for j in range(200):
            line=file.readline()
            t,rg,rg_0,rg_1,k2,ree = line.split()
            temptime.append(500*float(t))
            tempree.append(float(ree))
        time=temptime
        #print(tempree)
        reelist.append(tempree)
    #print(time)
    #print(reelist)
    reearr = np.asarray(reelist,dtype=float)
    print(reearr)
    # then want to average along the correct!! axis
    print(np.shape(reearr))
    reeav = np.average(reearr, axis=0)
    reestd = np.std(reearr, axis=0)
    reeerr = np.std(reearr, axis=0)/np.sqrt(np.shape(reearr)[0])
    print("should be 20: ",np.shape(reearr)[0])
    print(reeerr)
    print(reeav)

    # write out everything in case need to come back to
    outfile = open(cwd+"/endToEnd.dat", 'w')
    outfile.write("Average end to end distance with time \n")
    outfile.write("time\tree\tstd\tsterr\n")
    for t in range(len(time)):
        outfile.write("%e\t%e\t%e\t%e\n"%(time[t],reeav[t],reestd[t],reeerr[t]))	
    outfile.close()

    fig,ax = plt.subplots(1)
    #plt.plot(time, reeav)
    figdir = '/home/s1954660/Desktop/summer/analysis/endToEndVsTime'
    activity,chunks = gyrationTensor_2D.read(cwd)
    figname=str(figdir+"endToEnd_a"+str(activity)+"_c"+str(chunks)+".pdf")
    ed.errorbar_fill(np.asarray(time,dtype=float),reeav,yerr=reestd) #probs calculate an error
    title = str("End to end distance averaged over 20 runs, with standard devation. Act: "+str(activity)+". Chunks: "+str(chunks)+".")
    ax.set_title(title, fontsize=13)
    ax.set_xlabel(r"Time (simulation timesteps)",fontsize=FS)
    ax.set_ylabel(r"Average end to end distance (simulation cells)",fontsize=FS)
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    plt.savefig(figname,format='pdf', dpi = 'figure')
    plt.close(fig)
    
if test == True:
    myfunc()

else: #same but over all params (?!)
    # start in folder cluster01
    clustdir = os.getcwd()
    datedirlist = os.listdir(clustdir)
    for date in datedirlist: #go through date folders
        datedir = os.path.join(clustdir, date)
        #go through activities
        actdirlist = os.listdir(datedir)
        for act in actdirlist:
            if 'act' in act:
                actdir = os.path.join(datedir, act)
                #go through chunks
                chunkdirlist = os.listdir(actdir)
                for chunk in chunkdirlist:
                    chunkdir = os.path.join(actdir, chunk)
                    os.chdir(chunkdir)
                    myfunc()