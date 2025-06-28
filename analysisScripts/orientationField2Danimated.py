"""
	Orientation field (director) rendering script.
	Uses defect handler (from https://github.com/Shendruk-Lab/MPCDDefectLoader)

	Originally from Tyler N. Shendruk
	Modified by Timofey Kozhukhov
"""

from pylab import *
from subprocess import call
import sys
import os
import json
import argparse

from defectHandler import getDefectData

###########################################################
### Set up argsparse
###########################################################
parser = argparse.ArgumentParser(description='Orientation field rendering '
											 'script.')
parser.add_argument("dataname", type=str, help="Path to the data (should be "
											   "directorfield.dat)")
parser.add_argument("inputname", type=str, help="Path to input .json file")
parser.add_argument("start", type=int, help="Starting timestep for averaging")
parser.add_argument("finish", type=int, help="Finishing timestep for averaging")
parser.add_argument("--qx", type=int, help="Only show every qx arrow in x", default=1)
parser.add_argument("--qy", type=int, help="Only show every qy arrow in y",default=1)
parser.add_argument("avdim", type=str, help="Dimension to 'slice' over")
parser.add_argument("-c", "--length", type=float, help="Length of director lines", default=0.5)
parser.add_argument("-a", "--myAspect", type=str, help="'auto' or 'equal'", default="auto")
parser.add_argument("-k", "--keepFrames", type=int, help="0=don't keep (delete) frames; 1=keep frames", default=0)
parser.add_argument("-d", "--defectData", type=str, help="Path to defect data (if any)", default="")
args = parser.parse_args()

###########################################################
### Read arguments
###########################################################
print("Arguments:")
for arg, value in vars(args).items():
	print(f"\t{arg}: {value}")
dataName = args.dataname
inputName = args.inputname
start = args.start
finish = args.finish
qx = args.qx
qy = args.qy
avdim = args.avdim
c = args.length
myAspect = args.myAspect
keepFrames = args.keepFrames
defectData = args.defectData

###########################################################
### Format and style
###########################################################
# Use our custom style and colours
plt.style.use('shendrukGroupStyle')
import shendrukGroupFormat as ed
# Colour map to use
myMap=ed.plasma
# Adjust line width
myLW=1.0
# adjust line length
c *= (qx+qy)/2
#Animation stuff
bitrate=5000
framerate=12		#Number of frames per second in the output video
codec='libx264'		#Other options include mpeg4
suffix='.mp4'

###########################################################
### Initialize
###########################################################
if avdim=='x':
	dim=0
	d1=1
	d2=2
elif avdim=='y':
	dim=1
	d1=0
	d2=2
elif avdim=='z':
	dim=2
	d1=0
	d2=1
else:
	print( "avdim must be 'x', 'y' or 'z' - not %s"%avdim )
	exit()

###########################################################
### Read json
###########################################################
if not os.path.isfile(inputName):
	print("%s not found."%inputName)
	exit()
with open(inputName, 'r') as f:
  input = json.load(f)
xyzSize=array([30,30,1])
dim=2
if "domain" in input:
	xyzSize[0]=input['domain'][0]
	if(len(input['domain'])>1):
		xyzSize[1]=input['domain'][1]
	else:
		xyzSize[1]=1
	if(len(input['domain'])>2):
		xyzSize[2]=input['domain'][2]
	else:
		xyzSize[2]=1
# Data
XYZ = zeros(shape=(3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
DIR = zeros(shape=(3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
MEAN = zeros(shape=(3,xyzSize[d1],xyzSize[d2]),dtype=float)
MAG = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)
XY = zeros(shape=(2,xyzSize[d1],xyzSize[d2]),dtype=float)
S = zeros(shape=(xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
AVS = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)

### Setup the animation
# Figure
fig1 = plt.figure(1)
if myAspect == 'auto':
    shrink_factor = 1.0
else:
    shrink_factor = float(xyzSize[d2])/float(xyzSize[d1])
#Create the colorbar
CS3 = imshow(AVS.T,cmap=myMap,vmin=0, vmax=1,aspect=myAspect,extent=[0,xyzSize[d1],0,xyzSize[d2]])
cb=colorbar(CS3,shrink=shrink_factor,aspect=20*shrink_factor,pad=0.04)
cb.ax.set_ylabel(r'Scalar order, $S$')
# Make labels
if avdim=='x':
	labX='y'
	labY='z'
elif avdim=='y':
	labX='x'
	labY='z'
elif avdim=='z':
	labX='x'
	labY='y'

###########################################################
### defect handling if needed
###########################################################
LOADDEFECTS = False
defects = []
if defectData != "":
	print("Loading defects for rendering ...")
	LOADDEFECTS = True

	defContainer = getDefectData(defectData, np.array([xyzSize[0], xyzSize[1], xyzSize[2]]))
	for defList in defContainer:
		defects.append(defList.defectList)
	print("Finished loading defects")

###########################################################
### Read the data for animation
###########################################################
print( 'Reading data ...' )
if not os.path.isfile(dataName):
	print("%s not found."%dataName)
	exit()
infile = open(dataName,"r")
# Toss header
for i in range(13):
	line = infile.readline()
i=0
j=0
n=-1
while infile:
	i=i+1
	line = infile.readline()
	if( not line ):
		break
	else:
		if j>=start:
			t,Qx,Qy,Qz,Vx,Vy,Vz,s = line.split("\t",8)
			XYZ[0][int(Qx)][int(Qy)][int(Qz)] = float(Qx) + 0.5
			XYZ[1][int(Qx)][int(Qy)][int(Qz)] = float(Qy) + 0.5
			XYZ[2][int(Qx)][int(Qy)][int(Qz)] = float(Qz) + 0.5
			DIR[0][int(Qx)][int(Qy)][int(Qz)] = float(Vx)
			DIR[1][int(Qx)][int(Qy)][int(Qz)] = float(Vy)
			DIR[2][int(Qx)][int(Qy)][int(Qz)] = float(Vz)
			S[int(Qx)][int(Qy)][int(Qz)] = float(s)
	if i==xyzSize[0]*xyzSize[1]*xyzSize[2]:
		j=j+1
		if j>finish:
			break
		if j<start or j>finish:
			pass
		else:
			#Sum
			if avdim=='x':
				for x in range(xyzSize[0]):
					for y in range(xyzSize[1]):
						for z in range(xyzSize[2]):
							MEAN[0][y][z]=MEAN[0][y][z]+DIR[0][x][y][z]
							MEAN[1][y][z]=MEAN[1][y][z]+DIR[1][x][y][z]
							MEAN[2][y][z]=MEAN[2][y][z]+DIR[2][x][y][z]
							AVS[y][z]=AVS[y][z]+S[x][y][z]
				for y in range(xyzSize[1]):
					for z in range(xyzSize[2]):
						for i in range(3):
							MEAN[i][y][z]/=xyzSize[0]
						AVS[y][z]/=xyzSize[0]
			elif avdim=='y':
				for x in range(xyzSize[0]):
					for y in range(xyzSize[1]):
						for z in range(xyzSize[2]):
							MEAN[0][x][z]=MEAN[0][x][z]+DIR[0][x][y][z]
							MEAN[1][x][z]=MEAN[1][x][z]+DIR[1][x][y][z]
							MEAN[2][x][z]=MEAN[2][x][z]+DIR[2][x][y][z]
						AVS[x][z]=AVS[x][z]+S[x][y][z]
				for x in range(xyzSize[0]):
					for z in range(xyzSize[2]):
						for i in range(3):
							MEAN[i][x][z]/=xyzSize[1]
						AVS[x][z]/=xyzSize[1]
			elif avdim=='z':
				for x in range(xyzSize[0]):
					for y in range(xyzSize[1]):
						for z in range(xyzSize[2]):
							MEAN[0][x][y]=MEAN[0][x][y]+DIR[0][x][y][z]
							MEAN[1][x][y]=MEAN[1][x][y]+DIR[1][x][y][z]
							MEAN[2][x][y]=MEAN[2][x][y]+DIR[2][x][y][z]
							AVS[x][y]=AVS[x][y]+S[x][y][z]
				for x in range(xyzSize[0]):
					for y in range(xyzSize[1]):
						for i in range(3):
							MEAN[i][x][y]/=xyzSize[2]
						AVS[x][y]/=xyzSize[2]
			for x in range(xyzSize[d1]):
				for y in range(xyzSize[d2]):
					MAG[x][y]=sqrt( MEAN[0][x][y]**2+MEAN[1][x][y]**2+MEAN[2][x][y]**2 )

			#Save the instantaneous or current velocity field frame
			# Make Mesh
			if avdim=='x':
				for y in range(xyzSize[1]):
					for z in range(xyzSize[2]):
						XY[0][y][z]=XYZ[d1][0][y][z]
						XY[1][y][z]=XYZ[d2][0][y][z]
			elif avdim=='y':
				for x in range(xyzSize[0]):
					for z in range(xyzSize[2]):
						XY[0][x][z]=XYZ[d1][x][0][z]
						XY[1][x][z]=XYZ[d2][x][0][z]
			elif avdim=='z':
				for x in range(xyzSize[0]):
					for y in range(xyzSize[1]):
						XY[0][x][y]=XYZ[d1][x][y][0]
						XY[1][x][y]=XYZ[d2][x][y][0]
			# Save frame
			n=n+1
			fig1 = plt.figure(1)
			plt.cla()
			quiver( XY[0][::qx, ::qy], XY[1][::qx, ::qy], 
					c*MEAN[0][::qx, ::qy], c*MEAN[1][::qx, ::qy], 
					AVS[::qx, ::qy], cmap=myMap, clim=(0, 1), scale=50/c,
					width=0.005*myLW, headlength=0, headwidth=0, 
					headaxislength=0, pivot='middle')
			# load defects and draw them as necessary
			# FIXME: only works for 2d for now, doesnt take into account d1 or d2
			if LOADDEFECTS and (j < len(defects)):
				for defect in defects[j-1]: # j is not 0 indexed reeeeee
					defect.drawDefect(c*2, myLW*2)
			xlabel(r'$%s$'%labX)
			ylabel(r'$%s$'%labY)
			plt.axis(xmax=xyzSize[d1], xmin=0, ymax=xyzSize[d2], ymin=0)
			name='frame%04d.png'%(n)
			# uncomment below for snapshots
			plt.axis('off') 
			plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
			savefig( name, bbox_inches='tight', pad_inches=0 )
		#Zero matrix
		DIR= zeros( (3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float )
		MEAN = zeros(shape=(3,xyzSize[d1],xyzSize[d2]),dtype=float)
		AVS = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)
		i=0
infile.close()

###########################################################
### Visualize data
###########################################################
print( "Animating ..." )
name='orientation_%s%s'%(avdim,suffix)
myCommand="rm %s"%name
call(myCommand,shell=True)
myCommand = "ffmpeg -f image2 -r %d"%(framerate)+" -i frame%04d.png"+" -vcodec %s -b %dk -r %d %s"%(codec,bitrate,framerate,name)
call(myCommand,shell=True)
if not keepFrames:
	myCommand="rm frame*.png"
	call(myCommand,shell=True)
