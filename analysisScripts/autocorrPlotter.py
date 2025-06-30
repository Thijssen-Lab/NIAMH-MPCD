"""
	NIAMH-MPCD 
	Autocorrelation rendering script
	To average autocorrelations in time

	Uses shendrukGroupStyle for formatting (from https://github.com/Shendruk-Lab/MPCDDefectLoader)
	Must install it or remove calls to shendrukGroupFormat

	Created by Tyler Shendruk
"""

from pylab import *
from subprocess import call
from scipy import integrate
import os
import argparse

# Use our custom style and colours
plt.style.use('shendrukGroupStyle')
import shendrukGroupFormat as ed

###########################################################
### Set up argsparse
###########################################################
parser = argparse.ArgumentParser(description='Autocorrelation rendering script.')
parser.add_argument('dataName', type=str,
					help="Path to data (should be an autocorrelation output)")
parser.add_argument('yaxis', type=str,
					help="y axis label --- What is this the autocorrelation of?")
parser.add_argument("start", type=int, help="Starting timestep for averaging")
parser.add_argument("finish", type=int, help="Finishing timestep for averaging")
parser.add_argument("-m", "--makeMovie", type=int, default=0,
					help="Whether or not to animate temporal data (0 or 1)")
parser.add_argument("-k", "--keepFrames", type=int,
                    help="0=don't keep (delete) frames; 1=keep frames",
                    default=0)
parser.add_argument("-a", "--plotAv", type=int,
					help="Whether or not to plot average (0 or 1)",
					default=0)
args = parser.parse_args()

###########################################################
### Read arguments
###########################################################
print("Arguments:")
for arg, value in vars(args).items():
	print(f"\t{arg}: {value}")
dataName = args.dataName
yaxis = args.yaxis
start = args.start
finish = args.finish
makeMovie = args.makeMovie
keepFrames = args.keepFrames
plotAv = args.plotAv

###########################################################
### Assumed arguments that a user could change
###########################################################
avOutName=dataName.split(".dat")[0]+"_av"

###########################################################
### Format and style
###########################################################
#Animation stuff
bitrate=5000
framerate=12		#Number of frames per second in the output video
	# Note that the ideal for this is _very_ variable so play around with it
codec='libx264'		#Other options include mpeg4
suffix='.mp4'

###########################################################
### Number of bins
###########################################################
print( 'Reading file for bin size ...' )
numBins=101
if not os.path.isfile(dataName):
	print("%s not found."%dataName)
	exit()
infile = open(dataName,"r")
#toss header
for i in range(13):
	line = infile.readline()
# Read the time
line=infile.readline()
while line!="\n":
	# Read the first line of the data
	oldLine = line
	line = infile.readline()
oldLine = oldLine.split("\t")
numBins=int(oldLine[1])
print( "\tnumBins: %s"%numBins )

###########################################################
### Animation/Average
###########################################################
# Figure
fig1 = plt.figure(1)
print( 'Reading file for data ...' )
if not os.path.isfile(dataName):
	print("%s not found."%dataName)
	exit()
infile = open(dataName,"r")
#toss header
for i in range(13):
	line = infile.readline()
t=0
rad=linspace(0,numBins,num=numBins+1,endpoint=True,dtype=float)
corr=zeros(shape=(numBins+1),dtype=float)
corrAv=zeros(shape=(numBins+1),dtype=float)
corrStd=zeros(shape=(numBins+1),dtype=float)
minC=999999.
maxC=-999999.
while infile:
	line = infile.readline()
	if( not line ):
		break
	else:
		time=float(line) 			# Read the time
		for i in range(numBins+1):	# Read the correlation function data
			line = infile.readline()
			toss,x,c = line.split("\t")
			corr[i]=float(c)
		line = infile.readline()	# Read the empty line
		minC = minimum(minC, corr.min())
		maxC = maximum(maxC, corr.max())
		if(t>=start and t<=finish):
			for i in range(numBins+1):
				corrAv[i]+=corr[i]
				corrStd[i]+=corr[i]**2
		if(makeMovie):
			fig1 = plt.figure(1)
			plt.clf()
			title(r'Time, %s$\tau$'%(time))
			plot(rad,corr,'-o')
			xlabel(r'Radial distance, $r$')
			ylabel(r'Autocorr. %s'%yaxis)
			plt.axis(ymax=maxC, ymin=minC)
			name='frame%04d.png'%(t)
			savefig(name, bbox_inches='tight')
		t+=1
infile.close()
if(makeMovie):
	#Animate
	print( "\tAnimating ..." )
	name='%s_animation%s'%(avOutName,suffix)
	myCommand="rm %s"%name
	call(myCommand,shell=True)
	myCommand = "ffmpeg -f image2 -r %d"%(framerate)+" -i frame%04d.png"+" -vcodec %s -b %dk -r %d %s"%(codec,bitrate,framerate,name)
	call(myCommand,shell=True)
	if not keepFrames:
		myCommand="rm frame*.png"
		call(myCommand,shell=True)

###########################################################
### Average
###########################################################
if t<finish:
	finish=t
T=float(finish-start+1)
for i in range(numBins+1):
	if T>1:
		corrStd[i]=sqrt(corrStd[i]/T - (corrAv[i]/T)**2)
	corrAv[i]/=T
corrErr=corrStd/sqrt(T)

infile = open(dataName,"r")
# Keep header
header=[]
for i in range(12):
	header.append(infile.readline())
infile.close()
outfile=open(avOutName+".dat", "w")
for h in header:
	outfile.write(h)
outfile.write("dt\tav\tStd\tErr\n")
for i in range(numBins):
	outfile.write("%e\t%e\t%e\t%e\n"%(rad[i],corrAv[i],corrStd[i],corrErr[i]))
outfile.close()

if(plotAv):
	fig1 = plt.figure(1)
	plt.clf()
	ed.errorbar_fill(rad,corrAv,corrErr)
	xlabel(r'Radial distance, $r$')
	ylabel(r'Autocorr. %s'%yaxis)
	savefig('%s.png'%(avOutName), bbox_inches='tight', transparent=True)
	savefig('%s.pdf'%(avOutName), bbox_inches='tight', transparent=True)
	show()

exit()