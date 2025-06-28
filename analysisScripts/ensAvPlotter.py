"""
	Ensemble-average rendering script.
	To plot ensemble-averages in time.

	Created by Tyler Shendruk
"""

from pylab import *
from subprocess import call
from scipy import integrate
import os
import argparse

###########################################################
### Set up argsparse
###########################################################
parser = argparse.ArgumentParser(description='Ensemble-average rendering script.')
parser.add_argument('dataName', type=str,
					help="Path to data (any ensemble-averaged MPCD output)")
parser.add_argument('yaxis', type=str,help="y axis label")
parser.add_argument('-c','--columns', nargs='+', type=str,
					help="data --- What columns would you like to plot? 0 will be x-axis, 1 will be y-axis, etc. If you want to plot multiple columns, use -c 0 1 2 ...",
					default=['1'])
args = parser.parse_args()

###########################################################
### Format and style
###########################################################
# Use our custom style and colours
plt.style.use('shendrukGroupStyle')
import shendrukGroupFormat as ed

###########################################################
### Initialize
###########################################################
print("Arguments:")
for arg, value in vars(args).items():
	print(f"\t{arg}: {value}")
dataName = args.dataName
cols = args.columns
yaxis = args.yaxis
if not os.path.isfile(dataName):
	print("%s not found."%dataName)
	exit()
# Check that the columns are valid
numC = len(cols)
for i in range(numC):
	try:
		cols[i] = int(cols[i])
	except ValueError:
		print("Column %s is not an integer."%cols[i])
		exit()
	if cols[i] < 0:
		print("Column %s is negative."%cols[i])
		exit()
avOutName=dataName.split(".dat")[0]+"_av"

###########################################################
### Initialize/Read header
###########################################################
print( 'Reading file for header ...' )
if not os.path.isfile(dataName):
	print("%s not found."%dataName)
	exit()
infile = open(dataName,"r")
# Toss header
for i in range(13):
	line = infile.readline()
# Check that number of columns in the header matches the number of columns requested
headerCols = line.split("\t")
while '' in headerCols:
    headerCols.remove('')
lenHeader = len(headerCols)
for i in range(numC):
	if( int(cols[i])>=lenHeader ):
		print("Column %d requested, but only %d columns in header." % (int(cols[i]), lenHeader))
		infile.close()
		exit()

print( 'Reading file for data ...' )
time=[]
averages = []
for i in range(numC):
	averages.append([])
while infile:
	line = infile.readline()
	if( not line ):
		break
	else:
		line=line.split("\t")
		for i in range(numC):
			try:
				averages[i].append(float(line[cols[i]]))
			except IndexError:
				print("Column %d requested, but only %d columns in data." % (int(cols[i]), len(line)))
				infile.close()
				exit()
		time.append(float(line[0]))
time= np.array(time)
for i in range(numC):
	averages[i] = np.array(averages[i])

###########################################################
### Plot
###########################################################
for i in range(len(headerCols)):
	headerCols[i]=headerCols[i].replace("_"," ")

fig1 = plt.figure(1)
for i in range(numC):
	plot(time,averages[i],label="%s"%headerCols[cols[i]]) #,linewidth=2.0
legend(loc='best')
xlabel(r'Time, $t$')
ylabel(r'%s'%yaxis)
savefig('%s.png'%(avOutName), bbox_inches='tight', transparent=True)
savefig('%s.pdf'%(avOutName), bbox_inches='tight', transparent=True)
show()

exit()