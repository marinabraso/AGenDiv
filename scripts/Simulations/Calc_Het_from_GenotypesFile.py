#!/usr/bin/python3

################
## Debug & develop

import sys
import os
import subprocess
import re
import gzip



GenotypesFile = sys.argv[1]
Output = sys.argv[2]
chrlen = int(sys.argv[3])

###############
## General parameters
os.system("mkdir -p $(dirname "+ Output + ")")
os.system("rm " + Output + " 2> ~/null")
Output = Output[:-3] # remove .gz

###############
## Functions
def printHet(ginfo, clen, G, ofile):
	Hetvec=[]
	for i in range(0, len(ginfo)-1, 2):
		sumHet=0
		arr1=[*ginfo[i]]
		arr2=[*ginfo[i+1]]
		for n in range(len(arr1)):
			if(arr1[n] != arr2[n]):
				sumHet+=1
		Hetvec.append(sumHet/clen)
	print(Hetvec)
	ofile.write( G + "\t" + '\t'.join(str(x) for x in Hetvec) + "\n")

###############
## Read Mutations
inputfile = gzip.open(GenotypesFile,'rt')
outputfile = open(Output, 'w')
genotinfo = []
fst=0
for l in inputfile:
	if(re.match("Generation ", l)):
		if fst == 1:
			printHet(genotinfo, chrlen, Gen, outputfile)
			genotinfo = []
		Gen=l.rstrip("\n").split(" ")[1]
	else:
		fst=1
		genotinfo.append(l)
printHet(genotinfo, chrlen, Gen, outputfile)
inputfile.close()
outputfile.close()

os.system("gzip " + Output )



