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
def printPi(ginfo, clen, G, ofile):
	comp=0
	sumdif=0
	for i in range(len(ginfo)):
		for j in range(i+1, len(ginfo)):
			dif=0
			comp+=1
			iarr=[*ginfo[i]]
			jarr=[*ginfo[j]]
			for n in range(len(iarr)):
				if(iarr[n] != jarr[n]):
					dif+=1
			sumdif+=dif
	pi=sumdif/comp/clen
	print(pi)
	ofile.write( G + "\t" + "{:.20f}".format(pi) + "\n")

###############
## Read Mutations
inputfile = gzip.open(GenotypesFile,'rt')
outputfile = open(Output, 'w')
genotinfo = []
fst=0
for l in inputfile:
	if(re.match("Generation ", l)):
		if fst == 1:
			printPi(genotinfo, chrlen, Gen, outputfile)
			genotinfo = []
		Gen=l.rstrip("\n").split(" ")[1]
	else:
		fst=1
		genotinfo.append(l)
printPi(genotinfo, chrlen, Gen, outputfile)
inputfile.close()
outputfile.close()

os.system("gzip " + Output )



