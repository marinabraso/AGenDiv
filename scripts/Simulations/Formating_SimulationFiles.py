#!/usr/bin/python3

################
## Debug & develop

import sys
import os
import subprocess
import re
import gzip



GenotypesFile = sys.argv[1]
HetFile = sys.argv[2]
GenoFile = sys.argv[3]
BEDFile = sys.argv[4]
chrlen = int(sys.argv[5])

###############
## General parameters
os.system("mkdir -p $(dirname "+ HetFile + ")")
os.system("rm " + HetFile + " " + GenoFile + " " + BEDFile + " 2> ~/null")
HetFile = HetFile[:-3] # remove .gz
GenoFile = GenoFile[:-3] # remove .gz
BEDFile = BEDFile[:-3] # remove .gz

###############
## Functions
def print_genoinfo(ginfo, clen, ofile):
	for s in range(len(ginfo[0])):
		line="chr1\t1\t" + str(clen-1) + "\t" + str(s)
		for i in range(0,len(ginfo)-1,2):
			line = line + "\t" + str(ginfo[i][s]) + ":" + str(ginfo[i+1][s])
		ofile.write(line + "\n")

def print_hetinfo(ginfo, clen, ofile):
	for s in range(len(ginfo[0])):
		line="chr1\t1\t" + str(clen-1) + "\t" + str(s)
		for i in range(0,len(ginfo)-1,2):
			if ginfo[i][s] == str(ginfo[i+1][s]):
				line = line + "\t0"
			else:
				line = line + "\t1"
		ofile.write(line + "\n")

###############
## MAIN

# Dummy bed file
obedfile = open(BEDFile, 'w')
obedfile.write("chr1\t0\t" + str(chrlen) + "\n")
obedfile.close()
os.system("gzip " + BEDFile )

# Read mutations & print the other two files
inputfile = gzip.open(GenotypesFile,'rt')
ohetfile = open(HetFile, 'w')
ogenofile = open(GenoFile, 'w')
genotinfo = []
for l in inputfile:
	l=l.strip()
	if(re.match("Generation ", l)):
		genotinfo = []
	else:
		genotinfo.append(l)
print_hetinfo(genotinfo, chrlen, ohetfile)
print_genoinfo(genotinfo, chrlen, ogenofile)
inputfile.close()
ohetfile.close()
ogenofile.close()
os.system("gzip " + HetFile )
os.system("gzip " + GenoFile )



