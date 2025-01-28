#!/usr/bin/python3

################
## Debug & develop

import sys
import os
import subprocess
import re
import gzip



GenotypesFile = sys.argv[1]
OutputMF = sys.argv[2]
OutputNA = sys.argv[3]

###############
## General parameters
os.system("mkdir -p $(dirname "+ OutputMF + ")")
os.system("rm " + OutputMF + " 2> ~/null")
OutputMF = OutputMF[:-3] # remove .gz
print(OutputMF)
os.system("rm " + OutputNA + " 2> ~/null")
OutputNA = OutputNA[:-3] # remove .gz

###############
## Functions
def print_MajorFreq(ginfo, ofile):
	majfreq = []
	totalnumalleles = len(ginfo[0].split("\t")[4:len(ginfo[0].split("\t"))])
	print(totalnumalleles)
	for i in range(0, len(ginfo), 1):
		alleleslist = ginfo[i].split("\t")[4:len(ginfo[i].split("\t"))]
		max = 0
		for a in set(alleleslist):
			freq = alleleslist.count(a)
			if max < freq:
				max = freq
		if max < totalnumalleles:
			majfreq.append(max)
	sfs = []
	for i in range(totalnumalleles-1, int(totalnumalleles/2)-1, -1):
		sfs.append(str(majfreq.count(i)))
	ofile.write(ginfo[0].split("\t")[0] + "\t" + ginfo[0].split("\t")[1] + "\t" + ginfo[0].split("\t")[2] + "\t" + "\t".join(sfs) + "\n")
	#print(ginfo[0].split("\t")[0] + "\t" + ginfo[0].split("\t")[1] + "\t" + ginfo[0].split("\t")[2] + "\t" + "\t".join(sfs))

def print_NumAlleles(ginfo, ofile):
	numalleles = []
	for i in range(0, len(ginfo), 1):
		allelesset = set(ginfo[i].split("\t")[4:len(ginfo[i].split("\t"))])
		if len(allelesset) > 1:
			numalleles.append(len(allelesset))
	freqnumalleles = []
	numsitesaccounted = 0
	for i in range(2, 10, 1):
		numsitesaccounted += numalleles.count(i)
		freqnumalleles.append(str(numalleles.count(i)))
	freqnumalleles.append(str(len(numalleles)-numsitesaccounted))
	ofile.write(ginfo[0].split("\t")[0] + "\t" + ginfo[0].split("\t")[1] + "\t" + ginfo[0].split("\t")[2] + "\t" + "\t".join(freqnumalleles) + "\n")
	#print(ginfo[0].split("\t")[0] + "\t" + ginfo[0].split("\t")[1] + "\t" + ginfo[0].split("\t")[2] + "\t" + "\t".join(freqnumalleles))



###############
## Read Mutations
inputfile = gzip.open(GenotypesFile,'rt')
outputfileMF = open(OutputMF, 'w')
outputfileNA = open(OutputNA, 'w')
genotinfo = []
fst=0
reg=""
for l in inputfile:
	if fst == 0:
		reg = l.split("\t")[0] + "_" + l.split("\t")[1] + "_" + l.split("\t")[2]
		fst = 1
	if reg == l.split("\t")[0] + "_" + l.split("\t")[1] + "_" + l.split("\t")[2]:
		genotinfo.append(l.replace(":", "\t").replace("\n", ""))
	else:
		print_MajorFreq(genotinfo, outputfileMF)
		print_NumAlleles(genotinfo, outputfileNA)
		genotinfo = []
		genotinfo.append(l.replace(":", "\t").replace("\n", ""))
		reg = l.split("\t")[0] + "_" + l.split("\t")[1] + "_" + l.split("\t")[2]

inputfile.close()
print_MajorFreq(genotinfo, outputfileMF)
outputfileMF.close()
print_NumAlleles(genotinfo, outputfileNA)
outputfileNA.close()

os.system("gzip " + OutputMF )
os.system("gzip " + OutputNA )




