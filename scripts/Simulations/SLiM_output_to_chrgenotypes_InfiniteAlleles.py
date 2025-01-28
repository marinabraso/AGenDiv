#!/usr/bin/python3

################
## Debug & develop

import sys
import os
import subprocess
import re
import gzip


SLiMinput = sys.argv[1]
Output = sys.argv[2]

###############
## General parameters
os.system("mkdir -p $(dirname " + Output + ")")
os.system("rm " + Output + "* 2> ~/null" )
Output = Output[:-3] # remove .gz

charAlleles = [ "1","2","3","4","5","6","7","8","9","a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z" ]

###############
## Functions
def printGenotypes(dmut, dchr, ofile, chars):
	variantsIDinpos = {}
	numvariantsinpos = {}
	listmutpos = []
	for mut in dmut:
		if dmut[mut] in numvariantsinpos:
			numvariantsinpos[dmut[mut]]+=1
		else:
			numvariantsinpos[dmut[mut]]=1
			listmutpos.append(dmut[mut])
		if numvariantsinpos[dmut[mut]] > len(chars):
			print("Error: not enough character in charAlleles")
			exit(1)
		variantsIDinpos[mut]=chars[numvariantsinpos[dmut[mut]]]
	print(len(listmutpos))
	for c in dchr:
		pstr=""
		mutingenome={}
		for mut in dchr[c]:
			mutingenome[dmut[mut]]=mut
		for i in listmutpos:
			if i in mutingenome:
				pstr=pstr + str(variantsIDinpos[mutingenome[i]])
			else:
				pstr=pstr + "0"
		ofile.write(pstr + "\n")

###############
## Read Mutations

bash_command = "for i in $(ls " + SLiMinput + "); do zcat < $i; done"
result = subprocess.run(bash_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
inputfile = result.stdout.split("\n")
outputfile = open(Output, 'w')
mutinfo = {}
chrinfo = {}
gmode = 0
mmode = 0
for l in inputfile:
	l=l.strip()
	if re.match("# GEN ", l):
		if(gmode==1):
			printGenotypes(mutinfo, chrinfo, outputfile, charAlleles)
		gmode = 0
		mmode = 0
		mutinfo = {}
		chrinfo = {}
		Gen=l.rstrip("\n").split(" ")[2]
		outputfile.write("Generation " + Gen + "\n")
	elif re.match("Mutations", l):
		mmode = 1
	elif re.match("Genomes", l):
		gmode = 1
		mmode = 0
	elif mmode==1:
		mutinfo[int(l.split(" ")[0])] = int(l.split(" ")[3])
	elif gmode==1:
		listmut=[]
		if len(l.split(" ")) == 3:
			listmut.append(int(l.split(" ")[2]))
		elif len(l.split(" ")) > 3:
			listmut=list(map(int, l.split(" ")[2:]))
		if l.split(" ")[0] != "":
			chrinfo[l.split(" ")[0]] = listmut
printGenotypes(mutinfo, chrinfo, outputfile, charAlleles)
outputfile.close()

os.system("gzip " + Output )

