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


###############
## Functions
def printGenotypes(pmut, nmut, dchr, ofile, aseq):
	variantsIDinpos = {}
	numvariantsinpos = {}
	listmutpos = []
	for mut in pmut:
		if pmut[mut] in numvariantsinpos:
			numvariantsinpos[pmut[mut]]+=1
		else:
			numvariantsinpos[pmut[mut]]=1
			listmutpos.append(pmut[mut])
		variantsIDinpos[mut]=nmut[mut]
		if aseq[pmut[mut]] == nmut[mut]:
			print("Back mutation at position", pmut[mut], "! ", aseq[pmut[mut]], " -> ", aseq[pmut[mut]])
	for c in dchr:
		print(c)
		pstr=""
		mutingenome={}
		for mut in dchr[c]:
			mutingenome[pmut[mut]]=mut
		for i in listmutpos:
			if i in mutingenome:
				pstr=pstr + str(variantsIDinpos[mutingenome[i]])
			else:
				pstr=pstr + aseq[i]
		print("buu" + pstr)
		ofile.write(pstr + "\n")

###############
## Read Mutations
bash_command = "for i in $(ls " + SLiMinput + "); do zcat < $i; done"
result = subprocess.run(bash_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
inputfile = result.stdout.split("\n")
outputfile = open(Output, 'w')
mutposinfo = {}
mutnucinfo = {}
chrinfo = {}
ancestralseq = ""
amode = 0
gmode = 0
mmode = 0
for l in inputfile:
	l=l.strip()
	if re.match("# Ancestral sequence", l):
		amode = 1
	elif re.match("# GEN ", l):
		if(gmode==1):
			printGenotypes(mutposinfo, mutnucinfo, chrinfo, outputfile, ancestralseq)
		amode = 0
		gmode = 0
		mmode = 0
		mutposinfo = {}
		mutnucinfo = {}
		chrinfo = {}
		Gen=l.rstrip("\n").split(" ")[2]
		outputfile.write("Generation " + Gen + "\n")
		print("########## Generation " + Gen)
	elif re.match("Mutations", l):
		mmode = 1
	elif re.match("Genomes", l):
		gmode = 1
		mmode = 0
	elif amode==1:
		ancestralseq = ancestralseq + l
	elif mmode==1:
		mutposinfo[int(l.split(" ")[0])] = int(l.split(" ")[3])
		mutnucinfo[int(l.split(" ")[0])] = l.split(" ")[9]
	elif gmode==1:
		listmut=[]
		if len(l.split(" ")) == 3:
			listmut.append(int(l.split(" ")[2]))
		elif len(l.split(" ")) > 3:
			listmut=list(map(int, l.split(" ")[2:]))
		if l.split(" ")[0] != "":
			chrinfo[l.split(" ")[0]] = listmut
print(chrinfo)
printGenotypes(mutposinfo, mutnucinfo, chrinfo, outputfile, ancestralseq)
outputfile.close()

os.system("gzip " + Output )

