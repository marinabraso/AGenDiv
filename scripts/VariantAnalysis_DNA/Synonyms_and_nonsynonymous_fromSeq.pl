#!/bin/python3
import sys


###################################
# Functions
###################################

def retrieve_Syn_nonSyn(seq):     
    tableSNS = {
        'ATA':'NND', 'ATC':'NND', 'ATT':'NND', 'ATG':'NND',
        'ACA':'NNS', 'ACC':'NNS', 'ACG':'NNS', 'ACT':'NNS',
        'AAC':'NND', 'AAT':'NND', 'AAA':'NND', 'AAG':'NND',
        'AGC':'NND', 'AGT':'NND', 'AGA':'DND', 'AGG':'DND',                 
        'CTA':'DNS', 'CTC':'NNS', 'CTG':'DNS', 'CTT':'NNS',
        'CCA':'NNS', 'CCC':'NNS', 'CCG':'NNS', 'CCT':'NNS',
        'CAC':'NND', 'CAT':'NND', 'CAA':'NND', 'CAG':'NND',
        'CGA':'DNS', 'CGC':'NNS', 'CGG':'DNS', 'CGT':'NNS',
        'GTA':'NNS', 'GTC':'NNS', 'GTG':'NNS', 'GTT':'NNS',
        'GCA':'NNS', 'GCC':'NNS', 'GCG':'NNS', 'GCT':'NNS',
        'GAC':'NND', 'GAT':'NND', 'GAA':'NND', 'GAG':'NND',
        'GGA':'NNS', 'GGC':'NNS', 'GGG':'NNS', 'GGT':'NNS',
        'TCA':'NNS', 'TCC':'NNS', 'TCG':'NNS', 'TCT':'NNS',
        'TTC':'NND', 'TTT':'NND', 'TTA':'DND', 'TTG':'DND',
        'TAC':'NND', 'TAT':'NND', 'TAA':'NDD', 'TAG':'NND',
        'TGC':'NND', 'TGT':'NND', 'TGA':'NDD', 'TGG':'NND',
    }
    Syn_nonSyn =""
    #if len(seq)%3 == 0:
    for i in range(0, len(seq), 3):
    	if i+3 <= len(seq):
	        codon = seq[i:i + 3]
	        Syn_nonSyn+= tableSNS[codon]
    return Syn_nonSyn

def retrieve_Complementary(seq):     
    table = { 'A':'T', 'T':'A', 'G':'C', 'C':'G' }
    complement = ""
    for i in range(0, len(seq)):
        complement+= table[seq[i]]
    return complement

###################################
# Main
###################################
 
sequence=sys.argv[1]
sense=sys.argv[2]

#if len(sequence)%3 != 0:
#	print("Error, sequence length % 3 != 0 (sequence length =", len(sequence), ")")
#	exit()

if sense == "+" :
	print(retrieve_Syn_nonSyn(sequence))
else :
	#print(sequence)
	#print(sequence[::-1])
	#print(retrieve_Complementary(sequence[::-1]))
	print(retrieve_Syn_nonSyn(retrieve_Complementary(sequence[::-1]))[::-1])






