#!/bin/python3
import sys


###################################
# Functions
###################################

def retrieve_Translated(seq):
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'.', 'TAG':'.',
        'TGC':'C', 'TGT':'C', 'TGA':'.', 'TGG':'W',
    }
    protein =""
    #if len(seq)%3 == 0:
    for i in range(0, len(seq), 3):
    	if i+3 <= len(seq):
	        codon = seq[i:i + 3]
	        protein+= table[codon]
    return protein

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
	print(retrieve_Translated(sequence))
else :
	print(retrieve_Translated(retrieve_Complementary(sequence[::-1])))






