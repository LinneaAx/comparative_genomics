#!/usr/bin/env python2
import re
from collections import defaultdict
import math 
import numpy as np
from string import maketrans
import sys


def trans_aa(input_fasta): #tanslate ORF list and then comput_aa
    sequence = ''
    trans_dict = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
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
    'TAC':'Y', 'TAT':'Y', 'TGG':'W', 'TGT':'C',
    'TGC':'C'} #didn'include stop since not in ORF

    with open(input_fasta, 'r') as f:
        temp_list = []
        for line in f:
            if not line.startswith('>'):
                sequence = line 
                new = ''
                for i in range(0, len(sequence), 3):
                    if sequence[i:i+3] in trans_dict.keys():
                        new += trans_dict[sequence[i:i+3]] #if have to create own dictionary              
                temp_list.append(new)
        return ''.join(temp_list)
                
def compute_diaa(input_fasta): 
    
    trans = trans_aa(input_fasta)
    temp_dict = defaultdict(int)
    with open(input_fasta + 'di_aa', 'w') as w:
        for line in trans:   
            for item in range(len(trans)):
                temp_dict[trans[item:item+2]] +=1
                for k,v in sorted(temp_dict.items()):
                    total = sum(temp_dict.values())
                    result = float(v)/total
                    w.write(k+' : '+ str(result)+'\n')
                    break
                    
if __name__ == '__main__':
    input_fasta = sys.argv[1]
    compute_diaa(input_fasta)
    
