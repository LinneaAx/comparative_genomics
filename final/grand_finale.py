#!/usr/bin/env python2
import re
from collections import defaultdict
import math 
import numpy as np
from string import maketrans
import sys

##############################################################
######################### ORF FINDER #########################
##############################################################

def complementDNA(input_genome):
    s = ''
    with open(input_genome, 'r') as f:
        for line in f:
                if not line.startswith('>'):
                    s += line 
                    for ch in f:
                        if ch=='A':
                            s=s+'T'
                        elif ch=='T':
                            s=s+'A'
                        elif ch=='G':
                            s=s+'C'
                        else:
                            s=s+'G'
                    return s #complement
                    
def rev_complementDNA(input_genome):
    s = ''
    with open(input_genome, 'r') as f:
        for line in f:
                if not line.startswith('>'):
                    s += line 
                    for ch in f:
                        if ch=='A':
                            s=s+'T'
                        elif ch=='T':
                            s=s+'A'
                        elif ch=='G':
                            s=s+'C'
                        else:
                            s=s+'G'
                    return s[::-1] #reverse complement                    
    

def ORF_finder(input_genome):
#The input should be a genome file in FASTA format; the output should also be a file in FASTA format with separate entries for each ORF gene sequences and unique names identifying these ORFs

    complement = complementDNA(input_genome)
    reverse_complement = rev_complementDNA(input_genome)
    with open(input_genome + '_ORF_finder.fasta', 'w') as w:  
        for i in range(len(complement)-2): #iterates over all possible positions where a codon begin, so all except last 2
            
            correct_way = re.compile(r'(?=(ATG(?:...)*?)(?:TAG|TAA|TGA))') #starts with ATG then go to the next stop codon...so on
            
            #return set(correct_way.findall(complement))
            #i need to exclude the strings with start codons inside
            #read backwards solve ^ problem or $
            count = 0
            for i in correct_way.findall(complement):
                if len(i) > 100: #chose 100 because of Karlin et al. reference
                    #tack Kajetan format help
                    w.write('>ORF_{}\n{}\n'.format(count, ''.join(str(i)))) #took out set so takes gene duplication into account, since findall takes into account overlaps...way too many kept the set
                    count+=1
            
                   
            count = 0
            for i in correct_way.findall(reverse_complement):
                if len(i) > 100:       
                    w.write('>ORF_rev_{}\n{}\n'.format(count, ''.join(str(i))))
                    count+=1
                    
            break
            #REGEX FOR THE WIN
         


###################################################################
######################### STATISTICS TOOL #########################
###################################################################

def compute_gc(input_genome):
    sequence = ''
    with open(input_genome + '_GC_content_Frequency', 'w') as w:
        with open(input_genome, 'r') as f:
    #with open(filename, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    sequence += line
                    C=sequence.count('C') #count C
                    G=sequence.count('G') #count G
                    total=len(sequence)
                    res = float(C+G)/total #sum of GC count/total
            w.write('The GC content frequency is '+ str(res))

#look for GC content then sum/sequence

def compute_dinucleo(input_genome):
    sequence = ''
    item = 'AG'
    temp_dict = defaultdict(int)
    with open(input_genome + '_Dinucleotide_Frequency', 'w') as w:
        with open(input_genome, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    sequence += line
                for item in range(len(sequence)-1):
                    temp_dict[sequence[item:item+2]] +=1
                for k,v in sorted(temp_dict.items()):
                    w.write('Dinucleotide frequency of '+ k+' is '+ str(v)+'\n')

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
                
                    
def compute_aa(input_fasta):
    list_aa = ['A', 'G', 'I', 'L', 'P', 'V', 'F', 'W','Y', 'D', 'E', 'R', 'H', 'K', 'S', 'T', 'C', 'M', 'N', 'Q'] 
    trans = trans_aa(input_fasta)
    with open(input_fasta + '_amino_acid_frequency', 'w') as w:
        for item in list_aa:
            i = trans.count(item)
            res = float(i)/(len(trans)-1)
            w.write('Amino acid frequency of '+item+' is ' + str(res) + '\n')

       
if __name__ == '__main__':
    input_genome = sys.argv[1]
    compute_gc(input_genome)   
    compute_dinucleo(input_genome) 
    ORF_finder(input_genome)
