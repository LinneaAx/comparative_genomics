#!/usr/bin/env python2
import re
from collections import defaultdict
import math 
import numpy as np
from string import maketrans
import sys
#gc, dinucleotide, amino acid

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
    #stop_codon = ['TAG','TAA','TGA']
    with open(input_genome + '_ORF_finder.fasta', 'w') as w:  
        #codon_list = [] #forward strand list
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
                    w.write('>ORF_rev_{}\n{}\n'.format(count, ''.join(str(i)))) #write function sooooooooooooo slow don't know why, change to print to see results
                    count+=1
                    
            break
            #REGEX FOR THE WIN
         


###################################################################
######################### STATISTICS TOOL #########################
###################################################################

def compute_gc(input_genome):
    sequence = ''
    #ignore = re.compile('^[ \\t]*#.*', re.IGNORECASE)
    with open(input_genome + '_GC_content_Frequency', 'w') as w:
        with open(input_genome, 'r') as f:
    #with open(filename, 'r') as f:
            #f = f.readlines()
            for line in f:
                #print(line)
                if not line.startswith('>'):
                    sequence += line
                    #print(sequence)
                    C=sequence.count('C') #count C
                    G=sequence.count('G') #count G
                    total=len(sequence)
                    #print(total)
                    res = float(C+G)/total #sum of GC count/total
                 #print(res)
            w.write('The GC content frequency is '+ str(res))
            #return res
#look for GC content then sum/sequence

def compute_dinucleo(input_genome):
    sequence = ''
    #dinucleo = ['AG', 'AA', 'AC', 'AT','CG', 'CA', 'CC', 'CT','GG', 'GA', 'GC', 'GT', 'TG', 'TA', 'TC', 'TT']
    #count = 0
    item = 'AG'
    temp_dict = defaultdict(int)
    with open(input_genome + '_Dinucleotide_Frequency', 'w') as w:
        with open(input_genome, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    sequence += line
                    #total = len(sequence)-1
                for item in range(len(sequence)-1):
                    temp_dict[sequence[item:item+2]] +=1
                for k,v in sorted(temp_dict.items()):
                    w.write('Dinucleotide frequency of '+ k+' is '+ str(v)+'\n')
                    
                '''for item in dinucleo: #specify item
                    if item in sequence:
                        count += 1
                        total=len(sequence)-1
                        res=(count)/(total)'''
                        
    #print ('The dinucleotide frequency of AG is ', res, '\n') #do one at a time or print all
# AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT  
       
def compute_aa(input_fasta): #tanslate ORF list and then comput_aa
    sequence = ''
    list_aa = ['A', 'G', 'I', 'L', 'P', 'V', 'F', 'W','Y', 'D', 'E', 'R', 'H', 'K', 'S', 'T', 'C', 'M', 'N', 'Q'] #20aa
    with open(input_fasta + '_Amino_acid_Frequency', 'w') as w:
        with open(input_fasta, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    sequence += line       
            for item in list_aa:
                i = sequence.count(item)
                res = i/(len(sequence) - 1)
                w.write('Amino acid frequency of '+item+' is ' + str(res) + '\n')
    


######################################################################## 
######################### DISTANCE MATRIX TOOL #########################
########################################################################    

def distance_matrix(input_genomes):
#The tool should compute the distance between two genomes from the DNA statistic above. The distance matrix is then used to create a species tree

#matrix one against all
# numpy array of 2 vectors
    sequence=''  
    with open('Distance Matrix', 'w') as w:                  
        with open (input_genomes[0]) as f1:
            with open (input_genomes[1]) as f2:
                with open (input_genomes[2]) as f3:
                    with open (input_genomes[3]) as f4:
                        with open (input_genomes[4]) as f5: 
                            for line in f1:
                                if not line.startswith('>'):
                                    sequence += line
                                    CG=sequence.count('C')+sequence.count('G')
                                    total=len(sequence)
                                    res1 = float(CG)/total
                            for line in f2:
                                if not line.startswith('>'):
                                    sequence += line
                                    CG=sequence.count('C')+sequence.count('G')
                                    total=len(sequence)
                                    res2 = float(CG)/total
                            for line in f3:
                                if not line.startswith('>'):
                                    sequence += line
                                    CG=sequence.count('C')+sequence.count('G')
                                    total=len(sequence)
                                    res3 = float(CG)/total
                            for line in f4:
                                if not line.startswith('>'):
                                    sequence += line
                                    CG=sequence.count('C')+sequence.count('G')
                                    total=len(sequence)
                                    res4 = float(CG)/total
                            for line in f5:
                                if not line.startswith('>'):
                                    sequence += line
                                    CG=sequence.count('C')+sequence.count('G')
                                    total=len(sequence)
                                    res5 = float(CG)/total

        d_matrix = np.zeros((5,5))
        input_genomes = [res1, res2, res3, res4, res5]
        for x in range(0,len(input_genomes)):
            for y in range(0,len(input_genomes)):
                d_matrix[x,y] = math.sqrt((input_genomes[x]-input_genomes[y])**2)
        w.write(d_matrix)
    
if __name__ == '__main__':
    input_genome = sys.argv[1]
    compute_gc(input_genome)   
    compute_dinucleo(input_genome) 
    ORF_finder(input_genome)
    #compute_aa()
    input_genomes = sys.argv[1:]
    distance_matrix(input_genomes)
