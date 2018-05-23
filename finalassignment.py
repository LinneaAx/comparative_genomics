# -*- coding: utf-8 -*-
"""
Created on Wed May 23 15:51:57 2018

@author: u2365
"""

#makes the genome file into a string for counting
def makestring(input_file):
    with open(input_file, 'r') as genome:
        genome = genome.readlines()
        genome = genome[1:]
        
        data = str()
        for lines in genome:
            if lines.startswith('>') == False:
                data += lines
                return(data)


#function counting number of each nucleotide in a genome file
def countnucl(input_file):
        genome = makestring(input_file)
            
        gcontent = genome.count('G')
        ccontent = genome.count('C')
        acontent = genome.count('A')
        tcontent = genome.count('T')
        
        return(gcontent, ccontent, acontent, tcontent)

def countdinucl(input_file):
    genome = makestring(input_file)
        dilist = {'GC', 'CC', 'TC', 'AC', 'GG', 'CG', 'AG', 'TG', 'AT', 'TT', 'GT', 'CT', 'AA', 'TA', 'CA', 'GA'}
        for i in dilist:
                    
        

def seq_com(seq):
    
    complement = seq.replace('G', 'c').replace('C', 'g').replace('A','t').replace('T', 'a').upper()
    return(complement)
    
if __name__ == '__main__':
    string = makestring('03_test.txt')    
    print(string)
    result = countnucl('03_test.txt')
    print(result)