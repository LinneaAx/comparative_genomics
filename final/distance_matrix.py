#!/usr/bin/env python2
import re
from collections import defaultdict
import math 
import numpy as np
from string import maketrans
import sys

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
        print(d_matrix)
    
if __name__ == '__main__':
    input_genomes = sys.argv[1:]
    distance_matrix(input_genomes)
