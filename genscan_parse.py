
#for parsing out the nucleotide sequences 
def nucleotide_seq(file):
    with open(file, 'r') as o:
        write_file = open('30.both.txt_out_nucleotide', 'w')
        nucleotide_set = set('ctga\n') #the nucleotide sequences are all a subset of this
        for lines in o:
            if lines.endswith('bp\n'): #all the nucleotide sequence headers end with bp\n
                write_file.write(lines)
            if set(lines).issubset(nucleotide_set):
                if lines != '\n':
                    write_file.write(lines)

#for parsing out the protein sequences                     
def protein_seq(file):
    with open(file, 'r') as o:
        write_file = open('30.both.txt_out_proteins', 'w')
        for lines in o:
            if lines.startswith('>'): #all the protein sequence headers starts with >
                write_file.write(lines)
            if lines.strip('\n').isupper(): #all the proteins sequences are in upper case letters
                write_file.write(lines)
                
if __name__ == '__main__':
    nucleotide_seq('30.both.txt')
    protein_seq('30.both.txt')

