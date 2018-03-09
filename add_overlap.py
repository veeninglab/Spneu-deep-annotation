"""
Input: reference fasta file
Output: reference fasta file, with the first 1000 nucleotides copied at the end.
Command: python add_overlap.py [fasta_file] [output_folder]
Example: python add_overlap.py D39V.fna D39V_annotation
"""

import sys

fastafile = sys.argv[1]
out_folder = sys.argv[2]
fin = open(fastafile,'r')
lines = fin.readlines()
fin.close()

fout = open(out_folder+'/temp_files/genome_overlap.fna','w')
fout.write(lines[0])

sequence = ''
linecount = 1

for line in lines[1:]:
	fout.write(line)
	if linecount < 500:
		seq = line.split()[0]
		sequence += seq
	linecount += 1

lastline = sequence[:1000]
fout.write(lastline)
fout.close()