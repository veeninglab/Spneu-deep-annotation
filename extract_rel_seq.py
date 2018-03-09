"""
By Jelle Slager, j.slager@rug.nl

Output: sequences upstream and/or downstream of coordinates specified in a coordinate file.
Coordinate file should consist of two columns, without header: col1 = coordinate, col2 = strand
To run: python [this file] [genome file] [coordinate file] [start] [end] [output file]
For example, to extract 100 nts upstream of start sites coordinates:
python extract_rel_seq.py genome.fna start_coords.txt -100 -1 promoters.fna
"""

from Bio import SeqIO
import sys

genome_file = sys.argv[1]
coords_file = sys.argv[2]
start = int(sys.argv[3])
end = int(sys.argv[4])
out_file = sys.argv[5]

genome = SeqIO.read(genome_file,'fasta')
genome_size = len(genome)

fin = open(coords_file,'r')
lines = fin.readlines()
fin.close()

fout = open(out_file,'w')

for line in lines:
    data = line.split()
    pos = int(data[0])
    strand = data[1]
    if strand == '+':
        left_coord = pos + start
        right_coord = pos + end
    if strand == '-':
        left_coord = pos - end
        right_coord = pos - start
    if left_coord < 0: 
        leading_seq = genome.seq[((left_coord+genome_size)-1):] + genome.seq[:right_coord]
    elif right_coord > genome_size:
        leading_seq = genome.seq[(left_coord-1):] + genome.seq[:right_coord-genome_size]
    else:
        leading_seq = genome.seq[(left_coord-1):right_coord]
    if strand == '+':
        seq = leading_seq
    if strand == '-':
        seq = leading_seq.reverse_complement()
#    newline = '>{}_{}\n{}\n'.format(pos,strand,seq)
    newline = '{}\n'.format(seq)
    fout.write(newline)

fout.close()