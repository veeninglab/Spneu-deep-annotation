"""
Input: genome annotation (gff3 format)
Output: (i) coordinates of rRNAs, (ii) coordinates of all features smaller than 500, (iii) a genome.coords file (needed as input for TransTermHP later on)
Command: python process_annotation.py [gff_file] [output_folder]
Example: python process_annotation.py D39V.gff3 D39V_annotation
"""

import sys
import re

gff_path = sys.argv[1]
out_folder = sys.argv[2]
fin = open(gff_path,'r')
ann_lines = fin.readlines()
fin.close()

fout_rRNA = open(out_folder+'/temp_files/rRNAs.txt','w')
fout_smallfeats = open(out_folder+'/temp_files/small_features_1.txt','w')
fout_coords = open(out_folder+'/temp_files/genome.coords','w')

small_cut = 500
for n in range(len(ann_lines)):
	line = ann_lines[n]
	if line[0] == '#':
		continue
	data = re.split('\t|\n|\r',line)
	genome_name = data[0]
	feat_type = data[2]
	if feat_type == 'rRNA':
		rRNAline = '{}\t{}\t{}\n'.format(data[3],data[4],data[6])
		fout_rRNA.write(rRNAline)
	if feat_type == 'gene':
		next_feat = ann_lines[n+1].split()[2]
		left,right,strand = [int(data[3]),int(data[4]),data[6]]
		attributes = data[8].split(';')
		for attr in attributes:
			if attr.split('=')[0] == 'Name':
				feat_name = attr.split('=')[1]
		if right-left+1 < small_cut and 'pseudo=true' in attributes:
			smallfeatline = '{}\t{}\t{}\t{}\n'.format(left,right,strand,'pseudo')
			fout_smallfeats.write(smallfeatline)
		elif right-left+1 < small_cut:
			smallfeatline = '{}\t{}\t{}\t{}\n'.format(left,right,strand,next_feat)
			fout_smallfeats.write(smallfeatline)
		start = left if strand == '+' else right
		end = right if strand == '+' else left
		coordsline = '{}\t{}\t{}\t{}\n'.format(feat_name,start,end,genome_name)
		fout_coords.write(coordsline)
	if feat_type == 'pseudogene':
		left,right,strand = [int(data[3]),int(data[4]),data[6]]
		if right-left+1 < small_cut:
			smallfeatline = '{}\t{}\t{}\t{}\n'.format(left,right,strand,'pseudo')
			fout_smallfeats.write(smallfeatline)
		start = left if strand == '+' else right
		end = right if strand == '+' else left
		attributes = data[8].split(';')

		for attr in attributes:
			if attr.split('=')[0] == 'Name':
				feat_name = attr.split('=')[1]

		coordsline = '{}\t{}\t{}\t{}\n'.format(feat_name,start,end,genome_name)		

fout_rRNA.close()
fout_smallfeats.close()
fout_coords.close()