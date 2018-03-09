"""
Input: (i) genome sequence (Fasta); (ii) genome.coords file produced by 'process_annotation.py'
(iv) [output_folder]/temp_files/rRNAs.txt (produced by 'process_annotation.py'); (v) all normalized per-base count files in [output_folder]/output/perbase/ (produced by 'process_sam.py')
Output: (i) High-confidence TSSs; (ii) Lower-confidence TSSs; (iii) All TSS calls without clustering
Command: python TSS_analysis.py [output_folder] [fasta file]	(other required files are hardcoded in the script)
Example: python TSS_analysis.py D39V_annotation D39V.fna

Required file locations:
[output_folder]/temp_files/chromsizes.txt - produced by 'process_sam_nomulti (ran with iterno=1) - contains 1 line: [chromosome name][TAB][chromosome length]
[output_folder]/temp_files/genome.coords
[output_folder]/temp_files/rRNAs.txt
[output_folder]/output/perbase/norm_5pr_minus_coverage.txt
[output_folder]/output/perbase/norm_5pr_plus_coverage.txt
[output_folder]/output/perbase/norm_5pr_minus_start_count.txt
[output_folder]/output/perbase/norm_5pr_plus_start_count.txt
[output_folder]/output/perbase/norm_cont_minus_coverage.txt
[output_folder]/output/perbase/norm_cont_plus_coverage.txt
[output_folder]/output/perbase/norm_cont_minus_start_count.txt
[output_folder]/output/perbase/norm_cont_plus_start_count.txt
[output_folder]/output/perbase/norm_cont_minus_end_count.txt
[output_folder]/output/perbase/norm_cont_plus_end_count.txt
"""

from __future__ import division
import glob
import os
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO
np.seterr(divide='ignore', invalid='ignore')

out_folder = sys.argv[1]
fasta_file = sys.argv[2]

genome = SeqIO.read(fasta_file,'fasta')
chromsize_file = open(out_folder+'/temp_files/chromsizes.txt','r')
chromsize_lines = chromsize_file.readlines()
chromsize_file.close()

ann_file = open(out_folder+'/temp_files/genome.coords','r')
ann_lines = ann_file.readlines()
ann_file.close()

ann_list = pd.DataFrame(columns=['name','left','right','strand'])
for line in ann_lines:
	data = line.split()
	name = data[0]
	start = int(data[1])
	end = int(data[2])
	left = min(start,end)
	right = max(start,end)
	strand = '+' if start == left else '-'
	ann_list = ann_list.append({'name':name,'left':left,'right':right,'strand':strand},ignore_index=True)

ann_list = ann_list.sort_values(by='left')
ann_list = ann_list.reset_index(drop=True)

genome_name = chromsize_lines[0].split()[0]
genome_size = int(chromsize_lines[0].split()[1])

strands = ['plus','minus']
datanames = []
datasets = {}

for path in glob.glob(out_folder+'/output/perbase/norm*'):
	filename = os.path.basename(path)
	dataname = '_'.join(filename.split('.')[0].split('_')[1:4])
	datanames.append(dataname)
	file = open(path,'r')
	lines = file.readlines()
	file.close()
	datasets[dataname] = np.zeros((1,genome_size))
	for line in lines:
		nt = int(line.split()[0])
		count = float(line.split()[1])
		datasets[dataname][0,nt-1] += count

# Define rRNA regions
rRNA_file = open(out_folder+'/temp_files/rRNAs.txt','r')
rRNA_lines = rRNA_file.readlines()
rRNA_file.close()

rRNAs = {}
rRNAs['plus'] = np.ones((1,genome_size))
rRNAs['minus'] = np.ones((1,genome_size))

for line in rRNA_lines:
	buffer = 30
	specs = line.split()
	left = int(specs[0])
	right = int(specs[1])
	strand = 'plus' if specs[2] == '+' else 'minus'
	if strand == '+':
		left -= buffer
	if strand == '-':
		right += buffer
	if left-1 < 0:
		rRNAs[strand][0,genome_size+left-1:genome_size] = np.zeros((1,1-left))[0]
		rRNAs[strand][0,0:right] = np.zeros((1,right))[0]
	elif right > genome_size:
		rRNAs[strand][0,0:right] = np.zeros((1,right))[0]
		rRNAs[strand][0,(left-1):genome_size] = np.zeros((1,genome_size-left+1))[0]
	else:
		rRNAs[strand][0,(left-1):right] = np.zeros((1,right-left+1))[0]
		
### High confidence TSS: enriched in 5'-enriched sample
# Calculate ratios
ratio_TSS = {}
for strand in strands:
	ratio_TSS[strand] = abs(datasets['5pr_{}_start'.format(strand)][0] / datasets['cont_{}_start'.format(strand)][0])

TSS_cut_ratio = 2.5
cov_cut = 2
too_close = 10
TSSs = {}
TSSs_unclust = {}

for i in range(genome_size):
	for strand in strands:
		if rRNAs[strand][0,i] == 0:
			continue
		if ratio_TSS[strand][i] > TSS_cut_ratio and abs(datasets['5pr_{}_start'.format(strand)][0][i]) > cov_cut:  
			nt = i+1 # Convert from zero_based to normal counting
			value = datasets['5pr_{}_start'.format(strand)][0][i]
			# For unclustered TSSs
			if nt in TSSs_unclust:
				if abs(value) > abs(TSSs_unclust[nt]):
					TSSs_unclust[nt] = value
			else:
				TSSs_unclust[nt] = value
			# For clustered TSSs
			if i-too_close < 0:
				local_ratios = list(ratio_TSS[strand][i-too_close+genome_size:genome_size]) + list(ratio_TSS[strand][0:i+too_close+1])
				local_values = list(abs(datasets['5pr_{}_start'.format(strand)][0][i-too_close+genome_size:genome_size])) + list(abs(datasets['5pr_{}_start'.format(strand)][0][0:i+too_close+1]))
			if i+too_close+1 > genome_size:
				local_ratios = list(ratio_TSS[strand][i-too_close:genome_size]) + list(ratio_TSS[strand][0:i+too_close+1-genome_size])
				local_values = list(abs(datasets['5pr_{}_start'.format(strand)][0][i-too_close:genome_size])) + list(abs(datasets['5pr_{}_start'.format(strand)][0][0:i+too_close+1-genome_size]))
			else:
				local_ratios = list(ratio_TSS[strand][i-too_close:i+too_close+1])
				local_values = list(abs(datasets['5pr_{}_start'.format(strand)][0][i-too_close:i+too_close+1]))
			for j in range(2*too_close+1):
				if local_ratios[j] < TSS_cut_ratio:
					local_values[j] = 0
			local_max = max(local_values)
			max_index = [k for k, l in enumerate(local_values) if l == local_max]
			if (strand == 'plus' and min(max_index) == too_close) or (strand == 'minus' and max(max_index) == too_close):
				if nt in TSSs:
					if abs(value) > abs(TSSs[nt]):
						TSSs[nt] = value
				else:
					TSSs[nt] = value

# Write unclustered and clustered TSSs
coordlist_unclust = []
coordlist_clust = []

for TSS in TSSs_unclust:
	coordlist_unclust.append(TSS)

for TSS in TSSs:
	coordlist_clust.append(TSS)

coordlist_unclust.sort()
coordlist_clust.sort()

TSS_fout_unclust = {}
TSS_fout_clust = {}
TSS_fout_LC = {}
for strand in strands:
	TSS_fout_unclust[strand] = open(out_folder+'/output/bedgraph/TSS_unclust_{}.bedgraph'.format(strand),'w')
	TSS_fout_clust[strand] = open(out_folder+'/output/bedgraph/HC_TSS_{}.bedgraph'.format(strand),'w')
	TSS_fout_LC[strand] = open(out_folder+'/output/bedgraph/LC_TSS_{}.bedgraph'.format(strand),'w')

for coord in coordlist_unclust:
	nt = coord
	value = TSSs_unclust[nt]
	strand = 'plus' if value > 0 else 'minus'
	newline = '{}\t{}\t{}\t{}\n'.format(genome_name,nt-1,nt,value)
	TSS_fout_unclust[strand].write(newline)

all_TSS_list = pd.DataFrame(columns = ['conf','nt','strand','value'])
HC_list = []
for coord in coordlist_clust:
	nt = coord
	value = TSSs[nt]
	### If signal directly upstream in 5' sample is stronger, choose that as TSS
	index_temp = nt-1
	value_temp = value
	if value > 0:
		strand = 'plus'
		while datasets['5pr_plus_start'][0][index_temp-1] > value_temp:
			index_temp -= 1
			value_temp = datasets['5pr_plus_start'][0][index_temp]
	if value < 0:
		strand = 'minus'
		while datasets['5pr_minus_start'][0][index_temp+1] < value_temp:
			index_temp += 1
			value_temp = datasets['5pr_minus_start'][0][index_temp]
	nt=index_temp+1
	value = value_temp
	id = '{}_{}'.format(strand,nt)
	HC_list.append(id)
	all_TSS_list = all_TSS_list.append({'conf':'HC','nt':nt,'strand':strand,'value':abs(value*1e9)}, ignore_index=True) # Note that value is artificially higher to give HC preference over LC
	newline = '{}\t{}\t{}\t{}\n'.format(genome_name,nt-1,nt,value)
	TSS_fout_clust[strand].write(newline)

for strand in strands:
	TSS_fout_unclust[strand].close()				
	TSS_fout_clust[strand].close()

### Low confidence TSS: strong peak in control sample and a -10 TATAAT (max. 1 mismatch) 
starts_plus = pd.DataFrame({'nt': [i+1 for i in range(genome_size)],'strand': ['plus'] * genome_size, 'value': abs(datasets['cont_plus_start'][0])},index = range(genome_size))
starts_minus = pd.DataFrame({'nt': [i+1 for i in range(genome_size)],'strand': ['minus'] * genome_size, 'value': abs(datasets['cont_minus_start'][0])},index = range(genome_size,2*genome_size))
all_starts = pd.concat([starts_plus,starts_minus])

sorted_starts = all_starts.sort_values(by='value',ascending=0)
sorted_starts.index = range(0,len(sorted_starts))

motif = 'TATAAT'
def diff_bases(a,b):
    return sum ( a[i] != b[i] for i in range(len(a)))
LC_list = pd.DataFrame()
for index,row in sorted_starts.iterrows():
	nt = row['nt']
	ind = nt - 1
	strand = row['strand']
	value = row['value']
	if value < 10:
		break
	if '{}_{}'.format(strand,nt) in HC_list:
		continue
	if rRNAs[strand][0,nt-1] == 0:
		continue
	motif_hits = 0
	if strand == 'plus':
		upstream = str(genome.seq[nt-15-1:nt-5-1])
	if strand == 'minus':
		upstream = str(genome.seq[nt+5:nt+15].reverse_complement())
	for i in range(len(upstream)-(len(motif)-1)):
			subseq = upstream[i:i+len(motif)]
			mm_count = diff_bases(motif,subseq)
			if mm_count <=1:
				motif_hits += 1
	if motif_hits == 0:
		continue
	if ind-too_close < 0:
		local_values = list(abs(datasets['cont_{}_start'.format(strand)][0][ind-too_close+genome_size:genome_size])) + list(abs(datasets['cont_{}_start'.format(strand)][0][0:ind+too_close+1]))
	if ind+too_close+1 > genome_size:
		local_values = list(abs(datasets['cont_{}_start'.format(strand)][0][ind-too_close:genome_size])) + list(abs(datasets['cont_{}_start'.format(strand)][0][0:ind+too_close+1-genome_size]))
	else:
		local_values = list(abs(datasets['cont_{}_start'.format(strand)][0][ind-too_close:ind+too_close+1]))
	local_max = max(local_values)
	max_index = [k for k, l in enumerate(local_values) if l == local_max]
	if (strand == 'plus' and min(max_index) == too_close) or (strand == 'minus' and max(max_index) == too_close):	
		LC_list = LC_list.append(row)

LC_sorted = LC_list.sort_values(by='nt')
excl_list = [] # added at a later stage, based on manual curation. Leave empty for uncurated results
#excl_list = [65522,151367,208340,210647,210808,324116,402745,411547,443968,562663,566469,604839,644349,647231,694912,820542,826305,863299,1000776,1045181,1175129,1414496,1422051,1567017,1595904,1598396,1686995,1708145,1760037,1761661,1808008,1836434,1859450,2034063,2042192,2042222,2043657]

for index,row in LC_sorted.iterrows():
	nt = row['nt']
	ind = nt - 1
	strand = row['strand']
	value = row['value'] if strand == 'plus' else -1*row['value']	
	all_TSS_list = all_TSS_list.append({'conf':'LC','nt':nt,'strand':strand,'value':abs(value)}, ignore_index=True)
	if int(nt) not in excl_list: # Ugly way to get rid of filtered LC_TSS. Still need to script this
		newline = '{}\t{}\t{}\t{}\n'.format(genome_name,int(nt-1),int(nt),value)
		TSS_fout_LC[strand].write(newline)

for strand in strands:
	TSS_fout_LC[strand].close()