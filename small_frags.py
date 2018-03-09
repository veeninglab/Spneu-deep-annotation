"""
This script analyzes fragment size distributions to identify overrepresented fragments as putative sRNAs. 
Output: list of enriched fragments, with their genomic context. Additionally: a bed-file with the candidates visualization in Jbrowse in intermediate stage of analysis.
Command: python small_frags.py [output_folder] [GFF_file] iter=XXX
After the first iteration (iter=1), candidates are added to the 'small_features_N.txt' file.
Subsequently, a new library-wide distribution is produced by running 'process_sam_nomulti.py' again, followed by another iteration of this script (iter=2).

Required files:
[output_folder]/temp_files/chromsizes.txt - produced by 'process_sam_nomulti.py'
[output_folder]/temp_files/small_features_1.txt - produced by 'process_annotation.py'
[output_folder]/output/fragment_sizes/{}_fragsize_no_smallfeats.txt - produced by 'process_sam_nomulti.py'
[output_folder]/temp_files/cont_peaks_in_term.txt - produced by 'terminator_analysis.py'
"""

from __future__ import division
import numpy as np
import math
import sys
from collections import defaultdict

out_folder = sys.argv[1]

strands = ['plus','minus']
samples = ['cont']

binsize = 10
expand_rna = 10
iterno = int(sys.argv[3].split('=')[1])

# Read chromosome information
chromsize_file = open(out_folder+'/temp_files/chromsizes.txt','r')
chromsize_lines = chromsize_file.readlines()
chromsize_file.close()

genome_name = chromsize_lines[0].split()[0]
genome_size = int(chromsize_lines[0].split()[1])

## Read annotation
gff_file = open(sys.argv[2],'r')
ann_lines = gff_file.readlines()
gff_file.close()

pos_feat = defaultdict(lambda: defaultdict(set))

fcounts = {}
feat_types = ['tRNA','rRNA','CDS','ncRNA','pseudogene','repeat_region','repeat_unit']
for feat_type in feat_types:
	fcounts[feat_type] = 0

ann_details = {}
for ann_line in ann_lines:
	if ann_line[0] == '#':
		continue	
	line = ann_line.split('\t')
	ftype = line[2]
	strand = line[6]
	opp_strand = '-' if strand == '+' else '+'
	if ftype in feat_types:
		fcounts[ftype] += 1
		if 'Name=' in line[8]:
			fid = line[8].split('Name=')[1].split(';')[0]
		else:
			fid = line[8].split('ID=')[1].split(';')[0]
		left = int(line[3]) - expand_rna if (ftype in ['rRNA','tRNA'] and strand == '-') else int(line[3])
		right = int(line[4]) + expand_rna if (ftype in ['rRNA','tRNA'] and strand == '+') else int(line[4])
		ann_details[fid] = [ftype,left,right]
		for i in range(left, right+1): # Don't use zero-based coordinates here
			pos_feat[strand][i].add(fid)
			if ftype in ['repeat_region','repeat_unit']:
				pos_feat[opp_strand][i].add(fid)
				
pos_feat_list = {}
for strand, data in pos_feat.items():
	list_of_sets = []
	for k in range(genome_size):
		list_of_sets.append(list(data[k]))
	pos_feat_list[strand] = list_of_sets

# Read small_feature data
smallfeat_file = open(out_folder+'/temp_files/small_features_{}.txt'.format(iterno),'r')
smallfeat_lines = smallfeat_file.readlines()
smallfeat_file.close()
smallfeats = {}
smallfeats['+'] = np.zeros((1,genome_size))
smallfeats['-'] = np.zeros((1,genome_size))

smallfeat_ext_file = open(out_folder+'/temp_files/small_features_{}.txt'.format(iterno+1),'w')
for line in smallfeat_lines:
	smallfeat_ext_file.write(line)
	frontbuffer = 30
	specs = line.split()
	left = int(specs[0])
	right = int(specs[1])
	strand = specs[2]
	feat_type = specs[3]
	if feat_type == 'tRNA':
		afterbuffer = 10
	else:
		afterbuffer = 0
	
	if strand == '+':
		left -= frontbuffer
		right += afterbuffer
	if strand == '-':
		left -= afterbuffer
		right += frontbuffer	
	
	if left-1 < 0:
		smallfeats[strand][0,genome_size+left-1:genome_size] = np.ones((1,1-left))[0]
		smallfeats[strand][0,0:right] = np.ones((1,right))[0]
	elif right > genome_size:
		smallfeats[strand][0,0:right-genome_size] = np.ones((1,right-genome_size))[0]
		smallfeats[strand][0,(left-1):genome_size] = np.ones((1,genome_size-left+1))[0]
	else:
		smallfeats[strand][0,(left-1):(right)] = np.ones((1,right-left+1))[0]

# Begin analysis		
fragsizes = {}
end_counts = {}
for sample in samples:
	all_fragsize_file = open(out_folder+'/output/fragment_sizes/{}_fragsize_no_smallfeats.txt'.format(sample),'r')
	all_fragsize_lines = all_fragsize_file.readlines()
	all_fragsize_file.close()
	
	fout = open(out_folder+'/output/fragment_sizes/{}_small_RNA_elements.txt'.format(sample),'w')
	fout_bed = open(out_folder+'/output/bed/{}_sRNA_candidates.bed'.format(sample),'w')
	firstline = '\t'.join(['start','end','strand','enrichment','term_eff','occ_freq','classification','sense_overlap','antisense_overlap'])+'\n'
	fout.write(firstline)
	all_fragsizes = np.zeros((len(all_fragsize_lines)-1,1))
	for line in all_fragsize_lines[1:]: 
		data = line.split()
		fraglen = int(data[0])
		freq = int(float(data[1]))
		all_fragsizes[fraglen-1,0] = freq
		max_size = fraglen
	
	all_rel_freqs = all_fragsizes / sum(all_fragsizes)
	bin_no = int(math.ceil(max_size/binsize))
	binned_fragsizes = {}
	for l in range(bin_no):
		binpos = (l+0.5) * binsize + 0.5
		bincount = int(sum(all_fragsizes[l*binsize:(l+1)*binsize,0])) if (l+1)*binsize <= max_size - 1 else int(sum(all_fragsizes[l*binsize:,0]))
		binned_fragsizes[binpos] = bincount
		
	fragsizes[sample] = {}
	end_counts[sample] = {}
	
	for strand in strands:
		fragsizes[sample][strand] = {'per_end':{},'per_start':{}}
		per_start_file = open(out_folder+'/output/fragment_sizes/{}_per_start_{}.txt'.format(sample,strand),'r')
		per_start_lines = per_start_file.readlines()
		per_start_file.close()
		for line in per_start_lines:
			data = line.split()
			nt = int(data[0])
			sizes = map(int,data[1].split(','))
			fragsizes[sample][strand]['per_start'][nt] = sizes
		
		per_end_file = open(out_folder+'/output/fragment_sizes/{}_per_end_{}.txt'.format(sample,strand),'r')
		per_end_lines = per_end_file.readlines()
		per_end_file.close()
		for line in per_end_lines:
			data = line.split()
			nt = int(data[0])
			sizes = map(int,data[1].split(','))
			fragsizes[sample][strand]['per_end'][nt] = sizes
		end_counts[sample][strand] = np.genfromtxt(out_folder+'/output/perbase/raw_{}_{}_end_count.txt'.format(sample,strand),dtype=None,names=('nt','count'))
		
	peakfile = open(out_folder+'/temp_files/{}_peaks_in_term.txt'.format(sample),'r')
	peaklines = peakfile.readlines()
	peakfile.close()
	
	sRNA_count = 0
	for peak in peaklines:
		data = peak.split()
		strandsign = data[0]
		strand = 'plus' if strandsign == '+' else 'minus'
		start,end = map(int,data[1:3])
		in_term = data[3]
		counts = abs(end_counts[sample][strand]['count'][start-1:end])
		fragno = abs(sum(counts))
		max_term_nt = start + counts.argmax() # nucleotide with highest termination frequency
		
		starts = []
		for nt in range(start,end+1):
			if nt in fragsizes[sample][strand]['per_end']:
				if strand == 'plus':
					starts += list(nt - np.array(fragsizes[sample][strand]['per_end'][nt]) + 1)
				else:
					starts += list(nt + np.array(fragsizes[sample][strand]['per_end'][nt]) - 1)
		starts = np.array(starts) + 1000 # temporarily add 1000 to prevent that any start sites are negative
		abs_freqs = np.bincount(starts)
		rel_freqs = abs_freqs / fragno
		nonzeros = np.nonzero(rel_freqs)[0]
		dists = abs(nonzeros-1000 - max_term_nt) + 1 # Subtract the artificial 1000 again
		array = np.vstack((nonzeros-1000,dists,rel_freqs[nonzeros])).T # Note that rel_freqs is still built up with the nonzero coordinates that were increased by 1000
		index_max = array[:,2].argmax() # nucleotide from which the highest portion of reads terminated at this peak start
		nt_max = int(array[index_max,0])
		min_size = start - nt_max + 1 if strand == 'plus' else nt_max - end + 1
		max_size = end - nt_max + 1 if strand == 'plus' else nt_max - start + 1

		smallfeat_sum = np.sum(smallfeats[strandsign][0,start-1:end])
		if smallfeat_sum == 0: # Exclude small features when determining fragment size distribution
			all_other_fragsizes = np.array(all_fragsizes)
			all_other_rel_freqs = all_other_fragsizes / sum(all_other_fragsizes)
			enr = array[index_max,2] / all_other_rel_freqs[int(array[index_max,1])-1]
		else:
			enr = array[index_max,2] / all_rel_freqs[int(array[index_max,1])-1] # ratio between observed occurrence of this fragment size and the expected occurrence based on all fragments in the entire dataset
			
		if nt_max <= 0:
			nt_max += genome_size
		
		start_sizes = fragsizes[sample][strand]['per_start'][nt_max]
		frags_in = sum(size >= min_size for size in start_sizes)
		frags_term = sum(size >= min_size and size <= max_size for size in start_sizes)
		term_eff = 100 * frags_term / frags_in # fraction of reads starting at the chosen 'TSS' to be terminated in this peak
		
		### Cut-offs: 
		###	Termination eff. (term_eff): 30
		###	Enrichment (enr): 25
		### Fragment coverage (frags_term): 200 for non-HC-terminator peaks, 15 for HC-terminator peaks
		### Then classify based on annotation context: tRNA, in rRNA region, inside CDS etc.
		
		if (in_term == 'True' and enr >= 25 and term_eff >= 25 and frags_term >= 15) or (term_eff >= 25 and enr >= 25 and frags_term >= 200):
			left = min(nt_max,max_term_nt)
			right = max(nt_max,max_term_nt)
			opp_strand = '+' if strand == 'minus' else '-'
			
			flist = []
			for fset in pos_feat_list[strandsign][left:right+1]:
				for el in fset:
					flist.append(el)
			flist = list(set(flist))
			
			flist_as = []
			for fset in pos_feat_list[opp_strand][left:right+1]:
				for el in fset:
					flist_as.append(el)
			flist_as = list(set(flist_as))
			temp = list(flist_as)
			for feat in flist_as:
				if 'repeat' in feat:
					temp.remove(feat)
			flist_as = list(temp)
			
			classes = []
			if len(flist) == 0:
				classes.append('New')
			else:
				for feat in flist:
					ftype = ann_details[feat][0]
					fleft = ann_details[feat][1]
					fright = ann_details[feat][2]
					if ftype in ['tRNA','rRNA','ncRNA','pseudogene']:
						classes.append(ftype)
					if 'repeat' in ftype:
						classes.append('repeat')
					if ftype == 'CDS':
						if left <= fleft and right >= fright:
							classes.append('CDS_full')
						if left <= fleft and right < fright:
							classification = 'CDS_start' if strandsign == '+' else 'CDS_end'
							classes.append(classification)
						if left > fleft and right < fright:
							classes.append('CDS_inside')
						if left > fleft and right >= fright:
							classification = 'CDS_start' if strandsign == '-' else 'CDS_end'
							classes.append(classification)
			classes = list(set(classes))
			class_string = ','.join(classes)
			sense_partners = ','.join(flist)
			antisense_partners = ','.join(flist_as)
			sRNA_count += 1
			newline = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(left,right,strandsign,int(enr),int(round(term_eff)),frags_term,class_string,sense_partners,antisense_partners)
			fout.write(newline)
			bedline = '{}\t{}\t{}\tsRNA_{}\t{}\t{}\n'.format(genome_name,left-1,right,sRNA_count,int(round(term_eff)),strandsign)
			fout_bed.write(bedline)
			smallfeat_line = '{}\t{}\t{}\t{}\n'.format(left,right,strandsign,'ncRNA')
			smallfeat_ext_file.write(smallfeat_line)
	fout.close()
	fout_bed.close()

smallfeat_ext_file.close()