"""
This script loops over all mapped reads in a SAM file. In this version of the script (process_sam.py), multimappers are INcluded
Output: (i) chromsizes.txt, (ii) raw and normalized start counts, end counts (if paired-end) and coverage files, both per-base counts and bedgraph format; (iii) fragment size distribution of all reads that 
			do not map on small features, (iv) fragment sizes sorted by both start position and end position (e.g. fragment sizes of all fragments terminated at residue X).
Command: python bin/process_sam.py [output_folder] iter=1
Note that this version of the script needs to be run with 'iter=1'
 
Required files: 
[output_folder]/temp_files/genome_overlap.fna - genome sequence with first 1kb copied at the end (produced by 'add_overlap.py')
[output_folder]/temp_files/rRNAs.txt - rRNA coordinates (produced by 'process_annotation.py')
[output_folder]/temp_files/small_features_1.txt - coordinates of all features smaller than 500 bps (produced by 'process_annotation.py')
[output_folder]/output/sam/cont.sam - SAM file from mapping paired-end sequenced control sample
[output_folder]/output/sam/5pr.sam - SAM file from mapping single-end sequenced 5'-enriched sample
"""

from __future__ import division
import sys
import re
import numpy as np
import time
import math

out_folder = sys.argv[1]
iterno = int(sys.argv[2].split('=')[1])
# Input and output files
fasta_file = open(out_folder+'/temp_files/genome_overlap.fna','r')
fasta_lines= fasta_file.readlines()
fasta_file.close()

rRNA_file = open(out_folder+'/temp_files/rRNAs.txt','r')
rRNA_lines = rRNA_file.readlines()
rRNA_file.close()

smallfeat_file = open(out_folder+'/temp_files/small_features_{}.txt'.format(iterno),'r')
smallfeat_lines = smallfeat_file.readlines()
smallfeat_file.close()

if iterno == 1:
	samples = ['cont','5pr']
	#samples = ['cont'] # Temporary, for degradation analysis
else:
	samples = ['cont']

seq_mode = {'cont':'PE','3pr':'PE','5pr':'SE'}

filelength = {}
totallines = 0

samlines = {}
for sample in samples:
	samfile = open(out_folder+'/output/sam/{}.sam'.format(sample),'r')
	samlines[sample] = samfile.readlines()
	samfile.close()
	filelength[sample] = len(samlines[sample])
	totallines += filelength[sample]

out_subfolder = out_folder+'/output/'

# Extract genome size, name etc.
chr = fasta_lines[0].split()[0][1:]

seqcount = 0
for line in fasta_lines[1:]:
	seqcount += len(line.split()[0])
genome_size = seqcount - 1000

if iterno == 1:
	chromsizes = open(out_folder+'/temp_files/chromsizes.txt','w')
	newline = '\t'.join([chr,str(genome_size)]) + '\n'
	chromsizes.write(newline)
	chromsizes.close()

rRNAs = {}
rRNAs['+'] = np.ones((1,genome_size))
rRNAs['-'] = np.ones((1,genome_size))

smallfeats = {}
smallfeats['+'] = np.zeros((1,genome_size))
smallfeats['-'] = np.zeros((1,genome_size))

for line in rRNA_lines:
	buffer = 30
	specs = line.split()
	left = int(specs[0])
	right = int(specs[1])
	strand = specs[2]
	
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
	
for line in smallfeat_lines:
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

# Create dataframes
all_output = {}
fragment_sizes = {}
countall = 0
time_estimate = 0
for sample in samples:
	all_output[sample] = {'raw':{},'norm':{}}
	
	count = 0
	# Read and analyze data
	if seq_mode[sample] == 'SE':
		mode = 'SE'
		tracks = ['plus_start_count','minus_start_count','plus_coverage','minus_coverage']
	if seq_mode[sample] == 'PE':
		mode = 'PE'
		tracks = ['plus_start_count','plus_end_count','minus_start_count','minus_end_count','plus_coverage','minus_coverage']
		frag_tracks = ['per_start_plus','per_start_minus','per_end_plus','per_end_minus']
		fragment_sizes[sample] = {'all':np.zeros((1,1001))}
		for frag_track in frag_tracks:
			fragment_sizes[sample][frag_track] = {}
			for i in range(genome_size):
				fragment_sizes[sample][frag_track][i] = []
	for track in tracks:
		all_output[sample]['raw'][track] = np.zeros((1,genome_size))
	
	old_pair = ''
	
	perc_step = 10
	perc_count = perc_step
	
	start_time = time.time()
		
	print('Started reading "{}" sam-file'.format(sample))
	for line in samlines[sample]:
		count += 1
		countall += 1
		
		if line[0] == '@':
			continue
#		if count > 1500000: # for testing script
#			break # for testing script		
		
		# progress_counter
		if countall < 1e6:
			if 100.0 * count / filelength[sample] > perc_count:
				print(sample+': processed {} out of {} lines ({}%).'.format(count,filelength[sample],perc_count))
				perc_count += perc_step
				
		if countall == 1e6:
			timeperline = (time.time() - start_time) / 1e6
		
		if countall >= 1e6:
			if 100.0 * count / filelength[sample] > perc_count:
				rem_time = timeperline * (totallines-countall)
				print(sample+': processed {} out of {} lines ({}%). Estimated runtime left: {}'.format(count,filelength[sample],perc_count,time.strftime("%H:%M:%S", time.gmtime(rem_time))))
				perc_count += perc_step
		data = line.split()
		
		# Filter reads
		flag = "{0:012b}".format(int(data[1]))
		
		if mode == 'SE':
			if flag[11] == '1':
				print(sample +': data seems to be paired-end, while single-end was expected')
				break
			
			if flag[:3] != '000':
				print(line + sample +': check the meaning of the sam-flag on this line')
				continue	# some weird error took place
			
			if flag[9] == '1':
				continue	# read unmapped
			
		if mode == 'PE':
			if flag[11] == '0':
				print(sample + ': data seems to be single-end, while paired-end was expected')
				break	# not paired-end, should not happen
			
			if flag[10] == '0':
				continue	# not in proper pair	
			
			if flag[:3] != '000':
				print(line + sample +': check the meaning of the sam-flag on this line')
				continue	# some weird error took place
				
			if abs(int(data[8])) > 1000:
				continue	# mates of pair too far apart

		# Determine reliable read length
		cigar = re.findall(r'([0-9]+)([MIDNSHPX=])',data[5])
		readlength = 0
			
		for item in cigar:
			if item[1] in 'MX=DN':
				readlength += int(item[0])
	
		# Determine strand and location (zero-based)
		if mode == 'SE':
			if flag[7] == '1':
				strand = '-'
				start = int(data[3]) + (readlength - 1) - 1
				end = start - (readlength - 1)
			else:
				strand = '+'
				start = int(data[3]) - 1
				end = start + (readlength - 1)
			
			# Fold back ori-overlap
			if start > genome_size - 1:
				start = start - genome_size
			if end > genome_size - 1:
				end = end - genome_size
				
			# Counting
			if strand == '+':
				all_output[sample]['raw']['plus_start_count'][0,start] += 1
				if end > start:
					all_output[sample]['raw']['plus_coverage'][0,start:(end+1)] += np.ones((1,readlength))[0]
				if end < start:
					all_output[sample]['raw']['plus_coverage'][0,start:] += np.ones((1,genome_size - start))[0]
					all_output[sample]['raw']['plus_coverage'][0,:(end+1)] += np.ones((1,end + 1))[0]
			if strand == '-':
				all_output[sample]['raw']['minus_start_count'][0,start] -= 1
				if end < start:
					all_output[sample]['raw']['minus_coverage'][0,(end):(start+1)] -= np.ones((1,readlength))[0]
				if end > start:
					all_output[sample]['raw']['minus_coverage'][0,end:] -= np.ones((1,genome_size - end))[0]
					all_output[sample]['raw']['minus_coverage'][0,:(start+1)] -= np.ones((1,start+1))[0]
		if mode == 'PE':
			new_pair = data[0] # used to determine whether this read forms a pair with the previous one or with the next one
			
			# Determine strand and location (zero-based)
			if flag[5] == '1': # first read of pair
				if flag[7] == '1':
					strand = '-'
					start = int(data[3]) + (readlength - 1) - 1
				else:
					strand = '+'
					start = int(data[3]) - 1
			
			if flag[4] == '1': # second read of pair
				if flag[7] == '0':
					strand = '-' # mate read maps to opposite strand
					end = int(data[3]) - 1
				else:
					strand = '+'
					end = int(data[3]) + (readlength - 1) - 1
									
			if new_pair == old_pair: # Indicates that both start and end have been assigned for this read pair
				fragment_length = abs(end-start) + 1
				
				# Fold back ori-overlap
				if start > genome_size - 1:
					start = start - genome_size
				if end > genome_size - 1:
					end = end - genome_size

				if strand == '+':
					if end > start:
						smallfeat_sum = np.sum(smallfeats[strand][0,start:end+1])
					if end < start:
						smallfeat_sum = np.sum(smallfeats[strand][0,start:genome_size])
						smallfeat_sum += np.sum(smallfeats[strand][0,0:end+1])
				if strand == '-':
					if end < start:
						smallfeat_sum = np.sum(smallfeats[strand][0,end:start+1])
					if end > start:
						smallfeat_sum = np.sum(smallfeats[strand][0,end:genome_size])
						smallfeat_sum += np.sum(smallfeats[strand][0,0:start+1])
				
				if smallfeat_sum == 0: # Exclude small features when determining fragment size distribution
					fragment_sizes[sample]['all'][0,fragment_length] += 1

				# Counting
				if iterno == 1:
					if strand == '+':
						fragment_sizes[sample]['per_start_plus'][start].append(fragment_length)
						fragment_sizes[sample]['per_end_plus'][end].append(fragment_length)
						all_output[sample]['raw']['plus_start_count'][0,start] += 1
						all_output[sample]['raw']['plus_end_count'][0,end] += 1
						if end > start:
							all_output[sample]['raw']['plus_coverage'][0,start:(end+1)] += np.ones((1,fragment_length))[0]
						if end < start: # meaning the fragment overlaps with the origin
							all_output[sample]['raw']['plus_coverage'][0,start:] += np.ones((1,genome_size - start))[0]
							all_output[sample]['raw']['plus_coverage'][0,:(end+1)] += np.ones((1,end + 1))[0]
					if strand == '-':
						fragment_sizes[sample]['per_start_minus'][start].append(fragment_length)
						fragment_sizes[sample]['per_end_minus'][end].append(fragment_length)
						all_output[sample]['raw']['minus_start_count'][0,start] -= 1
						all_output[sample]['raw']['minus_end_count'][0,end] -= 1
						if end < start:
							all_output[sample]['raw']['minus_coverage'][0,(end):(start+1)] -= np.ones((1,fragment_length))[0]
						if end > start: # meaning the fragment overlaps with the origin
							all_output[sample]['raw']['minus_coverage'][0,end:] -= np.ones((1,genome_size - end))[0]
							all_output[sample]['raw']['minus_coverage'][0,:(start+1)] -= np.ones((1,start+1))[0]
						
			old_pair = new_pair
	
	print('Finished reading "{}" sam-file'.format(sample))
	
	# Normalization, leaving out reads that cover rRNAs
	if iterno == 1:
		norm_read_no = abs(np.dot(rRNAs['+'][0],all_output[sample]['raw']['plus_start_count'][0])) + abs(np.dot(rRNAs['-'][0],all_output[sample]['raw']['minus_start_count'][0]))

		for track in tracks:
			all_output[sample]['norm'][track] = 1e6 * all_output[sample]['raw'][track] / norm_read_no
	
		# Create output files
		print('{}: writing output files'.format(sample))
		for state in ['raw','norm']:
			for track in tracks:
				perbase = open(out_subfolder+'perbase/'+state+'_'+sample+'_'+track+'.txt','w')
				bedgraph = open(out_subfolder+'bedgraph/'+state+'_'+sample+'_'+track+'.bedgraph','w')
			
				count_start = 0
				prev_value = 0
			
				for i in range(genome_size):
					value = all_output[sample][state][track][0,i]
				
					# Perbase files
					perbaseline = str(i+1) + '\t' + str(value) +'\n'
					perbase.write(perbaseline)
				
					# Bedgraph files
					if value != prev_value and i != 0:
						newline = '\t'.join([chr,str(count_start),str(count_end),str(prev_value)])+'\n' 
						count_start = i
						bedgraph.write(newline)
					count_end = i+1
					prev_value = value
		
					if i == genome_size - 1:
						if value == prev_value:
							newline = '\t'.join([chr,str(count_start),str(count_end),str(prev_value)])+'\n'
							bedgraph.write(newline)
						if value != prev_value:
							newline = '\t'.join([chr,str(count_start),str(count_end),str(value)])+'\n'
							bedgraph.write(newline)
				bedgraph.close()
				perbase.close()
	
	if mode == 'PE':
		f_allfrag = open(out_subfolder+'fragment_sizes/'+sample+'_fragsize_no_smallfeats.txt','w')
		firstline = '\t'.join(['frag_size','frequency'])+'\n'
		f_allfrag.write(firstline)
		for i in range(1,1001):
			newline = '{}\t{}\n'.format(i,int(fragment_sizes[sample]['all'][0,i]))
			f_allfrag.write(newline)
		f_allfrag.close()
		
		if iterno == 1:
			for frag_track in frag_tracks:
				f_fragtrack = open(out_subfolder+'fragment_sizes/'+sample+'_'+frag_track+'.txt','w')
				for i in range(genome_size):
					if len(fragment_sizes[sample][frag_track][i]) > 0:
						newline = '{}\t{}\n'.format(i+1,','.join(map(str,fragment_sizes[sample][frag_track][i])))
						f_fragtrack.write(newline)
				f_fragtrack.close()