"""
This script performs termination peak-calling on fragment end counts from paired-end sequencing data. Additionally, high-confidence terminators 
are defined when such peaks are found in the poly-U tract of predicted stem-loop structures
Output: (i) putative termination peaks (ii) high-confidence terminators (iii) en passant also converts TransTermHP output to bed-format, more convenient for inspection.
Command: python terminator_analysis.py [output_folder]
Required files: 
[output_folder]/temp_files/chromsizes.txt - produced by 'process_sam_nomulti.py'
[output_folder]/output/perbase/raw_cont_plus_coverage.txt - produced by 'process_sam[_nomulti].py'
[output_folder]/output/perbase/raw_cont_minus_coverage.txt - produced by 'process_sam[_nomulti].py'
[output_folder]/output/perbase/raw_cont_plus_end_count.txt - produced by 'process_sam[_nomulti].py'
[output_folder]/output/perbase/raw_cont_minus_end_count.txt - produced by 'process_sam[_nomulti].py'
[output_folder]/temp_files/TransTerm_loose_output.txt - output file with predicted stem-loop structures by TransTermHP
"""

from __future__ import division
import sys
import numpy as np
import math
import pandas as pd
np.seterr(divide='ignore', invalid='ignore')
pd.options.mode.chained_assignment = None  # default='warn'

out_folder = sys.argv[1]

## Parameters to be used in the peak-calling routine ##
bordering_size = 8
min_peak_size = 3
max_peak_size = 12
window_sizes = range(min_peak_size,max_peak_size+1)
max_cut_off_ratio = 10 # in fact, the cut-off ratio for window size 3
drop_factor = 0.75
cut_off_ratios = [max_cut_off_ratio * (drop_factor + (1-drop_factor) * (3/i)) for i in window_sizes] # accommodate for the fact that the average of of a peak will go down when termination peaks are broader
cut_off_reads = 5
noise_min = 3

## Parameters to be used in the clustering routine for terminator predictions
cluster_cut = 5

## Parameters for combining peaks with predictions
peak_midpred_cut = 7 # half of minimum stem loop + loop region
peak_endpred_cut = 10 # maximum distance between end of stem loop and beginning of peak
#######################################################

chromsize_file = open(out_folder+'/temp_files/chromsizes.txt','r')
chromsize_lines = chromsize_file.readlines()
chromsize_file.close()

genome_name = chromsize_lines[0].split()[0]
genome_size = int(chromsize_lines[0].split()[1])

peak_samples = ['cont']
pred_samples = ['loose']
strands = ['plus','minus']

counts = {}
peaks = {}
coverage = {}
for sample in peak_samples:
	print('Analyzing peaks for sample {}'.format(sample))
	fout_peaks = open(out_folder+'/output/bed/{}_term_peaks.bed'.format(sample),'w')
	peak = 'OFF'
	peak_list = []
	counts[sample] = {}
	coverage[sample] = {}
	for strand in strands:
		counts[sample][strand] = np.zeros((genome_size+4 * max_peak_size,1)) # to make it easier to go over boundaries of circular genome
		fin_counts = open(out_folder+'/output/perbase/raw_{}_{}_end_count.txt'.format(sample,strand),'r')
		countlines = fin_counts.readlines()
		fin_counts.close()
		coverage[sample][strand] = np.genfromtxt(out_folder+'/output/perbase/raw_{}_{}_coverage.txt'.format(sample,strand),dtype=None,names=('nt','cov'))
		for line in countlines:
			data = line.split()
			nt = int(data[0])
			value = float(data[1])
			counts[sample][strand][2*max_peak_size+nt-1,0] += value
			counts[sample][strand][0:2*max_peak_size,0] = counts[sample][strand][genome_size:genome_size+2*max_peak_size,0]
			counts[sample][strand][genome_size+2*max_peak_size:,0] = counts[sample][strand][2*max_peak_size:4*max_peak_size,0]
		for j in range(len(window_sizes)):
			window_size = window_sizes[j]
			cut_off_ratio = cut_off_ratios[j]
			for i in range(genome_size):
#			for i in range(100000):
				if i % 100000 == 0:
					print('window_size {}, {}: {}'.format(window_size,strand,i))
				start = 2*max_peak_size + i
				end = start+window_size - 1
				prev_window = abs(counts[sample][strand][start-bordering_size:start,0])
				prev_av = np.sum(prev_window) / bordering_size
				window = abs(counts[sample][strand][start:end+1,0])
				window_av = np.sum(window) / window_size
				next_window = abs(counts[sample][strand][end+1:end+1+bordering_size,0])
				next_av = np.sum(next_window) / bordering_size
				if window_av > prev_av * cut_off_ratio and window_av > next_av * cut_off_ratio and window_av > cut_off_reads:
					is_peak = True
				else:
					is_peak = False
				if is_peak == True:	
					if peak == 'OFF':
						peakstart = start
					peakend = end
					peak = 'ON'
				if is_peak == False or (is_peak == True and i == genome_size - 1): # To also process the last peak if it stretches until the last nucleotide of the genome
					if peak == 'ON':
						while True:
							peak_window = abs(counts[sample][strand][peakstart:peakend+1,0])
							noise_level = max(noise_min,0.01*max(peak_window))
							if (abs(counts[sample][strand][peakstart,0]) <= 0.5 * max_cut_off_ratio * np.average(abs(counts[sample][strand][peakstart-bordering_size:peakstart,0])) and abs(counts[sample][strand][peakstart,0]) < 0.1 * max(peak_window)) or abs(counts[sample][strand][peakstart,0]) <= noise_level:
								peakstart += 1
							else:
								break
						while True:
							peak_window = abs(counts[sample][strand][peakstart:peakend+1,0])
							noise_level = max(noise_min,0.01*max(peak_window))
							if (abs(counts[sample][strand][peakend,0]) <= 0.5 * max_cut_off_ratio * np.average(abs(counts[sample][strand][peakend+1:peakend+1+window_size,0])) and abs(counts[sample][strand][peakend,0]) < 0.1 * max(peak_window)) or abs(counts[sample][strand][peakend,0]) <= noise_level:
								peakend -= 1
							else:
								break
						if peakend-peakstart + 1 >= min_peak_size:
							peak_list.append([peakstart,peakend,strand])
					peak = 'OFF'	
	
	# Write all potential peaks to temp_files
	all_output = open(out_folder+'/temp_files/{}_all_term_peaks.txt'.format(sample),'w')
	for peak in peak_list:
		all_output.write('{}\t{}\t{}\n'.format(peak[0],peak[1],peak[2]))
	all_output.close()
	
	peak_df = pd.DataFrame(peak_list,columns=['start','end','strand'])
	peak_df = peak_df.sort_values(by=['strand','start'])
	peak_df.index = range(0,len(peak_df))
	
	# Deal with overlapping peaks
	overlap_count = 1
	while overlap_count > 0:
		overlap_count = 0
		new_peak_list = []
		for it in range(len(peak_df)):
			curr = it + overlap_count
			next = curr + 1
			if curr > len(peak_df) - 1:
				break
			peak = list(peak_df.loc[curr,:])
			if curr == len(peak_df) - 1:
				new_peak_list.append(peak)
				break
			nextpeak = list(peak_df.loc[next,:])
			strand = peak[2]
			if peak[1] >= nextpeak[0] and strand == nextpeak[2]: # indication overlap between peaks
				newstart = peak[0]
				newend = nextpeak[1]
				if newstart == newend:
					print(newstart,newend)
				window = abs(counts[sample][strand][newstart:newend+1,0])
				noise_level = max(noise_min,0.01*max(window))
				while True:
					if (abs(counts[sample][strand][newstart,0]) <= 0.5 * max_cut_off_ratio * np.average(abs(counts[sample][strand][newstart-max_peak_size:newstart,0])) and abs(counts[sample][strand][newstart,0]) < 0.1 * max(window)) or abs(counts[sample][strand][newstart,0]) <= noise_level: # Drop nucleotides from edges of peak if they're not significantly higher than signals immediately to the left
						newstart += 1
					else:
						break
				while True:
					if (abs(counts[sample][strand][newend,0]) <= 0.5 * max_cut_off_ratio * np.average(abs(counts[sample][strand][newend+1:newend+1+max_peak_size,0])) and abs(counts[sample][strand][newend,0]) < 0.1 * max(window)) or abs(counts[sample][strand][newend,0]) <= noise_level: # Drop nucleotides from edges of peak if they're not significantly higher than signals immediately to the right
						newend -= 1
					else:
						break
				overlap_count += 1 # keeps repeating these steps for all peaks, until no overlapping peaks are found anymore
			else:
				newstart = peak[0]
				newend = peak[1]
			if newend-newstart + 1 >= min_peak_size:
				new_peak_list.append([newstart,newend,strand])
		peak_df = pd.DataFrame(new_peak_list,columns=['start','end','strand'])
	peaks[sample] = peak_df
	peaks[sample] = peaks[sample].sort_values(by='start')
	peaks[sample].index = range(len(peaks[sample]))
	peak_count = 0
	for index, row in peaks[sample].iterrows():
		peak_count += 1
		name = 'PEAK_' + str(peak_count)
		strandsign = '+' if row['strand'] == 'plus' else '-'
		score = max(abs(counts[sample][row['strand']][row['start']:row['end']+1,0]))			
		startbase = row['start']-2*max_peak_size+1 # Transform back to nucleotide count
		endbase = row['end']-2*max_peak_size+1 # Transform back to nucleotide count		
		newline = '\t'.join([genome_name,str(startbase-1),str(endbase),name,str(score),strandsign])+'\n'
		fout_peaks.write(newline)
	peaks[sample][['start','end']] = peaks[sample][['start','end']] - 2*max_peak_size + 1
	fout_peaks.close()

preds = {}
for sample in pred_samples:
	print('Analyzing predictions for sample {}'.format(sample))
	fin = open(out_folder+'/temp_files/TransTerm_{}_output.txt'.format(sample),'r')
	predlines = fin.readlines()
	fin.close()
	predlist = []
	for p in range(len(predlines)):
		line = predlines[p]
		data = line.split()
		if len(data) > 2:
			if data[0] == 'TERM':
				start = int(data[2])
				end = int(data[4])
				middle = (start+end)/2
				length = abs(end-start)+1
				strand = data[5]
				score = int(data[7])
				hairpin = float(data[8])
				nextline = predlines[p+1]
				nextdata = nextline.split()
				loop = len(nextdata[2])
				polyU = nextdata[4].count('T') if strand == '+' else nextdata[0].count('A')
				predlist.append([start,end,middle,length,strand,score,hairpin,loop,polyU])
	preds[sample] = pd.DataFrame(predlist,columns=['start','end','middle','length','strand','score','hairpin','loop','polyU'])
	preds[sample] = preds[sample].sort_values(by=['strand','middle'])
	preds[sample].index = range(len(preds[sample]))
	prev_mid = 1e42 # just a big number
	prev_strand = ''
	cluster_no = -1
	for i in range(len(preds[sample])):
		mid = preds[sample].loc[i,'middle']
		strand = preds[sample].loc[i,'strand']
		if abs(mid-prev_mid) > cluster_cut or strand != prev_strand:
			cluster_no += 1
		preds[sample].loc[i,'cluster'] = cluster_no
		prev_mid = mid
		prev_strand = strand
	cluster_no += 1	
	for m in range(cluster_no):
		cluster_df = preds[sample].loc[preds[sample].loc[:,'cluster'] == m]
		if len(cluster_df) > 1:
			cluster_max = cluster_df.loc[:,'score'].max(axis=0)
			for index, row in cluster_df.iterrows():
				if row['score'] < cluster_max: # Keep highest scoring terminator from cluster
					preds[sample] = preds[sample].drop(index)
					cluster_df = cluster_df.drop(index)
		if len(cluster_df) > 1:
			cluster_smallest = cluster_df.loc[:,'length'].min(axis=0)
			for index, row in cluster_df.iterrows():
				if row['length'] > cluster_smallest: # in case of a draw in score, keep the smallest in length
					preds[sample] = preds[sample].drop(index)
					cluster_df = cluster_df.drop(index)	
	preds[sample] = preds[sample].sort_values(by=['start'])
	preds[sample].index = range(len(preds[sample]))
	fout_preds = open(out_folder+'/output/bed/transterm_{}.bed'.format(sample),'w')
	for index, row in preds[sample].iterrows():
		start = row['start']
		end = row['end']
		name = 'PRED_'+str(index+1)
		score = row['score']
		strand = row['strand']
		newline = '\t'.join([genome_name,str(min(start,end)-1),str(max(start,end)),name,str(score),strand])+'\n'
		fout_preds.write(newline)
	fout_preds.close()	
	
preds_ref = preds['loose'][['start','end','strand','hairpin','loop','polyU']]
newcoords = pd.DataFrame({'start':[],'end':[]})
for index, row in preds_ref.iterrows():
	newcoords.loc[index,['start','end']]=[min(row['start'],row['end']),max(row['start'],row['end'])]
	
preds_ref['newstart'] = list(newcoords['start'])
preds_ref['newend'] = list(newcoords['end'])
preds_ref['middle'] = list(0.5 * preds_ref['start'] + 0.5 * preds_ref['end'])
preds_list = []
for index, row in preds_ref.iterrows():
	preds_list.append([row['strand'],row['newstart'],row['newend'],row['middle'],row['hairpin'],row['loop'],row['polyU']])

peak_lists = {}
peak_in_term = {}
for sample in peak_samples:
	peak_lists[sample] = []
	for index, row in peaks[sample].iterrows():
		strand = '+' if row['strand'] == 'plus' else '-'
		peak_lists[sample].append([strand,row['start'],row['end']])
	peak_in_term[sample] = [False] * len(peak_lists[sample])
	
term_HC = {}
term_HC_max = {}
term_HC_specs = {}
for sample in peak_samples:
	term_count = 0
	term_max_list = []
	print('Combining peaks with predictions for sample {}'.format(sample))
	term_HC[sample] = []
	term_HC_max[sample] = []
	term_HC_specs[sample] = []
	fout_term = open(out_folder+'/output/bed/HC_term_{}.bed'.format(sample),'w')
	fout_term_max = open(out_folder+'/output/bed/term_class/HC_term_max_{}.bed'.format(sample),'w')
	fout_term_specs = open(out_folder+'/output/bed/term_class/HC_term_specs_{}.txt'.format(sample),'w')
	fout_peaks_quick = open(out_folder+'/temp_files/{}_peaks_in_term.txt'.format(sample),'w')
	for pred in preds_list:
		[pred_strand,pred_start,pred_end,pred_mid,hairpin,loop,polyU] = pred
		if pred_strand == '+':
			for i in range(len(peak_lists[sample])-1,-1,-1):
				peak = peak_lists[sample][i]
				[peak_strand,peak_start,peak_end] = peak
				if peak_strand == '+':
					if peak_start - pred_end <= peak_endpred_cut and peak_start - pred_mid >= peak_midpred_cut:
						term_count += 1
						bef_cov = coverage[sample]['plus']['cov'][peak_start-1] # Coverage at the first base of the peak
						stopped = sum(counts[sample]['plus'][2*max_peak_size+peak_start-1:2*max_peak_size+peak_end,0]) # Number of fragments stopped by this terminator. Note that the counts dataframe has its counts shifted by 2*window_size
						score = int(round(100 * (stopped / bef_cov)))
						term_HC[sample].append([term_count,pred_start-1,peak_end,pred_strand,score])
						peak_counts = abs(counts[sample]['plus'][peak_start+2*max_peak_size-1:peak_end+2*max_peak_size,0])
						term_max = peak_start + np.argmax(peak_counts)
						term_HC_max[sample].append([term_count,term_max-1,term_max,pred_strand,score])
						term_HC_specs[sample].append([term_count,pred_strand,term_max,hairpin,loop,polyU])
						peak_in_term[sample][i] = True
						break
		if pred_strand == '-':	
			for i in range(len(peak_lists[sample])):
				peak = peak_lists[sample][i]
				[peak_strand,peak_start,peak_end] = peak
				if peak_strand == '-':
					if pred_start - peak_end <= peak_endpred_cut and pred_mid - peak_end >= peak_midpred_cut:
						term_count += 1
						bef_cov = abs(coverage[sample]['minus']['cov'][peak_end-1]) # Coverage at the first base of the peak
						stopped = abs(sum(counts[sample]['minus'][2*max_peak_size+peak_start-1:2*max_peak_size+peak_end,0])) # Number of fragments stopped by this terminator. Note that the counts dataframe has its counts shifted by 2*window_size
						score = int(round(100 * (stopped / bef_cov)))
						term_HC[sample].append([term_count,peak_start-1,pred_end,pred_strand,score])
						peak_counts = abs(counts[sample]['minus'][peak_start+2*max_peak_size-1:peak_end+2*max_peak_size,0])
						term_max = peak_start + np.argmax(peak_counts)
						term_HC_max[sample].append([term_count,term_max-1,term_max,pred_strand,score])
						term_HC_specs[sample].append([term_count,pred_strand,term_max,hairpin,loop,polyU])
						peak_in_term[sample][i] = True
						break
	
	for index,start,end,strand,score in term_HC[sample]:
		newline = '{}\t{}\t{}\tTERM_{}\t{}\t{}\n'.format(genome_name,int(start),int(end),index,score,strand)
		fout_term.write(newline)
	fout_term.close()
	
	for index,start,end,strand,score in term_HC_max[sample]:
		newline = '{}\t{}\t{}\tTERM_{}\t{}\t{}\n'.format(genome_name,int(start),int(end),index,score,strand)
		fout_term_max.write(newline)
	fout_term_max.close()
	
	fout_term_specs.write('index\tpos\thp_score\tloop_length\tU_freq\n')
	for index,strand,pos,hairpin,loop,polyU in term_HC_specs[sample]:
		newline = 'TERM_{}\t{}\t{}\t{}\t{}\t{}\n'.format(index,strand,pos,hairpin,loop,polyU)
		fout_term_specs.write(newline)
	fout_term_specs.close()
	
	for i in range(len(peak_lists[sample])):
		[peak_strand,peak_start,peak_end] = peak_lists[sample][i]
		in_term = peak_in_term[sample][i]
		newline = '{}\t{}\t{}\t{}\n'.format(peak_strand,peak_start,peak_end,in_term) 
		fout_peaks_quick.write(newline)
	fout_peaks_quick.close()