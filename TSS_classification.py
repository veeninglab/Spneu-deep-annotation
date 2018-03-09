"""
Output: 
Required files:
coords_[date].txt (here: coords_180302.txt):
		Tab-delimited with header: columns: feature_type, start, stop, strand
TSSs.txt:
	Tab-delimited with header: columns: nt, strand, conf, score [i.e. position, strand, confidence, score as reported by 'TSS_analysis.py') 
"""

from __future__ import division
import numpy as np
import pandas as pd
import re
np.seterr(divide='ignore', invalid='ignore')

date = '180302'

genome_size = 2046572 # change this to the length of the genome of interest
primary_cut = 300
inside_cut = 0 # negative when downstream of start codon is allowed
antisense_cut1 = 0
antisense_cut2 = 100 
antisense_cut3 = 300 # these cut-offs separate aTSSs close to and farther away from a feature
ann_start = 0
feat_TSS = {}
classes = []
pos_classes = ['primary','secondary','internal','aTSS-1','aTSS-2up','aTSS-3up','aTSS-2dn','aTSS-3dn','orphan']

ann_file = open('coords_{}.txt'.format(date),'r')
ann_lines = ann_file.readlines()
ann_file.close()
ann_list = pd.DataFrame(columns=['name','left','right','strand'])
for line in ann_lines[1:]:
    data = re.split('\t|\n',line)
    name = data[5]
    left = int(data[1])
    right = int(data[2])
    strand = data[3]
    ann_list = ann_list.append({'name':name,'left':left,'right':right,'strand':strand},ignore_index=True)
ann_list = ann_list.sort_values(by='left')
ann_list = ann_list.reset_index(drop=True)

TSSfile = open('TSSs.txt','r')
TSSlines = TSSfile.readlines()
TSSfile.close()

all_TSS_list = pd.DataFrame(columns = ['conf','nt','strand','value'])
for line in TSSlines[1:]:
    data = line.split()
    all_TSS_list = all_TSS_list.append({'conf':data[2],'nt':int(data[0]),'strand':data[1],'value':abs(float(data[3]))}, ignore_index=True)

all_TSS_list = all_TSS_list.sort_values(by='nt')
all_TSS_list = all_TSS_list.reset_index(drop=True)
all_TSS_list = all_TSS_list.set_index(all_TSS_list['nt'],verify_integrity=True)

class_loci = {}
for line in TSSlines[1:]:
    class_list = []
    TSS = int(line.split()[0])
    strand = line.split()[1]
    id = str(TSS)
    class_loci[id] = {}
    for cl in pos_classes:
        class_loci[id][cl] = []
    if genome_size-TSS < primary_cut:
        for ann_ind, ann_row in ann_list.iterrows():
            left = int(ann_row['left'])
            left_ol = left + genome_size
            right = int(ann_row['right'])
            right_ol = right + genome_size
            name = ann_row['name']
            if strand == '+':
                if ann_row['strand'] == '+' and ((left - TSS >= inside_cut and left - TSS <= primary_cut) or (left_ol - TSS >= inside_cut and left_ol - TSS <= primary_cut)):
                    if name not in feat_TSS:
                        feat_TSS[name] = [TSS]
                    else:
                        feat_TSS[name].append(TSS)
                if ann_row['strand'] == '+' and ((left - TSS < inside_cut and right >= TSS) or (left_ol - TSS < inside_cut and right_ol >= TSS)):
                    class_list.append('internal')
                    class_loci[id]['internal'].append(name)
                if ann_row['strand'] == '-':
                    if (TSS >= left - antisense_cut1 and TSS <= right + antisense_cut1) or (TSS >= left_ol - antisense_cut1 and TSS <= right_ol + antisense_cut1):
                        class_list.append('aTSS-1')
                        class_loci[id]['aTSS-1'].append(name)
                    if (TSS <= right + antisense_cut2 and TSS > right + antisense_cut1) or (TSS <= right_ol + antisense_cut2 and TSS > right_ol + antisense_cut1):
                        class_list.append('aTSS-2up')
                        class_loci[id]['aTSS-2up'].append(name)
                    if (TSS <= right + antisense_cut3 and TSS > right + antisense_cut2) or (TSS <= right_ol + antisense_cut3 and TSS > right_ol + antisense_cut2):
                        class_list.append('aTSS-3up')
                        class_loci[id]['aTSS-3up'].append(name)
                    if (TSS >= left - antisense_cut2 and TSS < left - antisense_cut1) or (TSS >= left_ol - antisense_cut2 and TSS < left_ol - antisense_cut1):
                        class_list.append('aTSS-2dn')
                        class_loci[id]['aTSS-2dn'].append(name)
                    if (TSS >= left - antisense_cut3 and TSS < left - antisense_cut2) or (TSS >= left_ol - antisense_cut3 and TSS < left_ol - antisense_cut2):
                        class_list.append('aTSS-3dn')
                        class_loci[id]['aTSS-3dn'].append(name)
            if strand == '-':
                if ann_row['strand'] == '-' and ((TSS - right >= inside_cut and TSS - right <= primary_cut) or (TSS - right_ol >= inside_cut and TSS - right_ol <= primary_cut)):
                    if name not in feat_TSS:
                        feat_TSS[name] = [TSS]
                    else:
                        feat_TSS[name].append(TSS)
                if ann_row['strand'] == '-' and ((TSS - right < inside_cut and TSS >= left) or (TSS - right_ol < inside_cut and TSS >= left_ol)):
                    class_list.append('internal')
                    class_loci[id]['internal'].append(name)
                if ann_row['strand'] == '+':
                    if (TSS >= left - antisense_cut1 and TSS <= right + antisense_cut1) or (TSS >= left_ol - antisense_cut1 and TSS <= right_ol + antisense_cut1):
                        class_list.append('aTSS-1')
                        class_loci[id]['aTSS-1'].append(name)
                    if (TSS <= right + antisense_cut2 and TSS > right + antisense_cut1) or (TSS <= right_ol + antisense_cut2 and TSS > right_ol + antisense_cut1):
                        class_list.append('aTSS-2dn')
                        class_loci[id]['aTSS-2dn'].append(name)
                    if (TSS <= right + antisense_cut3 and TSS > right + antisense_cut2) or (TSS <= right_ol + antisense_cut3 and TSS > right_ol + antisense_cut2):
                        class_list.append('aTSS-3dn')
                        class_loci[id]['aTSS-3dn'].append(name)
                    if (TSS >= left - antisense_cut2 and TSS < left - antisense_cut1) or (TSS >= left_ol - antisense_cut2 and TSS < left_ol - antisense_cut1):
                        class_list.append('aTSS-2up')
                        class_loci[id]['aTSS-3up'].append(name)
                    if (TSS >= left - antisense_cut3 and TSS < left - antisense_cut2) or (TSS >= left_ol - antisense_cut3 and TSS < left_ol - antisense_cut2):
                        class_list.append('aTSS-3up')
                        class_loci[id]['aTSS-3up'].append(name)
        classes.append(class_list)
        continue
    else:
        for ann_ind, ann_row in ann_list[ann_start:].iterrows():
            left = int(ann_row['left'])
            right = int(ann_row['right'])
            name = ann_row['name']
            if left - TSS > 2*max(primary_cut,antisense_cut3):
                break
            if TSS - right > 2*max(primary_cut,antisense_cut3):
                ann_start = ann_ind
            
            if strand == '+':
                if ann_row['strand'] == '+' and left - TSS >= inside_cut and left - TSS <= primary_cut:
                    if name not in feat_TSS:
                        feat_TSS[name] = [TSS]
                    else:
                        feat_TSS[name].append(TSS)
                if ann_row['strand'] == '+' and left - TSS < inside_cut and right >= TSS:
                    class_list.append('internal')
                    class_loci[id]['internal'].append(name)
                if ann_row['strand'] == '-':
                    if (TSS >= left - antisense_cut1 and TSS <= right + antisense_cut1):
                        class_list.append('aTSS-1')
                        class_loci[id]['aTSS-1'].append(name)
                    if (TSS <= right + antisense_cut2 and TSS > right + antisense_cut1):
                        class_list.append('aTSS-2up')
                        class_loci[id]['aTSS-2up'].append(name)
                    if (TSS <= right + antisense_cut3 and TSS > right + antisense_cut2):
                        class_list.append('aTSS-3up')
                        class_loci[id]['aTSS-3up'].append(name)
                    if (TSS >= left - antisense_cut2 and TSS < left - antisense_cut1):
                        class_list.append('aTSS-2dn')
                        class_loci[id]['aTSS-2dn'].append(name)
                    if (TSS >= left - antisense_cut3 and TSS < left - antisense_cut2):
                        class_list.append('aTSS-3dn')
                        class_loci[id]['aTSS-3dn'].append(name)
            if strand == '-':
                if ann_row['strand'] == '-' and TSS - right >= inside_cut and TSS - right <= primary_cut:
                    if name not in feat_TSS:
                        feat_TSS[name] = [TSS]
                    else:
                        feat_TSS[name].append(TSS)
                if ann_row['strand'] == '-' and TSS - right < inside_cut and TSS >= left:
                    class_list.append('internal')
                    class_loci[id]['internal'].append(name)
                if ann_row['strand'] == '+':
                    if (TSS >= left - antisense_cut1 and TSS <= right + antisense_cut1):
                        class_list.append('aTSS-1')
                        class_loci[id]['aTSS-1'].append(name)
                    if (TSS <= right + antisense_cut2 and TSS > right + antisense_cut1):
                        class_list.append('aTSS-2dn')
                        class_loci[id]['aTSS-2dn'].append(name)
                    if (TSS <= right + antisense_cut3 and TSS > right + antisense_cut2):
                        class_list.append('aTSS-3dn')
                        class_loci[id]['aTSS-3dn'].append(name)
                    if (TSS >= left - antisense_cut2 and TSS < left - antisense_cut1):
                        class_list.append('aTSS-2up')
                        class_loci[id]['aTSS-2up'].append(name)
                    if (TSS >= left - antisense_cut3 and TSS < left - antisense_cut2):
                        class_list.append('aTSS-3up')
                        class_loci[id]['aTSS-3up'].append(name)
        classes.append(class_list)

all_TSS_list['class'] = classes

primary_list = []
secondary_list = []
for feat in feat_TSS:
    cands = feat_TSS[feat]
    values = []
    if len(cands) == 1:
        primary_list.append(int(cands[0]))
        class_loci[str(cands[0])]['primary'].append(feat)
    if len(cands) > 1:
        for cand in cands:
            value = abs(all_TSS_list['value'][int(cand)])
            values.append(value)
        ind_max = np.argwhere(values == np.amax(values)).flatten().tolist()
        for n in range(len(cands)):
            if n in ind_max:
                primary_list.append(int(cands[n]))
                class_loci[str(cands[n])]['primary'].append(feat)
            else:
                secondary_list.append(int(cands[n]))
                class_loci[str(cands[n])]['secondary'].append(feat)
                
primary_list = list(set(primary_list))
secondary_list = list(set(secondary_list))
for TSS in primary_list:
    temp = all_TSS_list.loc[TSS,'class']
    temp.append('primary')
    all_TSS_list.at[TSS,'class'] = temp
for TSS in secondary_list:
    temp = all_TSS_list.loc[TSS,'class']
    temp.append('secondary')
    all_TSS_list.at[TSS,'class'] = temp

for TSS in all_TSS_list['nt']:
    if all_TSS_list.loc[TSS,'class'] == []:
        all_TSS_list.at[TSS,'class'] = ['orphan']

fout_class = open('TSS_classification_{}.txt'.format(date),'w')
firstline = '\t'.join(['conf','nt','strand','primary','secondary','internal','aTSS-1','aTSS-2up','aTSS-3up','aTSS-2dn','aTSS-3dn','orphan'])+'\n'
fout_class.write(firstline)
for index, row in all_TSS_list.iterrows():
    out_items = []
    out_items.append(row['conf'])
    out_items.append(str(int(row['nt'])))
    strand = row['strand']
    out_items.append(strand)
    for classif in pos_classes:
        if classif in row['class']:
            out_items.append('1')
        else:
            out_items.append('0')
    newline = '\t'.join(out_items)+'\n'
    fout_class.write(newline)
fout_class.close()