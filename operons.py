"""
Input: (i) genome sequence (fasta)
Command: python operons.py [fasta file]
Required files:
primary_TSS.txt - tab-delimited incl. header (nt, strand) - contains all primary TSSs as identified with 'TSS_classification.py'
strong_terminator.txt - tab-delimited no header (start, stop, strand) - contains all terminators with an efficiency above a certain threshold (we used 80%)
correlations.txt - tab-delimited incl. header (locus, type, start, end, strand, correlation) - contains correlation values of a gene with the next gene. Use correlation=1 for all features if no data available.
all_coords.txt - tab-delimited incl. header (type, start, stop, strand, old_locus_tag, locus_tag) - all annotated genes
"""

import re
from Bio import SeqIO
import sys

genome = SeqIO.read(sys.argv[1],'fasta')
genome_size = len(genome)

ext_genome = genome[:] + genome[:10000]

corr_cut = 0.75
term_overlap = 5
dist_cut = 600

pTSS_file = open('primary_TSS.txt','r')
pTSS_lines = pTSS_file.readlines()
pTSS_file.close()

term_file = open('strong_terminator.txt','r')
term_lines = term_file.readlines()
term_file.close()

corr_file = open('correlations.txt','r')
corr_lines = corr_file.readlines()
corr_file.close()

allcoords_file = open('all_coords.txt','r')
allcoords_lines = allcoords_file.readlines()
allcoords_file.close()

gene_coords = {}
for i in range(1,len(allcoords_lines)):
        featline = allcoords_lines[i]
        featdata = re.split('\t|\n',featline)
        locus = featdata[5]
        start = int(featdata[1])
        end = int(featdata[2])
        featstrand = featdata[3]
        gene_coords[locus] = [start,end,featstrand]

operons = {}
operon_count = 1
operon_list = []
op_coords = {}

firstbase = 1

for i in range(1,len(corr_lines)):
    line = corr_lines[i]
    data = line.split()
    locus = data[0]
    ftype = data[1]
    start = int(data[2])
    end = int(data[3])
    strand = data[4]
    corr = float(data[5])
    
    terminated = False
    nextterminated = False
    terminators = []
    nextterminators = []
    if i != len(corr_lines) - 1:
        nextline = corr_lines[i+1]
        nextdata = nextline.split()
        nextlocus = nextdata[0]
        nextftype = nextdata[1]
        nextstart = int(nextdata[2])
        nextend = int(nextdata[3])
        nextstrand = nextdata[4]
        
        for termline in term_lines[1:]:
            termdata = termline.split()
            termstart = int(termdata[0])
            termend = int(termdata[1])
            termstrand = termdata[2]
            if strand == '+':
                if termstrand == '+' and termend >= end and termstart <= nextstart:
                    terminated = True
                    terminators.append(termend)
            if nextstrand == '-':
                if termstrand == '-' and termend >= end and termstart <= nextstart:
                    nextterminated = True
                    nextterminators.append(termstart)
            
    op_name = 'operon_'+str(operon_count)
    operons[locus] = op_name
    
    if nextstrand != strand or corr < corr_cut or ftype == 'rRNA' or nextftype == 'rRNA' or i == len(corr_lines)-1 or terminated == True or nextterminated == True:
        if terminated == True:
            end = min(terminators)
        op_coords[op_name] = [firstbase,end,strand]
            
        if nextterminated == True:
            firstbase = max(nextterminators)
        else: 
            firstbase = nextstart
        
        operon_count += 1

op_genes = {}
for gene in operons:
    operon = operons[gene]
    if operon not in op_genes:
        op_genes[operon] = [gene]
    else:
        op_genes[operon].append(gene)

fout = open('primary_operons.txt','w')
feat_dict = {}
covered_operons = {}
full_op_count = 0
complete_operons = []
feat_coverage = []
genes_per_id = {}

gene_count = 0
pOpcount = 0
# Create operons from primary TSSs
for TSSline in pTSS_lines[1:]:
    TSSgenes = []    
    nt = int(TSSline.split()[0])
    TSSstrand = TSSline.split()[1]
    feats = {}
    #find next feature    
    for i in range(1,len(allcoords_lines)):
        featline = allcoords_lines[i]
        featdata = re.split('\t|\n',featline)
        locus = featdata[5]
        start = int(featdata[1])
        end = int(featdata[2])
        featstrand = featdata[3]
        feat_dict[locus] = [start,end]
        if TSSstrand == '+':
            dist = start - nt
            if featstrand == '+':
                if dist <= dist_cut and dist >= 0:
                    feats[locus] = dist
        if TSSstrand == '-':
            dist = nt - end
            if featstrand == '-':
                if dist <= dist_cut and dist >= 0:
                    feats[locus] = dist
    # Over genome boundaries
    if nt == 2046536: # specific to D39V
        closest = 'SPV_0001' # specific to D39V
    else:
        closest = min(feats, key=feats.get)
    operon = operons[closest]
    for gene in op_genes[operon]:
        genestart,geneend = feat_dict[gene]
        if TSSstrand == '+' and genestart >= nt:
            TSSgenes.append(gene)
            feat_coverage.append(gene)
        if TSSstrand == '-' and geneend <= nt:
            TSSgenes.append(gene)
            feat_coverage.append(gene)
    if operon not in covered_operons:
        covered_operons[operon] = 1
    else: 
        covered_operons[operon] += 1
    coords = op_coords[operon]
    if TSSstrand == '+':
        if nt <= coords[0]:
            full_op_count += 1
            complete_operons.append(operon)
            gene_count += len(op_genes[operon])
        coords[0] = nt
    if TSSstrand == '-':
        if nt >= coords[1]:
            full_op_count += 1
            complete_operons.append(operon)
            gene_count += len(op_genes[operon])
        coords[1] = nt
    operon_id = operon+"_"+str(covered_operons[operon])
    pOpcount += 1
    gff_opid = 'operon_'+str(pOpcount)
    genes_per_id[gff_opid] = TSSgenes
    newline = 'D39_JWV.0\tmanual\toperon\t{}\t{}\t.\t{}\t.\tID={};Name={}\n'.format(coords[0],coords[1],coords[2],gff_opid,gff_opid)
    fout.write(newline)
fout.close()