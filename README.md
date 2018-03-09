# Spneu-deep-annotation
Analysis scripts used for sRNA, TSS and terminator detection in Streptococcus pneumoniae D39V. Results are presented in [PneumoBrowse](https://veeninglab.com/pneumobrowse)

In summary, the pipeline consisted of the following steps:
1. add_overlap.py; extend reference genome to allow mapping over FASTA boundaries
2. bowtie2-build ([Manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) | [REF](https://www.ncbi.nlm.nih.gov/pubmed/22388286))
3. bowtie2 ([Manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) | [REF](https://www.ncbi.nlm.nih.gov/pubmed/22388286)); mapping FASTQ reads to extended genome ([REF](https://www.ncbi.nlm.nih.gov/pubmed/22388286))
4. process_annotation.py; extract coordinates of annotated features, etc.
5. process_sam_nomulti.py (iter=1); extract mapping parameters (NO multimappers), including fragment sizes of paired-end sequenced data
6. TSS_analysis.py - as support in sRNA analysis, repeated later on data including multimappers
7. TransTermHP ([Website](http://transterm.cbcb.umd.edu/) | [REF](https://www.ncbi.nlm.nih.gov/pubmed/17313685)); as support in sRNA analysis, repeated later
8. terminator_analysis.py; as support in sRNA analysis, repeated later on data including multimappers
9. small_frags.py (iter=1); find overrepresented small transcripts in paired-end data
10. process_sam_nomulti.py (iter=2); remove fragment sizes of just-predicted sRNAs from library-wide distribution
11. small_frags.py (iter=2); repeat analysis with improved library-wide distribution from step 10.
12. process_sam.py; extract mapping parameters (INCLUDING multimappers)
13. Repeat step 6
14. Repeat step 7
15. Repeat step 8
16. TSS_classification.py; classify TSSs as Primary, Secondary, Internal, Antisense or Orphan
17. operons.py; predict operons based on strong terminators (step 15) and expression correlation data from [PneumoExpress](https://veeninglab.com/pneumoexpress).
