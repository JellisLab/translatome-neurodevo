#!/bin/bash -x
#PBS -l nodes=1:ppn=1
#PBS -l gres=localhd:10
#PBS -l mem=10g
#PBS -l vmem=10g
#PBS -l walltime=03:30:00

# Load modules
module load bedops/2.4.3
module load python/3.7.1

# Folder to store temporary files
TMPDIR=/localhd/$PBS_JOBID

# Convert wig to bed
gunzip -c /path/to/data/chr${num}.phastCons46way.wigFix.gz | wig2bed - > /path/to/data/chr${num}.sorted.bed

# Define input. It's an output table from MOODS for RBP motifs from ATtRACT.
# It looks like this:
# ENST00000342066_SAMD11_hg19_chr1_879287_879955_+_utr_879533_879955,AGO1_ENSG00000092847_579.pfm,110,+,7.596969579470904,AATTTTAAA,
# 3'UTR_ID,motif_ID,position_in_3`UTR,strand,matching_score,sequence
all_unsorted_tab="/path/to/data/binding_elements_3utrs_ATTRACT_randomBG.tab"

# Define outputs. Sorted data from a chosen chromosome
unsorted_tab="/path/to/data/binding_elements_3utrs_chr${num}.tab"
unsorted_bed="/path/to/data/binding_elements_3utrs_chr${num}.bed"
sorted_bed="/path/to/data/binding_elements_3utrs_chr${num}.sorted.bed"

# Define outputs. Binding elements with average conservation score
phastConFn="/path/to/data/chr${phastnum}.sorted.bed"
quantified_regions="/path/to/data/regions_with_avg_phastcons_chr${num}.bed"
reformatted_regions="/path/to/data/regions_with_avg_phastcons_chr${num}_reformatted.bed"

# Subset of detected binding elements from a same chromosome
cat $all_unsorted_tab | grep chr${num} > $unsorted_tab


# Find and sort binding elements from a single chromosome
python /path/to/code/prepare_elements_bed.py $unsorted_tab
sort -k1,1 -k2,2n $unsorted_bed > $sorted_bed

# Calculate conservation of binding elements
time bedmap --echo --mean --chrom chr${phastnum} ${sorted_bed} ${phastConFn} > ${quantified_regions}

# Reformat bedmap output
python /path/to/code/add_conservation_score.py $quantified_regions

# Count binding elements in each transcript
python /path/to/code/count_elements.py $reformatted_regions

# Clean directory
rm $quantified_regions $unsorted_bed $sorted_bed
