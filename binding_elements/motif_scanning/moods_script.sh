#!/bin/bash

# Revert sequences on "-" strand

# Run code from MOODS{version info}/scripts directory
# We select BE only for human RBPs
# -m path to PFM files
# -s path sequences for scanning
# -p significance threshold (how does it work?)
# -R disable matching versus reverse complement strand

# Consult MOODS webpage for a proper installation:
# https://github.com/jhkorhonen/MOODS/wiki/Installation
# For example, install option 1 should do the job Dec, 2019 suggestion.


# Define files
moods=../../external_scripts/MOODS-python-1.9.4/scripts/moods-dna.py
moods=../../external_scripts/MOODS-python-1.9.4/scripts/moods-dna.py
seqs=../../data/qapa_hg19_3utrs.fa
output=../../binding_elements/binding_elements_3utrs.tab

# Processing ATtRACT RBP motif. Example RBP.
python $moods -m ../../data/HNRNPH2_ENSG00000126945_904.pfm -s $seqs -p 0.0001 -R > $output
