from itertools import groupby
import numpy as np
import pandas as pd
from pathlib import Path
import os

# Path to base folder from this code
git_path = str(Path(os.getcwd()).parents[1])
print(git_path)

# Define parser of ATtRACT position probability matrices (PPM's)
def iter_fasta(filename):
    with open(filename, "r") as f:
        iterator = groupby(f, lambda line: line[0] == ">")
        for is_header, group in iterator:
            if is_header:
                header = list(group)[0].split("\t")[0].split(">")[1]
            else:
                group_array = []
                for line in list(group):
                    line_array = []
                    for value in line.split("\t"):
                        # Convert to "pseudo" counts
                        counts = int(float(value)*1e4)
                        line_array.append(counts)
                    group_array.append(line_array)
                yield header, np.array(group_array).T


# Load info about ATtRACT database
attract_info = pd.read_table(git_path + "/data/ATtRACT_db.txt")


for header, matrix in iter_fasta(git_path + "/data/pwm.txt"):
    # Extract gene name and id for a given matrix id
    row_filter = attract_info['Matrix_id'] == header
    subset = attract_info[row_filter]
    subset = subset.dropna()
    genes = [x+"_"+y for x,y in zip(subset.Gene_name, subset.Gene_id) if "ENSG" in y]
    genes = list(set(genes))

    for gene in genes:
        gene_name = gene.split("_")[0]

        # Find top motifs. At least 0.9 of max Score for each RBP.
        thrld = 0.9
        subset = attract_info[attract_info.Gene_name == gene_name]
        subset = subset.dropna()
        subset.loc[:, 'Score'] = subset.Score.str.replace("[*]", "")
        best_motifs = subset.Matrix_id[subset.Score.astype('float') > thrld * max(subset.Score.astype('float'))].values

        if header in best_motifs:
            # There are multiple motifs for each gene
            path = git_path + "/data"
            filename = "{0}/{1}_{2}.pfm".format(path, gene, header)

            print("working with gene {0} ...".format(gene_name))

            # Save nucleotides counts in PFM format
            # Keep in mind that ATtRACT follows A C G T order
            # While MOODS expects A C G T order
            # We should switch row 3 and row 1 in a ATtRACT matrix (indices start from 0)
            # to match MOODS order
            with open(filename, "w") as f:
                text = ""
                for ind in [0, 1, 2, 3]:
                    row = matrix[ind]
                    for value in row:
                        text += str(value)
                        text += " "
                    # Remove space at the end of a line
                    text = text[:-1]
                    text += "\n"
                f.write(text[:-1])
    # Run the code only for 1 RBP
    break
