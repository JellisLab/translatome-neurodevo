import pandas as pd
import numpy as np
import time

def count_elements_in_genes(moods_output):
    """
    MOODS (python code) scan sequences with PFM (probabilistic representation of
    a binding motif)
    """

    # Note: Gene PWWP2A appeared from a match to pattern WWP2
    print "populating design matrix"

    with open(moods_output) as f:
        i = 0
        # Keep track of found RBPs
        found_rbps = []

        # Dicts to store design matrix
        binary_design_matrix = {}
        counts_design_matrix = {}
        st = time.time()
        for line in f:
            # Parse line
            tx_name = line.split(",")[0].split("_")[1]
            rbp_name = line.split("ENSG0")[1].split("_")[1].split(".")[0]

            # Add dict if RBP is new
            if not (rbp_name in found_rbps):
                found_rbps.append(rbp_name)
                binary_design_matrix[rbp_name] = {}
                counts_design_matrix[rbp_name] = {}

            # Store BE
            binary_design_matrix[rbp_name][tx_name] = 1
            try:
                counts_design_matrix[rbp_name][tx_name] += 1
            except:
                counts_design_matrix[rbp_name][tx_name] = 1

            # Update counter
            i += 1
            if i % 500000 == 0:
                print i
        print "elapsed time = ", time.time() - st

    return binary_design_matrix, counts_design_matrix

def count_elements_in_3UTRs(moods_output, conservation_threshold):
    """
    MOODS (python code) scan sequences with PFM (probabilistic representation of
    a binding motif)

    Count binding motifs in each 3'UTR isoforms
    """

    # Note: Gene PWWP2A appeared from a match to pattern WWP2
    print "populating design matrix"

    print "processing file", moods_output, "..."

    with open(moods_output, "r") as f:
        i = 0
        # Keep track of found RBPs
        found_rbps = []

        # Dicts to store design matrix
        binary_design_matrix = {}
        counts_design_matrix = {}
        st = time.time()
        for line in f:

            # Parse line
            utr3_name = line.split(",ENSG0")[0]
            rbp_name = line.split("ENSG0")[1].split("_")[1].split(".")[0]
            conservation_score = float(line.split(",")[-1])

            # Add dict if RBP is new
            if not (rbp_name in found_rbps):
                found_rbps.append(rbp_name)
                binary_design_matrix[rbp_name] = {}
                counts_design_matrix[rbp_name] = {}

            # Count BE if element is sufficiently conserved
            if conservation_score > conservation_threshold:
                try:
                    binary_design_matrix[rbp_name][utr3_name] = 1
                    counts_design_matrix[rbp_name][utr3_name] += 1
                except:
                    counts_design_matrix[rbp_name][utr3_name] = 1
            else:
                try:
                    binary_design_matrix[rbp_name][utr3_name] += 0
                    counts_design_matrix[rbp_name][utr3_name] += 0
                except:
                    binary_design_matrix[rbp_name][utr3_name] = 0
                    counts_design_matrix[rbp_name][utr3_name] = 0

            # Update counter
            i += 1
            if i % 1000000 == 0:
                print i

        print "elapsed time = ", time.time() - st

    return binary_design_matrix, counts_design_matrix

# Count occurences of binding elements in each 3UTR isoform
conservation_threshold = 0.9
binary_dict, counts_dict = count_elements_in_3UTRs("path/binding_elements_with_conservation_scores.tab",
                                                   conservation_threshold)

# Convert dictionary to DataFrame
print "converting dict to dataframe"
st = time.time()
binary_df = pd.DataFrame(binary_dict).fillna(0).astype(np.int)
counts_df = pd.DataFrame(counts_dict).fillna(0).astype(np.int)
print "elapsed time = ", time.time() - st

# Save dataframes
print "saving dataframes"
st = time.time()
binary_df.to_csv("path/binary_filtered_design.csv")
counts_df.to_csv("path/counts_filtered_design.csv")
print "elapsed time = ", time.time() - st
