import pandas as pd
import numpy as np
import sys

# Input file from shell
input_file = sys.argv[1]
input_base = input_file.split(".")[0]

def count_elements_in_3UTRs(moods_output, conservation_threshold):
    """
    MOODS (python code) scans sequences with PFM (probabilistic representation of
    a binding motif)

    Count binding motifs in each 3'UTR isoforms
    """

    # Note: Gene PWWP2A appeared from a match to pattern WWP2
    print("populating design matrix")
    print("processing file {0} ...".format(moods_output))

    with open(moods_output, "r") as f:
        i = 0

        # Dicts to store design matrix
        counts_design_matrix = {}
        for line in f:
            # Parse line
            # utr3_name = line.split(",M")[0]
            # rbp_name = "M" + line.split(",M")[1].split(".pfm")[0]
            # conservation_score = float(line.split("|")[-1])

            # Sequence ID of 3`UTR isoform
            l_split = line.split(",")
            sequence_id = ','.join(l_split[:-5])
            motif_id = ','.join(l_split[-5:])
            rbp_name = motif_id.split(".pfm")[0]
            conservation_score = float(line.split("|")[-1])

            counts_design_matrix[rbp_name] = counts_design_matrix.get(rbp_name, {})

            # Count BE if element is sufficiently conserved
            if conservation_score > conservation_threshold:
                counts_design_matrix[rbp_name][sequence_id] = counts_design_matrix[rbp_name].get(sequence_id, 0) + 1

            # Update counter
            i += 1
            if i % 1000000 == 0:
                print(i)

    return counts_design_matrix


# Count occurences of binding elements in each 3UTR isoform
conservation_threshold = 0.9
counts_dict = count_elements_in_3UTRs(input_file, conservation_threshold)

# Convert dictionary to DataFrame
counts_df = pd.DataFrame(counts_dict).fillna(0).astype(np.int)

# Save dataframes
counts_df.to_csv(input_base + "_conserv_thres_" + str(conservation_threshold) + ".csv")

