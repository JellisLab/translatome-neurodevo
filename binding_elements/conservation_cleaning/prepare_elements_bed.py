import sys
import pandas as pd

# Input binding elements from shell
input_file = sys.argv[1]
input_type = input_file.split(".")[-1]
input_len = len(input_file.split("/"))
assert input_type == "tab", "input file for found binding elements should be in tab format"
assert input_len != 1, "full path to an input file is required"
input_base = input_file.split(".")[0]

# Coordinates of QAPA 3`UTRs
inpath = "/path/to/qapa_hg19_3utrs.bed"
qapa_annotation = pd.read_csv(inpath, sep='\t', header = None)

bedfile = open(input_base + ".bed", "w")
with open(input_file, "r") as f:
    i = 0
    for l in f:
        # Sequence ID of 3`UTR isoform
        l_split = l.split(",")
        sequence_id = ','.join(l_split[:-6])
        motif_id = ','.join(l_split[-6:])

        # Strand
        strand = sequence_id.split(",")[-1].split("_")[6]
        if i % 1000 == 0:
            print(i)

        # 3`UTR information: chromosome, coordinates
        selector = qapa_annotation.iloc[:,3] == sequence_id
        utr_info = qapa_annotation[selector].values[0]
        chrom = utr_info[0]
        utr_start = int(utr_info[1])
        utr_end = int(utr_info[2])

        # RBP element information: start, width
        element_start = int(motif_id.split(",")[1])
        element_width = len(motif_id.split(",")[4])

        # Write coordinates of a binding element. Keep in mind
        # bed files are in 0-based format.
        if strand == "+":
            coordinate_A = utr_start + element_start
            coordinate_B = utr_start + element_start + element_width
            bedfile.write("{0} {1} {2} {3}".format(chrom, coordinate_A, coordinate_B, l))
        else:
            coordinate_A = utr_end - element_start - element_width
            coordinate_B = utr_end - element_start
            bedfile.write("{0} {1} {2} {3}".format(chrom, coordinate_A, coordinate_B, l))

        # Update index
        i += 1

bedfile.close()
