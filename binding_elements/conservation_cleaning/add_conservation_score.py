import sys

input_file = sys.argv[1]
input_type = input_file.split(".")[-1]
input_len = len(input_file.split("/"))
assert input_type == "bed", "input file for conservation scores should be in bed format"
assert input_len != 1, "full path to an input file is required"
input_base = input_file.split(".")[0]

cleaned_file = open(input_base + "_reformatted.bed", "w")
with open(input_file, "r") as f:
    i = 0 
    for l in f: 
        original_data = l.replace(" ", "\t").split("\t")[3]
        binding_info = original_data.split("|")[0]
        conservation_score = original_data.split("|")[1]

        if i % 100000 == 0:
            print(i)
        i += 1 
 
        cleaned_file.write("{0}|{1}".format(binding_info[:-1], conservation_score))

cleaned_file.close()
