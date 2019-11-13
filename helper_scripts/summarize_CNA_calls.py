import sys
import os
from collections import defaultdict

def parse_CNAs(file):
    with open(file) as f:
        f.readline()
        f.readline()
        f.readline()
        fdr = (f.readline().split(":")[-1]).strip()
        f.readline()
        ploidy = (f.readline().split(":")[-1]).strip()
        clonality = (f.readline().split(":")[-1]).strip()
    return(fdr, ploidy, clonality)

def clean_file(file, output_file, fdr_threshold):
    neutral = []
    lines = []
    output_of_sample = True
    with open(file) as f:
        for i in range(8):
            lines.append(f.readline())

        fdr = (lines[3].split(":")[-1]).strip()
        print(fdr)

        if fdr == "NA" :
            while f.readline():
                lines.append(f.readline())
        else:
            while True:
                line = f.readline()
                if not line: break
                if float(fdr) < fdr_threshold:
                    lines.append(line)
                else:
                    splitted_line = line.split("\t")
                    if len(splitted_line) > 1:
                        allele_balance = splitted_line[5] == splitted_line[6]
                        secondary_allele_balance = splitted_line[12] == splitted_line[13]
                        if int(splitted_line[5]) == 0 and not splitted_line[0].strip() == "chrY":
                            length_of_cnv = int(splitted_line[2]) - int(splitted_line[1])
                            if length_of_cnv > 10**7:
                                print("HOMOZYGOUS DELETION - QUITE LONG!!!!")
                        if not allele_balance or not secondary_allele_balance:
                            if not splitted_line[-2].strip() == "NA":
                                BAF_qval = float(splitted_line[-2])
                                if BAF_qval < 0.05:
                                    lines.append(line)
                                else:
                                    splitted_line = line.split("\t")
                                    neutral.append("\t".join([splitted_line[0], splitted_line[1], splitted_line[2], "2",
                                                              splitted_line[11]]))
                            else:
                                lines.append(line)
                        else:
                            lines.append(line)
    if output_of_sample:
        with open(output_file, "w") as f:
            for line in lines:
                f.write(line)
    return(neutral)




def main():
    in_directory = sys.argv[1]
    out_directory = sys.argv[2]
    if len(sys.argv) > 3:
        fdr_threshold = float(sys.argv[3])
    else:
        fdr_threshold = 0.05
    if len(sys.argv) > 4:
        QC_failed_threshold = float(sys.argv[4])
    else:
        QC_failed_threshold = 1.0
    print("FDR threshold: " + str(fdr_threshold))
    path = in_directory
    names_list = []
    fdr_list = []
    ploidy_list = []
    clonality_list = []
    neutral_lines = defaultdict(list)
    list_of_qc_failed_samples = [] # HERE YOU CAN PUT A LIST OF STRINGS WITH SAMPLE THAT DID NOT PASS YOUR QC FOR OTHER REASONS
    header = ""
    for r, d, f in os.walk(path):
        for file in f:
            if file.startswith("CNAs_"):
                sample_name = file[5:-4]
                if not sample_name in list_of_qc_failed_samples:
                    fdr, ploidy, clonality = parse_CNAs(r + "/" + file)
                    if not fdr == "NA":
                        if float(fdr) > QC_failed_threshold:
                            break
                    names_list.append(sample_name)
                    print(sample_name)
                    fdr_list.append(fdr)
                    ploidy_list.append(ploidy)
                    clonality_list.append(clonality)
                    neutral_regions = clean_file(r + "/" + file, out_directory + file, fdr_threshold)
                    neutral_lines[sample_name].extend(neutral_regions)
            if file.startswith("CNneutral"):
                sample_name = file[10:-4]
                if not sample_name in list_of_qc_failed_samples:
                    with open(r + "/" + file) as f:
                        header = f.readline().strip()
                        for line in f:
                            neutral_lines[sample_name].append(line.strip())

    for key in neutral_lines:
        with open(out_directory + "/neutral_" + key + ".txt", "w") as f:
            f.write(header + "\n")
            for elem in neutral_lines[key]:
                f.write(elem + "\n")


    with open(out_directory + "/" + "summarized_for_FDR.txt", "w") as f:
        for i, sample_n in enumerate(names_list):
            f.write(sample_n + "\t" + str(fdr_list[i]) + "\t" + str(ploidy_list[i]) + "\t" + str(clonality_list[i]) + "\n")




main()
