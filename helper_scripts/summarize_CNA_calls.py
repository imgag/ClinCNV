import sys
import os
from collections import defaultdict

DEFAULT_FDR_THRESHOLD = 0.05
DEFAULT_QC_FAILED_THRESHOLD = 1.0


def parse_CNAs(file):
    with open(file) as f:
        f.readline()
        f.readline()
        f.readline()
        fdr = (f.readline().split(":")[-1]).strip()
        f.readline()
        ploidy = (f.readline().split(":")[-1]).strip()
        clonality = (f.readline().split(":")[-1]).strip()
    return fdr, ploidy, clonality


def clean_file(file, output_file, fdr_threshold, sample_name):
    neutral = []
    lines = []
    output_of_sample = True
    homozygous_deletion_recall = []
    with open(file) as f:
        for i in range(8):
            lines.append(f.readline())

        fdr = (lines[3].split(":")[-1]).strip()
        print(fdr)

        if fdr == "NA":
            while True:
                line = f.readline()
                if not line:
                    break
                lines.append(line)
        else:
            while True:
                line = f.readline()
                if not line:
                    break
                splitted_line = line.split("\t")
                if int(splitted_line[5]) == 0 and not splitted_line[0].strip() == "chrY":
                    length_of_cnv = int(splitted_line[2]) - int(splitted_line[1])
                    if length_of_cnv > 10**7:
                        homozygous_deletion_recall.append(
                            [length_of_cnv, splitted_line[0], splitted_line[1], splitted_line[2]]
                        )
                if float(fdr) < fdr_threshold:
                    lines.append(line)
                else:
                    if len(splitted_line) > 1:
                        allele_balance = splitted_line[5] == splitted_line[6]
                        secondary_allele_balance = splitted_line[12] == splitted_line[13]
                        if not allele_balance or not secondary_allele_balance:
                            if not splitted_line[-2].strip() == "NA":
                                BAF_qval = float(splitted_line[-2])
                                if BAF_qval < 0.05:
                                    lines.append(line)
                                else:
                                    splitted_line = line.split("\t")
                                    neutral.append("\t".join([
                                        splitted_line[0], splitted_line[1], splitted_line[2], "2", splitted_line[11]
                                    ]))
                            else:
                                lines.append(line)
                        else:
                            lines.append(line)

    if len(homozygous_deletion_recall) > 0:
        longest_cnvs = sorted(homozygous_deletion_recall, reverse=True)[0]
        print(sorted(homozygous_deletion_recall))
        print("Length of the largest homozygous variant: ", longest_cnvs[0] / 10**6, " MBs")
        print("You may want to recall your samples with")
        print("--guideBaseline {}:{}-{} --reanalyseCohort --tumorSample {} --normalSample {}".format(
            longest_cnvs[1], longest_cnvs[2], longest_cnvs[3], sample_name.split("-")[0], sample_name.split("-")[1]
        ))
    if output_of_sample:
        with open(output_file, "w") as f:
            for line in lines:
                f.write(line)
    return neutral


def process_directory(in_directory, out_directory, fdr_threshold, qc_failed_threshold):
    path = in_directory
    names_list = []
    fdr_list = []
    ploidy_list = []
    clonality_list = []
    neutral_lines = defaultdict(list)
    # HERE YOU CAN PUT A LIST OF STRINGS WITH SAMPLE THAT DID NOT PASS YOUR QC FOR OTHER REASONS
    list_of_qc_failed_samples = []
    header = ""
    for r, d, f in os.walk(path):
        for file in f:
            if file.startswith("CNAs_"):
                sample_name = file[5:-4]
                if not sample_name in list_of_qc_failed_samples:
                    fdr, ploidy, clonality = parse_CNAs(r + "/" + file)
                    if not fdr == "NA":
                        if float(fdr) > qc_failed_threshold:
                            break
                    names_list.append(sample_name)
                    print(sample_name)
                    fdr_list.append(fdr)
                    ploidy_list.append(ploidy)
                    clonality_list.append(clonality)
                    neutral_regions = clean_file(r + "/" + file, out_directory + file, fdr_threshold, sample_name)
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
            f.write("{}\t{}\t{}\t{}\n".format(sample_n, fdr_list[i], ploidy_list[i], clonality_list[i]))


def main() -> None:
    in_directory = sys.argv[1]
    out_directory = sys.argv[2]
    fdr_threshold = float(sys.argv[3]) if len(sys.argv) > 3 else DEFAULT_FDR_THRESHOLD
    qc_failed_threshold = float(sys.argv[4]) if len(sys.argv) > 4 else DEFAULT_QC_FAILED_THRESHOLD
    print("FDR threshold:", fdr_threshold)
    process_directory(in_directory, out_directory, fdr_threshold, qc_failed_threshold)


if __name__ == "__main__":
    main()
