import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import NamedTuple, List

# HERE YOU CAN PUT A LIST OF STRINGS WITH SAMPLE THAT DID NOT PASS YOUR QC FOR OTHER REASONS
QC_FAILED_SAMPLES = []

DEFAULT_FDR_THRESHOLD = 0.05
DEFAULT_QC_FAILED_THRESHOLD = 1.0


SampleInfo = NamedTuple("SampleInfo", [
    ("fdr", str),
    ("ploidy", str),
    ("clonality", str),
])


Sample = NamedTuple("Sample", [
    ("name", str),
    ("info", SampleInfo),
])


def parse_CNAs(file: str) -> SampleInfo:
    with open(file) as f:
        f.readline()
        f.readline()
        f.readline()
        fdr = (f.readline().split(":")[-1]).strip()
        f.readline()
        ploidy = (f.readline().split(":")[-1]).strip()
        clonality = (f.readline().split(":")[-1]).strip()
    return SampleInfo(fdr, ploidy, clonality)


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


def write_summarized_for_fdr(out_directory: Path, samples: List[Sample]):
    with open(str(out_directory / "summarized_for_FDR.txt"), "w") as f:
        for sample in samples:
            f.write("{s.name}\t{s.info.fdr}\t{s.info.ploidy}\t{s.info.clonality}\n".format(s=sample))


def process_directory(
        in_directory: Path, out_directory: Path, fdr_threshold: float, qc_failed_threshold: float) -> None:
    samples = []
    neutral_lines = defaultdict(list)
    header = ""
    for r, _, f in os.walk(str(in_directory)):
        for file in f:
            if file.startswith("CNAs_"):
                sample_name = file[5:-4]
                if sample_name not in QC_FAILED_SAMPLES:
                    sample_info = parse_CNAs(r + "/" + file)
                    if sample_info.fdr != "NA":
                        if float(sample_info.fdr) > qc_failed_threshold:
                            break
                    samples.append(Sample(sample_name, sample_info))
                    print(sample_name)
                    neutral_regions = clean_file(r + "/" + file, out_directory / file, fdr_threshold, sample_name)
                    neutral_lines[sample_name].extend(neutral_regions)
            if file.startswith("CNneutral"):
                sample_name = file[10:-4]
                if sample_name not in QC_FAILED_SAMPLES:
                    with open(r + "/" + file) as f:
                        header = f.readline().strip()
                        for line in f:
                            neutral_lines[sample_name].append(line.strip())
    for key in neutral_lines:
        with open(str(out_directory / ("neutral_" + key + ".txt")), "w") as f:
            f.write(header + "\n")
            for elem in neutral_lines[key]:
                f.write(elem + "\n")
    write_summarized_for_fdr(out_directory, samples)


def main() -> None:
    in_directory = Path(sys.argv[1])
    out_directory = Path(sys.argv[2])
    fdr_threshold = float(sys.argv[3]) if len(sys.argv) > 3 else DEFAULT_FDR_THRESHOLD
    qc_failed_threshold = float(sys.argv[4]) if len(sys.argv) > 4 else DEFAULT_QC_FAILED_THRESHOLD
    print("FDR threshold:", fdr_threshold)
    process_directory(in_directory, out_directory, fdr_threshold, qc_failed_threshold)


if __name__ == "__main__":
    main()
