import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, NamedTuple, Optional

# HERE YOU CAN PUT A LIST OF STRINGS WITH SAMPLE THAT DID NOT PASS YOUR QC FOR OTHER REASONS
QC_FAILED_SAMPLES = []

DEFAULT_FDR_THRESHOLD = 0.05
DEFAULT_QC_FAILED_THRESHOLD = 1.0

ABNORMALITIES_FILE_PREFIX = 'CNAs_'
NEUTRAL_FILE_PREFIX = 'CNneutral_'

SAMPLE_NAME_REGEX = re.compile(r'(?:{}|{})(.+)'.format(ABNORMALITIES_FILE_PREFIX, NEUTRAL_FILE_PREFIX))

SampleInfo = NamedTuple("SampleInfo", [
    ("fdr", str),
    ("ploidy", str),
    ("clonality", str),
])


Sample = NamedTuple("Sample", [
    ("name", str),
    ("info", SampleInfo),
])


def get_sample_name(path: Path) -> Optional[str]:
    return SAMPLE_NAME_REGEX.match(path.stem).group(1)


def parse_CNAs(path: Path) -> SampleInfo:
    with path.open() as f:
        f.readline()
        f.readline()
        f.readline()
        fdr = (f.readline().split(":")[-1]).strip()
        f.readline()
        ploidy = (f.readline().split(":")[-1]).strip()
        clonality = (f.readline().split(":")[-1]).strip()
    return SampleInfo(fdr, ploidy, clonality)


def clean_file(path: Path, output_path: Path, fdr_threshold: float, sample_name: str) -> List[str]:
    neutral = []
    lines = []
    output_of_sample = True
    homozygous_deletion_recall = []
    with path.open() as input_file:
        for i in range(8):
            lines.append(input_file.readline())

        fdr = lines[3].rpartition(":")[-1].strip()
        print(fdr)

        if fdr == "NA":
            lines.extend(input_file.readlines())
        else:
            for line in input_file:
                fields = line.split("\t")
                length_of_cnv = int(fields[2]) - int(fields[1])
                if int(fields[5]) == 0 and fields[0].strip() != "chrY" and length_of_cnv > 10**7:
                    homozygous_deletion_recall.append([length_of_cnv, fields[0], fields[1], fields[2]])

                if float(fdr) < fdr_threshold:
                    lines.append(line)
                    continue

                if len(fields) == 1:
                    continue

                allele_balance = fields[5] == fields[6]
                secondary_allele_balance = fields[12] == fields[13]
                if allele_balance and secondary_allele_balance or fields[-2].strip() == "NA" or float(fields[-2]) < 0.05:
                    lines.append(line)
                    continue

                fields = line.split("\t")
                neutral.append("\t".join([fields[0], fields[1], fields[2], "2", fields[11]]))

    if homozygous_deletion_recall:
        longest_cnvs = max(homozygous_deletion_recall)
        print(sorted(homozygous_deletion_recall))
        print("Length of the largest homozygous variant:", longest_cnvs[0] / 10**6, "MBs")
        print("You may want to recall your samples with")
        print("--guideBaseline {}:{}-{} --reanalyseCohort --tumorSample {} --normalSample {}".format(
            longest_cnvs[1], longest_cnvs[2], longest_cnvs[3], sample_name.split("-")[0], sample_name.split("-")[1]
        ))
    if output_of_sample:
        with output_path.open("w") as output_file:
            for line in lines:
                output_file.write(line)
    return neutral


def write_neutral_file(output_path: Path, header: str, lines: List[str]) -> None:
    with output_path.open("w") as output:
        print(header, file=output)
        for line in lines:
            print(line, file=output)


def write_neutral_files(directory: Path, neutral_lines: Dict[str, List[str]], header: str) -> None:
    for sample_name, lines in neutral_lines.items():
        output_path = directory / ("neutral_" + sample_name + ".txt")
        write_neutral_file(output_path, header, lines)


def write_summarized_for_fdr(path: Path, samples: List[Sample]):
    with path.open("w") as f:
        for s in samples:
            print(s.name, s.info.fdr, s.info.ploidy, s.info.clonality, sep='\t', file=f)


def get_all_files(directory: Path) -> Iterable[Path]:
    for dir_path, _, file_names in os.walk(str(directory)):
        for file_name in file_names:
            yield Path(dir_path) / file_name


def process_directory(
        in_directory: Path, out_directory: Path, fdr_threshold: float, qc_failed_threshold: float) -> None:
    samples = []
    neutral_lines = defaultdict(list)
    header = ""
    for file_path in get_all_files(in_directory):
        sample_name = get_sample_name(file_path)
        if not sample_name or sample_name in QC_FAILED_SAMPLES:
            continue

        if file_path.name.startswith(ABNORMALITIES_FILE_PREFIX):
            sample_info = parse_CNAs(file_path)
            if sample_info.fdr != "NA" and float(sample_info.fdr) > qc_failed_threshold:
                continue
            samples.append(Sample(sample_name, sample_info))
            print(sample_name)

            neutral_regions = clean_file(file_path, out_directory / file_path.name, fdr_threshold, sample_name)
            neutral_lines[sample_name].extend(neutral_regions)

        if file_path.name.startswith(NEUTRAL_FILE_PREFIX):
            with file_path.open() as neutral_file:
                header = neutral_file.readline().strip()
                for line in neutral_file:
                    neutral_lines[sample_name].append(line.strip())

    write_neutral_files(out_directory, neutral_lines, header)
    write_summarized_for_fdr(out_directory / "summarized_for_FDR.txt", samples)


def main() -> None:
    in_directory = Path(sys.argv[1])
    out_directory = Path(sys.argv[2])
    fdr_threshold = float(sys.argv[3]) if len(sys.argv) > 3 else DEFAULT_FDR_THRESHOLD
    qc_failed_threshold = float(sys.argv[4]) if len(sys.argv) > 4 else DEFAULT_QC_FAILED_THRESHOLD
    print("FDR threshold:", fdr_threshold)
    process_directory(in_directory, out_directory, fdr_threshold, qc_failed_threshold)


if __name__ == "__main__":
    main()
