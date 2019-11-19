import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, NamedTuple, Optional, IO

# HERE YOU CAN PUT A LIST OF STRINGS WITH SAMPLE THAT DID NOT PASS YOUR QC FOR OTHER REASONS
QC_FAILED_SAMPLES = []

DEFAULT_FDR_THRESHOLD = 0.05
DEFAULT_QC_FAILED_THRESHOLD = 1.0

ABNORMALITIES_FILE_PREFIX = 'CNAs_'
NEUTRAL_FILE_PREFIX = 'CNneutral_'

SAMPLE_NAME_REGEX = re.compile(r'(?:{}|{})(.+)'.format(ABNORMALITIES_FILE_PREFIX, NEUTRAL_FILE_PREFIX))

CNA_HEADER_SIZE = 8
CNA_HEADER_FDR_LINE = 3
CNA_HEADER_PLOIDY_LINE = 5
CNA_HEADER_CLONALITY_LINE = 6

CNA_CHROMOSOME_FIELD = 0
CNA_START_FIELD = 1
CNA_END_FIELD = 2
CNA_MAJOR_CN_ALLELE_FIELD = 5
CNA_MINOR_CN_ALLELE_FIELD = 6
CNA_NUMBER_OF_REGIONS_FIELD = 11
CNA_MAJOR_CN_ALLELE2_FIELD = 12
CNA_MINOR_CN_ALLELE2_FIELD = 13
CNA_BAF_QVAL_FDR_FIELD = 22


CNAHeaderInfo = NamedTuple("CNAHeaderInfo", [
    ("header", List[str]),
    ("fdr", str),
    ("ploidy", str),
    ("clonality", str),
])


Sample = NamedTuple("Sample", [
    ("name", str),
    ("info", CNAHeaderInfo),
])


def get_sample_name(path: Path) -> Optional[str]:
    return SAMPLE_NAME_REGEX.match(path.stem).group(1)


def get_header_value(line: str) -> str:
    return line.rpartition(":")[-1].strip()


def parse_CNAs_header(file: IO[str]) -> CNAHeaderInfo:
    header = [next(file) for _ in range(CNA_HEADER_SIZE)]
    fdr = get_header_value(header[CNA_HEADER_FDR_LINE])
    ploidy = get_header_value(header[CNA_HEADER_PLOIDY_LINE])
    clonality = get_header_value(header[CNA_HEADER_CLONALITY_LINE])
    return CNAHeaderInfo(header, fdr, ploidy, clonality)


def clean_file(file: IO[str], output_path: Path, fdr_threshold: float, sample: Sample) -> List[str]:
    neutral = []
    lines = list(sample.info.header)
    output_of_sample = True
    homozygous_deletion_recall = []

    print(sample.info.fdr)
    if sample.info.fdr == "NA":
        lines.extend(file.readlines())
    else:
        for line in file:
            fields = line.split("\t")
            start = int(fields[CNA_START_FIELD])
            end = int(fields[CNA_END_FIELD])
            length_of_cnv = end - start
            chromosome = fields[CNA_CHROMOSOME_FIELD].strip()
            major_cn_allele = fields[CNA_MAJOR_CN_ALLELE_FIELD].strip()
            minor_cn_allele = fields[CNA_MINOR_CN_ALLELE_FIELD].strip()
            major_cn_allele2 = fields[CNA_MAJOR_CN_ALLELE2_FIELD].strip()
            minor_cn_allele2 = fields[CNA_MINOR_CN_ALLELE2_FIELD].strip()
            BAF_qval_fdr = fields[CNA_BAF_QVAL_FDR_FIELD].strip()

            if major_cn_allele == "0" and chromosome != "chrY" and length_of_cnv > 10**7:
                homozygous_deletion_recall.append([length_of_cnv, chromosome, start, end])

            if float(sample.info.fdr) < fdr_threshold:
                lines.append(line)
                continue

            if len(fields) == 1:
                continue

            allele_balance = major_cn_allele == minor_cn_allele
            secondary_allele_balance = major_cn_allele2 == minor_cn_allele2
            if allele_balance and secondary_allele_balance or BAF_qval_fdr == "NA" or float(BAF_qval_fdr) < 0.05:
                lines.append(line)
                continue

            neutral.append("\t".join([chromosome, str(start), str(end), "2", fields[CNA_NUMBER_OF_REGIONS_FIELD]]))

    if homozygous_deletion_recall:
        longest_cnvs = max(homozygous_deletion_recall)
        length, chromosome, start, end = longest_cnvs
        print(sorted(homozygous_deletion_recall))
        print("Length of the largest homozygous variant:", length / 10**6, "MBs")
        print("You may want to recall your samples with")
        print("--guideBaseline {}:{}-{} --reanalyseCohort --tumorSample {} --normalSample {}".format(
            chromosome, start, end, *sample.name.split("-")
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

        with file_path.open() as file:
            if file_path.name.startswith(ABNORMALITIES_FILE_PREFIX):
                header_info = parse_CNAs_header(file)
                if header_info.fdr != "NA" and float(header_info.fdr) > qc_failed_threshold:
                    continue
                sample = Sample(sample_name, header_info)
                samples.append(sample)
                print(sample_name)

                neutral_regions = clean_file(file, out_directory / file_path.name, fdr_threshold, sample)
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
