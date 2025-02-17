__author__ = "gdemidov"

import gzip
import sys
from pathlib import Path
from typing import Callable, Dict, IO, List, Optional

CHROMOSOME_COLUMN = 0
POSITION_COLUMN = 1
REFERENCE_BASES_COLUMN = 3
ALLELIC_BASES_COLUMN = 4
QUALITY_COLUMN = 5
SAMPLE_FIELDS_COLUMN = 8
SAMPLES_START_COLUMN = 9

GENOTYPE_FIELD = "GT"
DEPTH_FIELD = "DP"
OTHER_DEPTH_FIELD = "AO"
ALLELE_FREQUENCY_FIELD = "AF"


def parse_sample_data(fields: List[str], sample: str) -> Dict[str, str]:
    return dict(zip(fields, sample.strip().split(":")))


def parse_vcf_line(line: str, minimum_quality: float) -> Optional[List]:
    if line.startswith("#"):
        return None

    vcf_fields = line.split("\t")
    chromosome = vcf_fields[CHROMOSOME_COLUMN]
    start = vcf_fields[POSITION_COLUMN]
    end = vcf_fields[POSITION_COLUMN]  # is this right?

    reference_variant = vcf_fields[REFERENCE_BASES_COLUMN]
    if len(reference_variant) > 1:
        return None

    allelic_bases = vcf_fields[ALLELIC_BASES_COLUMN].split(",")
    if len(allelic_bases) > 2 or len(allelic_bases[0]) > 1:  # why compare allelic_bases with 2?
        return None
    allelic_variant = allelic_bases[0]

    quality = float(vcf_fields[QUALITY_COLUMN])
    if quality < minimum_quality:
        return None

    sample_fields = vcf_fields[SAMPLE_FIELDS_COLUMN].split(":")

    samples = vcf_fields[SAMPLES_START_COLUMN:]
    baf_depth = parse_sample_data(sample_fields, samples[0])
    if baf_depth[GENOTYPE_FIELD] != "0/1":
        return None

    depth = int(baf_depth[DEPTH_FIELD])
    if depth <= 1:
        return None

    if ALLELE_FREQUENCY_FIELD in baf_depth:
        frequency = baf_depth[ALLELE_FREQUENCY_FIELD]
    elif OTHER_DEPTH_FIELD in baf_depth:
        frequency = float(baf_depth[OTHER_DEPTH_FIELD]) / depth
    else:
        sys.exit("No AO or AF fields in VCF format column. Quit.")

    return [
        chromosome, start, end, reference_variant, allelic_variant, chromosome + "_" + start, frequency, depth
    ]


def write_filtered_lines(input_file: IO[str], output_file: IO[str], minimum_quality: float) -> None:
    for line in input_file:
        parsed_line = parse_vcf_line(line, minimum_quality)
        if parsed_line:
            output_file.write("\t".join(map(str, parsed_line)) + "\n")


def write_header(output_file: IO[str]) -> None:
    output_file.write("#chr\tstart\tend\tref\tobs\tid\tfreq\tdepth\n")


def get_file_reader(input_path: Path) -> Callable[[str], IO[str]]:
    suffix = input_path.suffix.lower()
    if suffix == ".gz":
        return gzip.open
    if suffix == ".vcf":
        return open
    raise Exception('Unknown file format', suffix)


def parse_vcf(input_path: Path, output_path: Path, minimum_quality: float) -> None:
    file_reader = get_file_reader(input_path)
    with file_reader(str(input_path)) as input_file, open(str(output_path), "w") as output_file:
        write_header(output_file)
        write_filtered_lines(input_file, output_file, minimum_quality)


def main() -> None:
    input_file_name = Path(sys.argv[1])
    output_file_name = Path(sys.argv[2])
    minimum_quality = float(sys.argv[3])
    parse_vcf(input_file_name, output_file_name, minimum_quality)


if __name__ == '__main__':
    main()
