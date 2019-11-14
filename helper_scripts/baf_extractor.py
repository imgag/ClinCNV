__author__ = "gdemidov"

import gzip
import sys
from typing import Dict, List, Optional

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
    end = vcf_fields[POSITION_COLUMN]

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


def parse_vcf(input_file_name: str, output_file_name: str, minimum_quality: float) -> None:
    if input_file_name.endswith(".gz"):
        input_file = gzip.open(input_file_name)
    elif input_file_name.endswith((".vcf", ".VCF")):
        input_file = open(input_file_name)
    header = "#chr\tstart\tend\tref\tobs\tid\tfreq\tdepth\n"
    with open(output_file_name, "w") as output_file:
        output_file.write(header)
        for line in input_file:
            parsed_line = parse_vcf_line(line, minimum_quality)
            if parsed_line:
                output_file.write("\t".join(map(str, parsed_line)) + "\n")


def main() -> None:
    input_file_name = sys.argv[1]
    output_file_name = sys.argv[2]
    minimum_quality = float(sys.argv[3])
    parse_vcf(input_file_name, output_file_name, minimum_quality)


if __name__ == '__main__':
    main()
