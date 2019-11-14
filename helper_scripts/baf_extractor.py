__author__ = 'gdemidov'

import gzip
import sys
from typing import List, Optional


def vcf_parse_one_line(line: str, minimum_quality: float) -> Optional[List]:
    if line.startswith("#"):
        return None

    sample_data_types = line.split("\t")[8].split(":")
    genotype_index = sample_data_types.index("GT")
    depth_index = sample_data_types.index("DP")
    try:
        other_depth_index = sample_data_types.index("AO")
    except ValueError:
        other_depth_index = -1
    try:
        frequency_index = sample_data_types.index("AF")
    except ValueError:
        frequency_index = -1

    vcf_fields = line.split("\t")
    chromosome = vcf_fields[0]
    start = vcf_fields[1]
    end = vcf_fields[1]
    reference_variant = vcf_fields[3]
    allelic_variant = vcf_fields[4].split(",")
    if len(allelic_variant[0]) > 1 or len(reference_variant) > 1 or len(allelic_variant) > 2:
        return None
    quality = float(vcf_fields[5])
    if quality < minimum_quality:
        return None
    baf_depth = vcf_fields[9].strip().split(":")
    if baf_depth[genotype_index] != "0/1":
        return None
    depth = baf_depth[depth_index]
    if frequency_index > 0:
        freq = baf_depth[frequency_index]
    elif depth_index > 0:
        freq = float(baf_depth[other_depth_index]) / int(depth)
    else:
        sys.exit("No AO or AF fields in VCF format column. Quit.")
    if int(depth) > 1:
        return [chromosome, start, end, reference_variant, allelic_variant[0], chromosome + "_" + start,
                freq, depth]


def parse_vcf(input_file_name: str, output_file_name: str, minimum_quality: float) -> None:
    if input_file_name.endswith(".gz"):
        input_file = gzip.open(input_file_name)
    elif input_file_name.endswith((".vcf", ".VCF")):
        input_file = open(input_file_name)
    header = "#chr\tstart\tend\tref\tobs\tid\tfreq\tdepth\n"
    with open(output_file_name, "w") as output_file:
        output_file.write(header)
        for line in input_file:
            parsed_line = vcf_parse_one_line(line, minimum_quality)
            if parsed_line:
                output_file.write("\t".join(map(str, parsed_line)) + "\n")


def main() -> None:
    input_file_name = sys.argv[1]
    output_file_name = sys.argv[2]
    minimum_quality = float(sys.argv[3])
    parse_vcf(input_file_name, output_file_name, minimum_quality)


if __name__ == '__main__':
    main()
