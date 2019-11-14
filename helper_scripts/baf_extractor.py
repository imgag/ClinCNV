__author__ = 'gdemidov'

import gzip
import sys
from typing import List, Optional


def vcf_parse_one_line(line: str, minimum_quality: float) -> Optional[List]:
    if line.startswith("#"):
        return None

    format_fields = line.split("\t")[8].split(":")
    index_for_genotype = format_fields.index("GT")
    index_for_depth = format_fields.index("DP")
    try:
        index_for_other_depth = format_fields.index("AO")
    except ValueError:
        index_for_other_depth = -1
    try:
        index_for_frequency = format_fields.index("AF")
    except ValueError:
        index_for_frequency = -1

    splitted_line = line.split("\t")
    chrom = splitted_line[0]
    start = splitted_line[1]
    end = splitted_line[1]
    reference_variant = splitted_line[3]
    allelic_variant = splitted_line[4].split(",")
    if len(allelic_variant[0]) > 1 or len(reference_variant) > 1 or len(allelic_variant) > 2:
        return None
    quality = float(splitted_line[5])
    if quality < minimum_quality:
        return None
    baf_depth = splitted_line[9].strip().split(":")
    if baf_depth[index_for_genotype] == "0/1":
        depth = baf_depth[index_for_depth]
        if index_for_frequency > 0:
            freq = baf_depth[index_for_frequency]
        elif index_for_depth > 0:
            freq = float(baf_depth[index_for_other_depth]) / int(depth)
        else:
            sys.exit("No AO or AF fields in VCF format column. Quit.")
        if int(depth) > 1:
            return [chrom, start, end, reference_variant, allelic_variant[0], chrom + "_" + start,
                    freq, depth]
    return None


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
