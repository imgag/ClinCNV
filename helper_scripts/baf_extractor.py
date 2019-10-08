__author__ = 'gdemidov'
import sys
import gzip

def vcf_parse_one_line(line, minimum_qual):
    if not line.startswith("#"):
        splitted_line = line.split("\t")
        chrom = splitted_line[0]
        start = splitted_line[1]
        end = splitted_line[1]
        reference_variant = splitted_line[3]
        allelic_variant = splitted_line[4].split(",")
        if len(allelic_variant[0]) > 1 or len(reference_variant) > 1 or len(allelic_variant) > 2:
            return 0
        quality = float(splitted_line[5])
        if quality < minimum_qual:
            return 0
        bafDepth = splitted_line[9].split(":")
        if bafDepth[0] == "0/1":
            depth = bafDepth[3]
            freq = bafDepth[1]
            if int(depth) > 30:
                return [chrom, start, end, reference_variant, allelic_variant[0], chrom + "_" + start,
                        freq, depth]
        return 0
    return 0


def parse_vcf(input_file_name, output_file_name, minimum_qual):
    with gzip.open(input_file_name, 'r') as f:
        header = "#chr\tstart\tend\tref\tobs\tid\tfreq\tdepth\n"
        with open(output_file_name, "w") as f1:
            f1.write(header)
            for line in f:
                parsed_line = vcf_parse_one_line(line, minimum_qual)
                if parsed_line != 0:
                    f1.write("\t".join(map(str, parsed_line)) + "\n")


def main():
    input_file_name = sys.argv[1]
    output_file_name = sys.argv[2]
    minimum_qual = float(sys.argv[3])
    parse_vcf(input_file_name, output_file_name, minimum_qual)

main()