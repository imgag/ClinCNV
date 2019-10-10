__author__ = 'gdemidov'
import sys
import gzip

def vcf_parse_one_line(line, minimum_qual):
    index_for_genotype = 0
    index_for_depth = 0
    index_for_other_depth = -1
    index_for_frequency = -1
    if not line.startswith("#"):
        format_fields = line.split("\t")[8].split(":")
        index_for_genotype = format_fields.index("GT")
        index_for_depth = format_fields.index("DP")
        try:
            index_for_other_depth = format_fields.index("AO")
        except:
            pass
        try:
            index_for_frequency = format_fields.index("AF")
        except:
            pass
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
        bafDepth = splitted_line[9].strip().split(":")
        if bafDepth[index_for_genotype] == "0/1":
            depth = bafDepth[index_for_depth]
            if index_for_frequency > 0:
                freq = bafDepth[index_for_frequency]
            elif index_for_depth > 0:
                freq = float(bafDepth[index_for_other_depth]) / int(depth)
            else:
                print("No AO or AF fields in VCF format column. Quit.")
                quit()
            if int(depth) > 1:
                return [chrom, start, end, reference_variant, allelic_variant[0], chrom + "_" + start,
                        freq, depth]
        return 0
    return 0


def parse_vcf(input_file_name, output_file_name, minimum_qual):
    if input_file_name.endswith(".gz"):
        f = gzip.open(input_file_name)
    elif input_file_name.endswith(".vcf") or input_file_name.endswith(".VCF"):
        f = open(input_file_name)
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
