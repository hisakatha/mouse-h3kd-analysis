#!/usr/bin/env python3
import sys

count_prefix = "Aligned.sortedByCoord.out.uniq_collated.template.sorted.nooverhang_outside_gene2_slop"
count_suffix_chr_all = ".bed.count"
count_suffix_chr_main = ".bed.count_main"
total_count_file_chr_all = "Aligned.sortedByCoord.out.uniq_collated.template.bed.count"
total_count_file_chr_main = "Aligned.sortedByCoord.out.uniq_collated.template.bed.count_main"
output_path = "template_count_outside_gene2.summarize.csv"
output = open(output_path, 'w', encoding='utf-8')

widths = ["0", "100", "500", "1k", "5k", "10k", "50k", "100k", "500k", "1m"]
widths_numeric = {"0":0, "100":100, "500":500, "1k":1000, "5k":5000, "10k":10000, "50k":50000, "100k":100000, "500k":500000, "1m":1000000}

sample_names = {
    "NoRibs_MII201209": "MII #1",
    "NoRibs_MII201218": "MII #2",
    "NoRibs_2c_18hpi201209": "2C 18hpi #1",
    "NoRibs_2c_18hpi201218": "2C 18hpi #2",
    "NoRibs_2c_Noinjection201209": "2C No Injection #1",
    "NoRibs_2c_Noinjection201218": "2C No Injection #2",
    "NoRibs_2c_Control201209": "2C Control #1",
    "NoRibs_2c_Control201218": "2C Control #2",
    "NoRibs_2c_H3_1_3_2KD201209": "2C H3.1/3.2 KD #1",
    "NoRibs_2c_H3_1_3_2KD201218": "2C H3.1/3.2 KD #2",
    "NoRibs_2cAphidicolin201209": "2C DNA(-) #1",
    "NoRibs_2c_Aphidicoilin201218": "2C DNA(-) #2",
}
sample_dirs = list(sample_names)
sample_categories = {
    "NoRibs_MII201209": "MII",
    "NoRibs_MII201218": "MII",
    "NoRibs_2c_18hpi201209": "2C 18hpi",
    "NoRibs_2c_18hpi201218": "2C 18hpi",
    "NoRibs_2c_Noinjection201209": "2C No Injection",
    "NoRibs_2c_Noinjection201218": "2C No Injection",
    "NoRibs_2c_Control201209": "2C Control",
    "NoRibs_2c_Control201218": "2C Control",
    "NoRibs_2c_H3_1_3_2KD201209": "2C H3.1/3.2 KD",
    "NoRibs_2c_H3_1_3_2KD201218": "2C H3.1/3.2 KD",
    "NoRibs_2cAphidicolin201209": "2C DNA(-)",
    "NoRibs_2c_Aphidicoilin201218": "2C DNA(-)",
}

print("sample_dir,sample_name,sample_category,chr_set,gene_extension_width,gene_extension_width_text,intergenic_read_count,total_read_count,intergenic_read_ratio",
    file=output)
output_format = "{:s},{:s},{:s},{:s},{:d},{:s},{:d},{:d},{:.3g}"

for width in widths:
    for sample_dir in sample_dirs:
        count_path_chr_all = sample_dir + "/" + count_prefix + width + count_suffix_chr_all
        with open(count_path_chr_all, 'r') as f:
            lines = [s.strip() for s in f.readlines()]
            if len(lines) > 1:
                print("[ERROR] Unexpected lines in {}".format(f.name), file=sys.stderr)
                sys.exit(1)
            else:
                count_chr_all = int(lines[0])

        total_count_path_chr_all = sample_dir + "/" + total_count_file_chr_all
        with open(total_count_path_chr_all, 'r') as f:
            lines = [s.strip() for s in f.readlines()]
            if len(lines) > 1:
                print("[ERROR] Unexpected lines in {}".format(f.name), file=sys.stderr)
                sys.exit(1)
            else:
                total_count_chr_all = int(lines[0])

        print(output_format.format(sample_dir, sample_names[sample_dir], sample_categories[sample_dir], "all",
            widths_numeric[width], width, count_chr_all, total_count_chr_all, count_chr_all / total_count_chr_all), file=output)

        count_path_chr_main = sample_dir + "/" + count_prefix + width + count_suffix_chr_main
        with open(count_path_chr_main, 'r') as f:
            lines = [s.strip() for s in f.readlines()]
            if len(lines) > 1:
                print("[ERROR] Unexpected lines in {}".format(f.name), file=sys.stderr)
                sys.exit(1)
            else:
                count_chr_main = int(lines[0])
        
        total_count_path_chr_main = sample_dir + "/" + total_count_file_chr_main
        with open(total_count_path_chr_main, 'r') as f:
            lines = [s.strip() for s in f.readlines()]
            if len(lines) > 1:
                print("[ERROR] Unexpected lines in {}".format(f.name), file=sys.stderr)
                sys.exit(1)
            else:
                total_count_chr_main = int(lines[0])

        print(output_format.format(sample_dir, sample_names[sample_dir], sample_categories[sample_dir], "main",
            widths_numeric[width], width, count_chr_main, total_count_chr_main, count_chr_main / total_count_chr_main), file=output)

output.close()
