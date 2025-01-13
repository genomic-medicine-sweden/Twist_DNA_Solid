
import os

vcf_dna_filenames = snakemake.input.id_snp_vcf_dna
vcf_rna_filenames = snakemake.input.id_snp_vcf_rna
report = open(snakemake.output.mixup_report, "w")
match_cutoff = float(snakemake.params.match_cutoff)

vcf_dict_dna = {}
dna_samples = {}
vcf_dict_rna = {}
rna_samples = {}


def read_vcf(vcf_filename, vcf_dict, samples):
    vcf = open(vcf_filename)
    sample_type = os.path.split(vcf_filename)[1].split(".")[0]
    vcf_dict[sample_type] = {}
    samples[sample_type] = {}
    header = True
    for line in vcf:
        if header:
            if line[:6] == "#CHROM":
                header = False
            continue
        columns = line.strip().split("\t")
        chrom = columns[0][3:]
        if chrom == "X":
            continue
        pos = columns[1]
        key = chrom + "_" + pos
        FORMAT = columns[8].split(":")
        DATA = columns[9].split(":")
        GT_index = 0
        i = 0
        for f in FORMAT:
            if f == "GT":
                GT_index = i
            i += 1
        GT = DATA[GT_index]
        vcf_dict[sample_type][key] = GT
    vcf.close()


for vcf_filename in vcf_dna_filenames:
    read_vcf(vcf_filename, vcf_dict_dna, dna_samples)
for vcf_filename in vcf_rna_filenames:
    read_vcf(vcf_filename, vcf_dict_rna, rna_samples)

for rna_sample in rna_samples:
    for dna_sample in dna_samples:
        rna_samples[rna_sample][dna_sample] = 0
        for key in vcf_dict_rna[rna_sample]:
            GT_rna = vcf_dict_rna[rna_sample][key]
            if key in vcf_dict_dna[dna_sample]:
                GT_dna = vcf_dict_dna[dna_sample][key]
                if GT_rna == GT_dna:
                    rna_samples[rna_sample][dna_sample] += 1

report.write("RNA_sample\tDNA_sample\tnr_matches\t%_match\tmatch\n")
for rna_sample in rna_samples:
    best_dna_sample = ""
    best_gt_match = 0
    for dna_sample in rna_samples[rna_sample]:
        if rna_samples[rna_sample][dna_sample] > best_gt_match:
            best_dna_sample = dna_sample
            best_gt_match = rna_samples[rna_sample][dna_sample]
    p_match = round(best_gt_match * 100 / 42.0, 1)
    report.write(f"{rna_sample}\t{best_dna_sample}\t{best_gt_match}\t{p_match}%\t")
    if p_match > match_cutoff * 100:
        report.write(f"yes\n")
    else:
        report.write(f"no\n")
