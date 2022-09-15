
fusions = open(snakemake.input.fusions)
report = open(snakemake.output.report, "w")
min_unique_reads = snakemake.params.min_unique_reads

report.write("Gene1\tGene2\tNr_unique_reads\tGene_region1\tBreak_point1\tTranscript1\tGene_region2\tBreak_point2\tTranscript2\n")

FP_gene_pairs = ["NPM1_ALK", "CLTC_NTRK3", "MSH2_ALK", "MSH2_HIP1"]
Noisy_genes_pairs = {"LMNA_EZR": 9, "ABL1_STRN": 7, "EZR_ALK": 8, "RSPO2_BRAF": 8, "LMNA:HIP1" : 12}

for line in fusions:
    if line[:8] != "#Fusion:":
        continue
    gene1 = line.split(" ")[1].split("_")[0]
    gene2 = line.split("___")[1].split("_")[0]
    unique_reads = int(line.split("unique:")[1].split(")")[0])
    gene_region1 = ""
    if "intron" in line.split("___")[0]:
        gene_region1 = "intron " + line.split("___")[0].split("intron:")[1].split("|")[0]
    if "exon" in line.split("___")[0]:
        gene_region1 = "exon " + line.split("___")[0].split("exon:")[1].split("|")[0]
    gene_region2 = ""
    if "intron" in line.split("___")[1]:
        gene_region2 = "intron " + line.split("___")[1].split("intron:")[1].split("|")[0]
    if "exon" in line.split("___")[1]:
        gene_region2 = "exon " + line.split("___")[1].split("exon:")[1].split("|")[0]
    bp1 = "chr" + line.split("___")[0].split("chr")[1]
    bp2 = "chr" + line.split("___")[1].split("chr")[1].split("  (")[0]
    transcript1 = line.split("_")[1].split(":")[0]
    transcript2 = line.split("___")[1].split("_")[1].split(":")[0]
    min_u_r = min_unique_reads
    key = gene1 + "_" + gene2
    if key in Noisy_genes_pairs:
        min_u_r = Noisy_genes_pairs[key]
    if unique_reads < min_u_r:
        continue
    if key in FP_gene_pairs:
        continue
    report.write(
        f"{gene1}\t{gene2}\t{unique_reads}\t{gene_region1}\t{bp1}\t{transcript1}\t{gene_region2}\t{bp2}\t{transcript2}\n"
    )
report.close()
