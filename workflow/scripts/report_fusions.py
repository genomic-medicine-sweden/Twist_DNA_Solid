
import operator
import gzip

input_bed = open(snakemake.input.bed)
input_bed_extra_annotation = open(snakemake.input.bed_extra_annotation)
dedup_coverage_filename = snakemake.input.dedup_coverage
input_arriba = open(snakemake.input.arriba)
input_starfusion = open(snakemake.input.star_fusion)
input_fusioncatcher = open(snakemake.input.fusioncatcher)
input_bam = snakemake.input.bam
output_fusions = open(snakemake.output.fusions, "w")
star_fusion_flag_low_support = snakemake.params.star_fusion_flag_low_support
star_fusion_low_support = snakemake.params.star_fusion_low_support
star_fusion_low_support_inframe = snakemake.params.star_fusion_low_support_inframe
fusioncatcher_flag_low_support = snakemake.params.fusioncatcher_flag_low_support
fusioncatcher_low_support = snakemake.params.fusioncatcher_low_support
fusioncatcher_low_support_inframe = snakemake.params.fusioncatcher_low_support_inframe
fp_fusions_filename = snakemake.params.fp_fusions

output_fusions.write("callers\tgene1\tgene2\texon1\texon2\tconfidence\tFC-callers\tpredicted_effect\tbreakpoint1\tbreakpoint2\t")
output_fusions.write("dedup_coverage\tA_split_reads\tA_spanning_pairs\tA_total_supporting_reads\t")
output_fusions.write("SF_split_reads\tSF_spanning_pairs\tSF_total_supporting_reads\tFC_split_reads\tFC_spanning_pairs\t")
output_fusions.write("FC_total_supporting_reads\tAll_total_supporting_reads\n")

# Filter noisy genes and housekeeping genes involved in fusions
housekeeping_genes = {}
artefact_gene_dict = {}
if fp_fusions_filename != "":
    with open(fp_fusions_filename) as fp_fusions:
        for line in fp_fusions:
            if line.startswith("#"):
                continue
            columns = line.strip().split("\t")
            gene1 = columns[0]
            gene2 = columns[1]
            read_limit_SF = int(columns[2])
            read_limit_FC = int(columns[3])
            af_limit = float(columns[4])
            if gene2 == "housekeeping":
                housekeeping_genes[gene1] = [read_limit_SF, read_limit_FC, af_limit, 0]
            if gene1 not in artefact_gene_dict:
                artefact_gene_dict[gene1] = {}
            artefact_gene_dict[gene1][gene2] = [read_limit_SF, read_limit_FC, af_limit, 0]

# Only keep fusions with one gene that are in the design
design_genes = {}
for line in input_bed:
    lline = line.strip().split("\t")
    chrom = lline[0]
    start = int(lline[1])
    end = int(lline[2])
    exon = lline[3]
    gene = lline[3].split("_")[0]
    if gene in design_genes:
        design_genes[gene].append([chrom, start, end, exon])
    else:
        design_genes[gene] = [[chrom, start, end, exon]]

# Extra annotation of partner genes
annotation_genes = {}
for line in input_bed_extra_annotation:
    lline = line.strip().split("\t")
    chrom = lline[0]
    start = int(lline[1])
    end = int(lline[2])
    exon = lline[3]
    gene = lline[3].split("_")[0]
    if gene in annotation_genes:
        annotation_genes[gene].append([chrom, start, end, exon])
    else:
        annotation_genes[gene] = [[chrom, start, end, exon]]

# Deduplicated coverage of fusion regions
dedup_coverage_list = []
with gzip.open(dedup_coverage_filename, 'rt') as dedup_coverage:
    for line in dedup_coverage:
        columns = line.strip().split("\t")
        chrom = columns[0]
        start_pos = int(columns[1])
        end_pos = int(columns[2])
        coverage = round(float(columns[4]))
        dedup_coverage_list.append([chrom, start_pos, end_pos, coverage])

fusion_dict = {}

# Arriba fusions
header = True
for line in input_arriba:
    if header:
        header = False
        continue
    lline = line.strip().split("\t")
    gene1 = lline[0]
    gene2 = lline[1]
    # Only keep fusions with one gene that are in the design
    if (gene1 not in design_genes and gene2 not in design_genes):
        continue
    confidence = lline[14]
    breakpoint1 = lline[4]
    breakpoint2 = lline[5]
    split_reads1 = lline[9]
    split_reads2 = lline[10]
    total_split_reads = int(split_reads1) + int(split_reads2)
    discordant_mates = lline[11]
    coverage1 = lline[12]
    coverage2 = lline[13]
    predicted_effect = lline[15]
    # Compare fusion coverage with coverage in breakpoints
    chrom1 = "chr" + breakpoint1.split(":")[0]
    pos1 = breakpoint1.split(":")[1]
    chrom2 = "chr" + breakpoint2.split(":")[0]
    pos2 = breakpoint2.split(":")[1]
    # Get exon name if it is in design
    exon1 = ""
    exon2 = ""
    if gene1 in design_genes:
        for region in design_genes[gene1]:
            if int(pos1) >= region[1] and int(pos1) <= region[2]:
                exon1 = region[3]
    elif gene1 in annotation_genes:
        for region in annotation_genes[gene1]:
            if int(pos1) >= region[1] and int(pos1) <= region[2]:
                exon1 = region[3]
    if gene2 in design_genes:
        for region in design_genes[gene2]:
            if int(pos2) >= region[1] and int(pos2) <= region[2]:
                exon2 = region[3]
    elif gene2 in annotation_genes:
        for region in annotation_genes[gene2]:
            if int(pos2) >= region[1] and int(pos2) <= region[2]:
                exon2 = region[3]
    total_supporting_reads = int(total_split_reads) + int(discordant_mates)
    break_points = breakpoint1 + "_" + breakpoint2
    if break_points not in fusion_dict:
        fusion_dict[break_points] = {}
    fusion_dict[break_points]["Arriba"] = [gene1, gene2, exon1, exon2, confidence, "", predicted_effect, breakpoint1, breakpoint2,
                                           total_split_reads, discordant_mates, total_supporting_reads]
    # dedup coverage for potential filtering in StarFusion and FusionCatcher
    exon_coverage1 = 0
    exon_coverage2 = 0
    max_exon_coverage = 0
    for exon in dedup_coverage_list:
        if chrom1 == exon[0] and pos1 >= exon[1] and pos1 <= exon[2]:
            exon_coverage1 = exon[3]
        if chrom2 == exon[0] and pos2 >= exon[1] and pos2 <= exon[2]:
            exon_coverage2 = exon[3]
    max_exon_coverage = max(exon_coverage1, exon_coverage2)
    if gene1 in artefact_gene_dict and gene2 in artefact_gene_dict[gene1]:
        artefact_gene_dict[gene1][gene2][3] = max_exon_coverage
    if gene2 in artefact_gene_dict and gene1 in artefact_gene_dict[gene2]:
        artefact_gene_dict[gene2][gene1][3] = max_exon_coverage


# Star-fusions
header = True
for line in input_starfusion:
    if header:
        header = False
        continue
    lline = line.strip().split("\t")
    gene1 = lline[0].split("--")[0]
    gene2 = lline[0].split("--")[1]
    # Only keep fusions with one gene that are in the design
    if (gene1 not in design_genes and gene2 not in design_genes):
        continue
    Junction_read_count = lline[1]
    Spanning_Frag_count = lline[2]
    predicted_effect = lline[21]
    # Flag fusions with junction_read_count < 15 and Spanning_Frag_count < 2
    confidence = ""
    if int(Junction_read_count) < star_fusion_flag_low_support:
        confidence = "Low support"
    # Remove Fusions with very weak read support
    if int(Junction_read_count) <= star_fusion_low_support_inframe and predicted_effect != "INFRAME":
        continue
    if int(Junction_read_count) <= star_fusion_low_support:
        continue
    # Higher demand of read support for genes with frequent FP gene fusions and house keeping genes
    if (gene1 in artefact_gene_dict and gene2 in artefact_gene_dict[gene1]):
        if int(Junction_read_count) < artefact_gene_dict[gene1][gene2][0]:
            continue
    if (gene2 in artefact_gene_dict and gene1 in artefact_gene_dict[gene2]):
        if int(Junction_read_count) < artefact_gene_dict[gene2][gene1][0]:
            continue
    if gene1 in housekeeping_genes:
        if int(Junction_read_count) < housekeeping_genes[gene1][0]:
            continue
    if gene2 in housekeeping_genes:
        if int(Junction_read_count) < housekeeping_genes[gene2][0]:
            continue
    # Min AF for frequent FP gene fusions and housekeeping gene
    if (gene1 in artefact_gene_dict and gene2 in artefact_gene_dict[gene1] and artefact_gene_dict[gene1][gene2][3] > 0):
        if int(Junction_read_count) / artefact_gene_dict[gene1][gene2][3] < artefact_gene_dict[gene1][gene2][2]:
            continue
    if (gene2 in artefact_gene_dict and gene1 in artefact_gene_dict[gene2] and artefact_gene_dict[gene2][gene1][3] > 0):
        if int(Junction_read_count) / artefact_gene_dict[gene2][gene1][3] < artefact_gene_dict[gene2][gene1][2]:
            continue
    if gene1 in housekeeping_genes and housekeeping_genes[gene1][3] > 0:
        if int(Junction_read_count) / housekeeping_genes[gene1][3] < housekeeping_genes[gene1][2]:
            continue
    if gene2 in housekeeping_genes and housekeeping_genes[gene2][3] > 0:
        if int(Junction_read_count) / housekeeping_genes[gene2][3] < housekeeping_genes[gene2][2]:
            continue
    breakpoint1 = lline[7][:-2]
    breakpoint2 = lline[9][:-2]
    FFPM = lline[11]
    DBs = lline[16]
    # Compare fusion coverage with coverage in breakpoints
    chrom1 = breakpoint1.split(":")[0]
    pos1 = breakpoint1.split(":")[1]
    chrom2 = breakpoint2.split(":")[0]
    pos2 = breakpoint2.split(":")[1]
    # Get exon name if it is in design
    exon1 = ""
    exon2 = ""
    if gene1 in design_genes:
        for region in design_genes[gene1]:
            if int(pos1) >= region[1] and int(pos1) <= region[2]:
                exon1 = region[3]
    elif gene1 in annotation_genes:
        for region in annotation_genes[gene1]:
            if int(pos1) >= region[1] and int(pos1) <= region[2]:
                exon1 = region[3]
    if gene2 in design_genes:
        for region in design_genes[gene2]:
            if int(pos2) >= region[1] and int(pos2) <= region[2]:
                exon2 = region[3]
    elif gene2 in annotation_genes:
        for region in annotation_genes[gene2]:
            if int(pos2) >= region[1] and int(pos2) <= region[2]:
                exon2 = region[3]
    total_supporting_reads = int(Junction_read_count) + int(Spanning_Frag_count)
    break_points = breakpoint1 + "_" + breakpoint2
    if break_points not in fusion_dict:
        fusion_dict[break_points] = {}
    fusion_dict[break_points]["StarFusion"] = [gene1, gene2, exon1, exon2, confidence, "", predicted_effect, breakpoint1,
                                               breakpoint2, Junction_read_count, Spanning_Frag_count, total_supporting_reads]


# FusionCatcher
header = True
for line in input_fusioncatcher:
    if header:
        header = False
        continue
    lline = line.strip().split("\t")
    gene1 = lline[0]
    gene2 = lline[1]
    # Only keep fusions with one gene that are in the design
    if (gene1 not in design_genes and gene2 not in design_genes):
        continue
    fp_filters = lline[2].split(",")
    DBs = lline[2]
    common_mapping = lline[3]
    Spanning_pairs = lline[4]
    Spanning_reads_unique = lline[5]
    Fusion_finding_method = lline[7]
    breakpoint1 = f"chr{lline[8][:-2]}"
    breakpoint2 = f"chr{lline[9][:-2]}"
    predicted_effect = lline[15]
    # Flag fusions with Spanning_reads_unique < 5
    confidence = ""
    if int(Spanning_reads_unique) < fusioncatcher_flag_low_support:
        confidence = "Low support"
    # Filter fusions with very low support
    if int(Spanning_reads_unique) <= fusioncatcher_low_support_inframe and predicted_effect != "in-frame":
        continue
    if int(Spanning_reads_unique) <= fusioncatcher_low_support:
        continue
    # Higher demand of read support for genes with frequent FP gene fusions and house keeping genes
    if (gene1 in artefact_gene_dict and gene2 in artefact_gene_dict[gene1]):
        if int(Spanning_reads_unique) < artefact_gene_dict[gene1][gene2][1]:
            continue
    if (gene2 in artefact_gene_dict and gene1 in artefact_gene_dict[gene2]):
        if int(Spanning_reads_unique) < artefact_gene_dict[gene2][gene1][1]:
            continue
    if gene1 in housekeeping_genes:
        if int(Spanning_reads_unique) < housekeeping_genes[gene1][1]:
            continue
    if gene2 in housekeeping_genes:
        if int(Spanning_reads_unique) < housekeeping_genes[gene2][1]:
            continue
    # Min AF for frequent FP gene fusions and housekeeping gene
    if (gene1 in artefact_gene_dict and gene2 in artefact_gene_dict[gene1] and artefact_gene_dict[gene1][gene2][3] > 0):
        if int(Spanning_reads_unique) / artefact_gene_dict[gene1][gene2][3] < artefact_gene_dict[gene1][gene2][2]:
            continue
    if (gene2 in artefact_gene_dict and gene1 in artefact_gene_dict[gene2] and artefact_gene_dict[gene2][gene1][3] > 0):
        if int(Spanning_reads_unique) / artefact_gene_dict[gene2][gene1][3] < artefact_gene_dict[gene2][gene1][2]:
            continue
    if gene1 in housekeeping_genes and housekeeping_genes[gene1][3] > 0:
        if int(Spanning_reads_unique) / housekeeping_genes[gene1][3] < housekeeping_genes[gene1][2]:
            continue
    if gene2 in housekeeping_genes and housekeeping_genes[gene2][3] > 0:
        if int(Spanning_reads_unique) / housekeeping_genes[gene2][3] < housekeeping_genes[gene2][2]:
            continue
    # Flag fusions annotated that are fusions with very high probability
    fp_db = [
        "banned", "bodymap2", "cacg", "1000genomes", "conjoing", "cortex", "distance1000bp", "ensembl_fully_overlapping",
        "ensembl_same_strand_overlapping", "gtex", "hpa", "mt", "paralogs", "refseq_fully_overlapping",
        "refseq_same_strand_overlapping", "rrna", "similar_reads", "similar_symbols", "ucsc_fully_overlapping",
        "ucsc_same_strand_overlapping"
    ]
    fp_found = ""
    for fp in fp_db:
        if fp in fp_filters:
            fp_found = "FP"
    # Compare fusion coverage with coverage in breakpoints
    pos1 = "0"
    pos2 = "0"
    if len(breakpoint1.split(":")) == 2 and len(breakpoint2.split(":")) == 2:
        chrom1 = "chr" + breakpoint1.split(":")[0]
        pos1 = breakpoint1.split(":")[1]
        chrom2 = "chr" + breakpoint2.split(":")[0]
        pos2 = breakpoint2.split(":")[1]
    # Get exon name if it is in design
    exon1 = ""
    exon2 = ""
    if gene1 in design_genes:
        for region in design_genes[gene1]:
            if int(pos1) >= region[1] and int(pos1) <= region[2]:
                exon1 = region[3]
    elif gene1 in annotation_genes:
        for region in annotation_genes[gene1]:
            if int(pos1) >= region[1] and int(pos1) <= region[2]:
                exon1 = region[3]
    if gene2 in design_genes:
        for region in design_genes[gene2]:
            if int(pos2) >= region[1] and int(pos2) <= region[2]:
                exon2 = region[3]
    elif gene2 in annotation_genes:
        for region in annotation_genes[gene2]:
            if int(pos2) >= region[1] and int(pos2) <= region[2]:
                exon2 = region[3]
    total_supporting_reads = int(Spanning_pairs) + int(Spanning_reads_unique)
    break_points = breakpoint1 + "_" + breakpoint2
    if break_points not in fusion_dict:
        fusion_dict[break_points] = {}
    fusion_dict[break_points]["FusionCatcher"] = [gene1, gene2, exon1, exon2, confidence, Fusion_finding_method,
                                                  predicted_effect, breakpoint1, breakpoint2, Spanning_pairs,
                                                  Spanning_reads_unique, total_supporting_reads]


merged_fusions = []
i = 0
for break_points in fusion_dict:

    exon_coverage1 = 0
    exon_coverage2 = 0
    max_exon_coverage = 0
    caller = ""
    if "Arriba" in fusion_dict[break_points]:
        caller = "Arriba"
    elif "StarFusion" in fusion_dict[break_points]:
        caller = "StarFusion"
    elif "FusionCatcher" in fusion_dict[break_points]:
        caller = "FusionCatcher"
    breakpoint1 = fusion_dict[break_points][caller][7]
    breakpoint2 = fusion_dict[break_points][caller][8]
    flag_coverage1 = False
    flag_coverage2 = False
    if len(breakpoint1.split(":")) == 2:
        chrom1 = breakpoint1.split(":")[0]
        pos1 = int(breakpoint1.split(":")[1])
    else:
        flag_coverage1 = True
    if len(breakpoint2.split(":")) == 2:
        chrom2 = breakpoint2.split(":")[0]
        pos2 = int(breakpoint2.split(":")[1])
    else:
        flag_coverage2 = True
    for exon in dedup_coverage_list:
        if chrom1 == exon[0] and pos1 >= exon[1] and pos1 <= exon[2] and not flag_coverage1:
            exon_coverage1 = exon[3]
        if chrom2 == exon[0] and pos2 >= exon[1] and pos2 <= exon[2] and not flag_coverage2:
            exon_coverage2 = exon[3]
    max_exon_coverage = max(exon_coverage1, exon_coverage2)

    first = True
    if "Arriba" in fusion_dict[break_points]:
        data = fusion_dict[break_points]["Arriba"]
        merged_fusions.append([data[0], data[1], data[2], data[3], data[4], "", data[6], data[7], data[8], data[9], data[10],
                               data[11], "", "", "", "", "", "", int(data[11]), 1, "Arriba", max_exon_coverage])
        first = False
    if "StarFusion" in fusion_dict[break_points]:
        data = fusion_dict[break_points]["StarFusion"]
        if first:
            merged_fusions.append([data[0], data[1], data[2], data[3], data[4], "", data[6], data[7], data[8], "", "", "",
                                   data[9], data[10], data[11], "", "", "", int(data[11]), 1, "StarFusion", max_exon_coverage])
            first = False
        else:
            merged_fusions[i][4] += f", {data[4]}"
            merged_fusions[i][6] += f", {data[6]}"
            merged_fusions[i][12] = data[9]
            merged_fusions[i][13] = data[10]
            merged_fusions[i][14] = data[11]
            merged_fusions[i][18] += int(data[11])
            merged_fusions[i][19] += 1
            merged_fusions[i][20] += ", StarFusion"
    if "FusionCatcher" in fusion_dict[break_points]:
        data = fusion_dict[break_points]["FusionCatcher"]
        if first:
            merged_fusions.append([data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[8], "", "", "",
                                   "", "", "", data[9], data[10], data[11], int(data[11]), 1, "FusionCatcher", max_exon_coverage])
            first = False
        else:
            merged_fusions[i][4] += f", {data[4]}"
            merged_fusions[i][5] += data[5]
            merged_fusions[i][6] += f", {data[6]}"
            merged_fusions[i][15] = data[9]
            merged_fusions[i][16] = data[10]
            merged_fusions[i][17] = data[11]
            merged_fusions[i][18] += int(data[11])
            merged_fusions[i][19] += 1
            merged_fusions[i][20] += ", FusionCatcher"
    i += 1

merged_fusions.sort(key=operator.itemgetter(18), reverse=True)

for fusion in merged_fusions:
    output_fusions.write(f"{fusion[20]}\t{fusion[0]}\t{fusion[1]}\t{fusion[2]}\t{fusion[3]}\t{fusion[4]}\t{fusion[5]}\t")
    output_fusions.write(f"{fusion[6]}\t{fusion[7]}\t{fusion[8]}\t{fusion[21]}\t{fusion[9]}\t{fusion[10]}\t{fusion[11]}\t")
    output_fusions.write(f"{fusion[12]}\t{fusion[13]}\t{fusion[14]}\t{fusion[15]}\t{fusion[16]}\t{fusion[17]}\t")
    output_fusions.write(f"{fusion[18]}\n")
