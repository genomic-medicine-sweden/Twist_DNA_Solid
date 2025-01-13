import logging
import statistics
from pysam import VariantFile
from hydra_genetics.utils.io import utils

log = logging.getLogger()


def check_fp(chrom, start, end, gatk_cnr_dict, cn, max_cnv_fp_size):
    FP_flag = "-"
    cnv_length = end - start + 1
    if cnv_length > max_cnv_fp_size:
        return "-"
    gatk_data = gatk_cnr_dict[chrom]
    i = 0
    j = 0
    cnv_region = {"surrounding_region": [], "region": []}
    for region in gatk_data:
        if region[1] < start:
            i += 1
            j += 1
        elif region[0] > end:
            break
        else:
            cnv_region["region"].append(region[2])
            j += 1
    i1 = max(0, i-11)
    i2 = max(0, i-1)
    while i1 < i2:
        cnv_region["surrounding_region"].append(gatk_data[i1][2])
        i1 += 1
    j1 = min(len(gatk_data), j+1)
    j2 = min(len(gatk_data), j+11)
    while j1 < j2:
        cnv_region["surrounding_region"].append(gatk_data[j1][2])
        j1 += 1

    median_surrounding_region = statistics.median(cnv_region["surrounding_region"])
    stdev_surrounding_region = 0.0
    if len(cnv_region["surrounding_region"]):
        stdev_surrounding_region = statistics.stdev(cnv_region["surrounding_region"])
    median_region = statistics.median(cnv_region["region"])
    stdev_region = 0.0
    if len(cnv_region["region"]) >= 3:
        stdev_region = statistics.stdev(cnv_region["region"])

    if cn < 1.5 and median_region > -0.2 and median_region < 0.2:
        if median_region + stdev_region > median_surrounding_region - stdev_surrounding_region:
            FP_flag = "FP"
    if cn > 2.5 and median_region > -0.2 and median_region < 0.2:
        if median_region - stdev_region < median_surrounding_region + stdev_surrounding_region:
            FP_flag = "FP"

    return FP_flag


def create_tsv_report(
    input_vcfs, input_org_vcfs, input_del, input_amp, in_chrom_arm_size, in_gatk_cnr, amp_cn_limit,
    output_txt, out_additional_only, out_tsv_chrom_arms, out_vcf_filename, del_1p19q_cn, del_1p19q_chr_arm_fraction,
    chr_arm_fraction, del_chr_arm_cn_limit, amp_chr_arm_cn_limit, normal_cn_lower_limit, normal_cn_upper_limit,
    normal_baf_lower_limit, normal_baf_upper_limit, baseline_fraction_limit, polyploidy_fraction_limit, TC, max_cnv_fp_size,
):
    chrom_arm_size = {}
    chrom_arm_del = {}
    chrom_arm_loh = {}
    chrom_arm_amp = {}
    baseline = [0, 0]
    polyploidy = 0
    genome_size = 0
    with open(in_chrom_arm_size) as chrom_arm_size_file:
        next(chrom_arm_size_file)
        for line in chrom_arm_size_file:
            columns = line.strip().split("\t")
            chrom = columns[0]
            p_size = int(columns[1])
            q_size = int(columns[2])
            chrom_arm_size[chrom] = [[0, p_size, p_size], [p_size, p_size + q_size, q_size]]
            chrom_arm_del[chrom] = [0, 0]
            chrom_arm_loh[chrom] = [0, 0]
            chrom_arm_amp[chrom] = [0, 0]
            genome_size += p_size + q_size

    gatk_cnr_dict = {}
    with open(in_gatk_cnr) as gatk_cnr:
        header = True
        for line in gatk_cnr:
            if header:
                if line[:6] == "CONTIG":
                    header = False
                continue
            columns = line.strip().split("\t")
            chrom = columns[0]
            start_pos = int(columns[1])
            end_pos = int(columns[2])
            log2ratio = float(columns[3])
            if chrom not in gatk_cnr_dict:
                gatk_cnr_dict[chrom] = []
            gatk_cnr_dict[chrom].append([start_pos, end_pos, log2ratio])

    gene_all_dict = {}
    nr_writes = 0
    log.info(f"Opening output tsv file: {output_txt}")
    with open(output_txt, "w") as writer:
        writer.write("gene(s)\tchrom\tregion\tcaller\tfreq_in_db\tcopy_number\tFP_flag")
        out_additional_only.write("gene(s)\tchrom\tregion\tcaller\tfreq_in_db\tcopy_number")
        file1 = True
        for input_org_vcf in input_org_vcfs:
            del_1p19q = {
                "1p_cnvkit": [0, 0], "19q_cnvkit": [0, 0], "1p_gatkcnv": [0, 0], "19q_gatkcnv": [0, 0],
                "1p": [0, 125000000, 125000000], "19q": [26500000, 59128983, 32628983],
            }
            log.info(f"Opening vcf file: {input_org_vcf}")
            variants = VariantFile(input_org_vcf)
            samples = list(variants.header.samples)
            if len(samples) > 1:
                raise Exception(f"Unable to process vcf with more then one sample: {samples}")
            for variant in variants:
                genes = utils.get_annotation_data_info(variant, "Genes")
                log.debug(f"Processing variant: {variant}")
                if isinstance(genes, tuple):
                    genes = ",".join(genes)
                chr = variant.chrom
                start = variant.pos
                end = variant.pos + int(utils.get_annotation_data_info(variant, "SVLEN")) - 1
                size = end - start + 1
                caller = utils.get_annotation_data_info(variant, "CALLER")
                cn = utils.get_annotation_data_info(variant, "CORR_CN")
                AF = utils.get_annotation_data_info(variant, "Twist_AF")
                BAF = utils.get_annotation_data_info(variant, "BAF")
                p_size = 0
                q_size = 0
                if start >= chrom_arm_size[chr][0][0] and start <= chrom_arm_size[chr][0][1]:
                    if end >= chrom_arm_size[chr][0][0] and end <= chrom_arm_size[chr][0][1]:
                        p_size = size
                    else:
                        p_size = chrom_arm_size[chr][0][1] - start + 1
                        q_size = end - chrom_arm_size[chr][0][1]
                else:
                    q_size = size
                if BAF:
                    BAF = float(BAF)
                if AF is None:
                    AF = 0.0
                # 1p19q deletion
                if cn < del_1p19q_cn and chr == "chr1" and start >= del_1p19q["1p"][0] and start <= del_1p19q["1p"][1]:
                    if caller == "cnvkit":
                        del_1p19q["1p_cnvkit"][0] += p_size
                        del_1p19q["1p_cnvkit"][1] += p_size * cn
                    elif caller == "gatk":
                        del_1p19q["1p_gatkcnv"][0] += p_size
                        del_1p19q["1p_gatkcnv"][1] += p_size * cn
                if cn < del_1p19q_cn and chr == "chr19" and start >= del_1p19q["19q"][0] and start <= del_1p19q["19q"][1]:
                    if caller == "cnvkit":
                        del_1p19q["19q_cnvkit"][0] += q_size
                        del_1p19q["19q_cnvkit"][1] += q_size * cn
                    elif caller == "gatk":
                        del_1p19q["19q_gatkcnv"][0] += q_size
                        del_1p19q["19q_gatkcnv"][1] += q_size * cn
                # Large chromosome CNV
                if file1:
                    if caller == "cnvkit":
                        if cn < del_chr_arm_cn_limit:
                            chrom_arm_del[chr][0] += p_size
                            chrom_arm_del[chr][1] += q_size
                        if (
                            cn > normal_cn_lower_limit and cn < normal_cn_upper_limit and
                            BAF and (BAF < normal_baf_lower_limit or BAF > normal_baf_upper_limit)
                        ):
                            chrom_arm_loh[chr][0] += p_size
                            chrom_arm_loh[chr][1] += q_size
                        if cn > amp_chr_arm_cn_limit:
                            chrom_arm_amp[chr][0] += p_size
                            chrom_arm_amp[chr][1] += q_size
                # Baseline check
                if file1:
                    if cn > normal_cn_lower_limit and cn < normal_cn_upper_limit:
                        if caller == "cnvkit":
                            baseline[0] += size
                        else:
                            baseline[1] += size
                # Poliploidy check
                if file1:
                    if caller == "cnvkit":
                        if (
                            cn > normal_cn_lower_limit and cn < normal_cn_upper_limit and
                            BAF and (BAF < normal_baf_lower_limit or BAF > normal_baf_upper_limit)
                        ):
                            polyploidy += size
                        elif (
                            cn > amp_chr_arm_cn_limit and
                            BAF and (BAF > normal_baf_lower_limit and BAF < normal_baf_upper_limit)
                        ):
                            polyploidy += size
                        elif (
                            cn < del_chr_arm_cn_limit and
                            BAF and (BAF > normal_baf_lower_limit and BAF < normal_baf_upper_limit)
                        ):
                            polyploidy += size

                if genes is not None:
                    for gene in genes.split(","):
                        if gene not in gene_all_dict:
                            gene_all_dict[gene] = []
                            gene_all_dict[gene].append([chr, start, end, caller, cn, AF])
                        else:
                            duplicate = False
                            for cnv_variant in gene_all_dict[gene]:
                                if (
                                    chr == cnv_variant[0] and start == cnv_variant[1] and
                                    end == cnv_variant[2] and caller == cnv_variant[3]
                                ):
                                    duplicate = True
                                    break
                            if not duplicate:
                                gene_all_dict[gene].append([chr, start, end, caller, cn, AF])
            fraction_1p_cnvkit = del_1p19q["1p_cnvkit"][0] / del_1p19q["1p"][2]
            fraction_19q_cnvkit = del_1p19q["19q_cnvkit"][0] / del_1p19q["19q"][2]
            if (fraction_1p_cnvkit > del_1p19q_chr_arm_fraction and fraction_19q_cnvkit > del_1p19q_chr_arm_fraction):
                if nr_writes < 2:
                    avg_cn = ((del_1p19q["1p_cnvkit"][1] + del_1p19q["19q_cnvkit"][1]) /
                              (del_1p19q["1p_cnvkit"][0] + del_1p19q["19q_cnvkit"][0]))
                    writer.write(f"\n1p19q\t1p19q\t")
                    writer.write(f"{fraction_1p_cnvkit*100:.0f}%,{fraction_19q_cnvkit*100:.0f}%")
                    writer.write(f"\tcnvkit\tNA\t{avg_cn:.2f}\t-")
                    out_additional_only.write(f"\n1p19q\t1p19q\tNA\tcnvkit\tNA\tNA")
                    nr_writes += 1
            fraction_1p_gatkcnv = del_1p19q["1p_gatkcnv"][0] / del_1p19q["1p"][2]
            fraction_19q_gatkcnv = del_1p19q["19q_gatkcnv"][0] / del_1p19q["19q"][2]
            if (fraction_1p_gatkcnv > del_1p19q_chr_arm_fraction and fraction_19q_gatkcnv > del_1p19q_chr_arm_fraction):
                if nr_writes < 2:
                    avg_cn = ((del_1p19q["1p_gatkcnv"][1] + del_1p19q["19q_gatkcnv"][1]) /
                              (del_1p19q["1p_gatkcnv"][0] + del_1p19q["19q_gatkcnv"][0]))
                    writer.write(f"\n1p19q\t1p19q\t")
                    writer.write(f"{fraction_1p_gatkcnv*100:.0f}%,{fraction_19q_gatkcnv*100:.0f}%")
                    writer.write(f"\tgatk_cnv\tNA\t{avg_cn:.2f}\t-")
                    out_additional_only.write(f"\n1p19q\t1p19q\tNA\tgatk_cnv\tNA\tNA")
                    nr_writes += 1

            file1 = False

        file_nr = 0
        for input_vcf in input_vcfs:
            gene_variant_dict = {}
            log.info(f"Opening vcf file: {input_vcf}")
            variants = VariantFile(input_vcf)
            samples = list(variants.header.samples)
            if file_nr == 1:
                header = variants.header
                out_vcf = VariantFile(out_vcf_filename, "w", header=header)
            if len(samples) > 1:
                raise Exception(f"Unable to process vcf with more then one sample: {samples}")
            counter = 0
            for variant in variants:
                genes = utils.get_annotation_data_info(variant, "Genes")
                log.debug(f"Processing variant: {variant}")
                if isinstance(genes, tuple):
                    genes = ",".join(genes)
                chr = variant.chrom
                start = variant.pos
                end = variant.pos + int(utils.get_annotation_data_info(variant, "SVLEN")) - 1
                caller = utils.get_annotation_data_info(variant, "CALLER")
                cn = utils.get_annotation_data_info(variant, "CORR_CN")
                AF = utils.get_annotation_data_info(variant, "Twist_AF")
                if AF is None:
                    AF = 0.0
                both_callers = False
                for gene in genes.split(","):
                    if gene not in gene_variant_dict:
                        gene_variant_dict[gene] = []
                    gene_variant_dict[gene].append([chr, start, end, caller, cn, AF])
                    if gene in gene_all_dict:
                        nr_callers = {"cnvkit": 0, "gatk": 0, "jumble": 0}
                        for cnv in gene_all_dict[gene]:
                            if cnv[4] > 2.5 or cnv[4] < 1.5:
                                nr_callers[cnv[3]] += 1
                        if nr_callers["cnvkit"] > 0 and nr_callers["gatk"] > 0:
                            both_callers = True
                FP_flag = "-"
                if caller == "cnvkit" and not both_callers:
                    FP_flag = check_fp(chr, start, end, gatk_cnr_dict, cn, max_cnv_fp_size)
                if file_nr == 1:
                    variant.info["FP_FLAG"] = FP_flag
                    out_vcf.write(variant)
                writer.write(f"\n{genes}\t{chr}\t{start}-{end}\t{caller}\t{AF:.2f}\t{cn:.2f}\t{FP_flag}")
                counter += 1

            for gene in gene_variant_dict:
                if len(gene_variant_dict[gene]) < len(gene_all_dict[gene]):
                    org_callers = []
                    for caller_info in gene_variant_dict[gene]:
                        org_callers.append(caller_info[3])
                    for cnv in gene_all_dict[gene]:
                        if cnv[3] not in org_callers:
                            chr = cnv[0]
                            start = cnv[1]
                            end = cnv[2]
                            new_caller = cnv[3]
                            cn = cnv[4]
                            AF = cnv[5]
                            if (
                                (start >= gene_variant_dict[gene][0][1] and start <= gene_variant_dict[gene][0][2]) or
                                (end >= gene_variant_dict[gene][0][1] and end <= gene_variant_dict[gene][0][2]) or
                                (gene_variant_dict[gene][0][1] >= start and gene_variant_dict[gene][0][1] <= end) or
                                (gene_variant_dict[gene][0][2] >= end and gene_variant_dict[gene][0][2] <= start)
                            ):
                                writer.write(f"\n{gene}\t{chr}\t{start}-{end}\t{new_caller}\t{AF:.2f}\t{cn:.2f}\t-")

            if file_nr == 1:
                out_vcf.close()
            file_nr += 1
        log.info(f"Processed {counter} variants")

        deletions = open(input_del)
        header_list = next(deletions).split("\t")
        for deletion in deletions:
            columns = {k: v for k, v in zip(header_list, deletion.strip().split("\t"))}
            gene = columns['Gene(s)']
            chr = columns['Chromosome']
            start = columns['Gene_start']
            end = columns['Gene_end']
            caller = "small_deletion"
            AF = "NA"
            log_odds_ratio = float(columns['Median_L2R_deletion'])
            cn = 2*pow(2, float(log_odds_ratio))
            ccn = cn
            if TC > 0.0:
                ccn = round(2 + (cn - 2) * (1/float(TC)), 2)
            writer.write(f"\n{gene}\t{chr}\t{start}-{end}\t{caller}\t{AF}\t{ccn:.2f}\t-")
            out_additional_only.write(f"\n{gene}\t{chr}\t{start}-{end}\t{caller}\t{AF}\t{ccn:.2f}")
        amplifications = open(input_amp)
        header_list = next(amplifications).split("\t")
        for amplification in amplifications:
            columns = {k: v for k, v in zip(header_list, amplification.strip().split("\t"))}
            gene = columns['Gene(s)']
            chr = columns['Chromosome']
            start = columns['Gene_start']
            end = columns['Gene_end']
            caller = "small_amplification"
            AF = "NA"
            log_odds_ratio = float(columns['Median_L2R_amplification'])
            cn = 2*pow(2, float(log_odds_ratio))
            ccn = cn
            if TC > 0.0:
                ccn = round(2 + (cn - 2) * (1/float(TC)), 2)
            if ccn > amp_cn_limit:
                writer.write(f"\n{gene}\t{chr}\t{start}-{end}\t{caller}\t{AF}\t{ccn:.2f}\t-")
                out_additional_only.write(f"\n{gene}\t{chr}\t{start}-{end}\t{caller}\t{AF}\t{ccn:.2f}")

    with open(out_tsv_chrom_arms, "w") as writer:
        writer.write("chrom\tarm\tcaller\ttype\tfraction")
        if baseline[0] / genome_size < baseline_fraction_limit:
            writer.write(f"\nWarning: baseline of GATK CNV might be shifted!\t\t\t\t")
            writer.write(f"{baseline[0] * 100 / genome_size:.1f}% on baseline")
        if baseline[1] / genome_size < baseline_fraction_limit:
            writer.write(f"\nWarning: baseline of CNVkit might be shifted!\t\t\t\t")
            writer.write(f"{baseline[1] * 100 / genome_size:.1f}% on baseline")
        if polyploidy / genome_size > polyploidy_fraction_limit:
            writer.write(f"\nWarning: potential polyploidy detected!\t\t\t\t")
            writer.write(f"{polyploidy * 100 / genome_size:.1f}% polyploid regions")
        for chrom in chrom_arm_del:
            if chrom_arm_del[chrom][0] / chrom_arm_size[chrom][0][2] > chr_arm_fraction:
                writer.write(f"\n{chrom}\tp\tcnvkit\tdeletion\t")
                writer.write(f"{chrom_arm_del[chrom][0] * 100 / chrom_arm_size[chrom][0][2]:.1f}%")
            if chrom_arm_amp[chrom][0] / chrom_arm_size[chrom][0][2] > chr_arm_fraction:
                writer.write(f"\n{chrom}\tp\tcnvkit\tduplication\t")
                writer.write(f"{chrom_arm_amp[chrom][0] * 100 / chrom_arm_size[chrom][0][2]:.1f}%")
            if chrom_arm_loh[chrom][0] / chrom_arm_size[chrom][0][2] > chr_arm_fraction:
                writer.write(f"\n{chrom}\tp\tcnvkit\tloh\t{chrom_arm_loh[chrom][0] * 100 / chrom_arm_size[chrom][0][2]:.1f}%")
            if chrom_arm_del[chrom][1] / chrom_arm_size[chrom][1][2] > chr_arm_fraction:
                writer.write(f"\n{chrom}\tq\tcnvkit\tdeletion\t")
                writer.write(f"{chrom_arm_del[chrom][1] * 100 / chrom_arm_size[chrom][1][2]:.1f}%")
            if chrom_arm_amp[chrom][1] / chrom_arm_size[chrom][1][2] > chr_arm_fraction:
                writer.write(f"\n{chrom}\tq\tcnvkit\tduplication\t")
                writer.write(f"{chrom_arm_amp[chrom][1] * 100 / chrom_arm_size[chrom][1][2]:.1f}%")
            if chrom_arm_loh[chrom][1] / chrom_arm_size[chrom][1][2] > chr_arm_fraction:
                writer.write(f"\n{chrom}\tq\tcnvkit\tloh\t{chrom_arm_loh[chrom][1] * 100 / chrom_arm_size[chrom][1][2]:.1f}%")


if __name__ == "__main__":
    in_vcfs = snakemake.input.vcfs
    in_org_vcfs = snakemake.input.org_vcfs
    in_del = snakemake.input.deletions
    in_amp = snakemake.input.amplifications
    in_chrom_arm_size = snakemake.input.chrom_arm_size
    in_gatk_cnr = snakemake.input.gatk_cnr
    amp_cn_limit = snakemake.params.call_small_amplifications_cn_limit
    out_tsv = snakemake.output.tsv
    out_tsv_chrom_arms = snakemake.output.tsv_chrom_arms
    out_vcf_filename = snakemake.output.vcf_del
    del_1p19q_cn = snakemake.params.del_1p19q_cn_limit
    del_1p19q_chr_arm_fraction = snakemake.params.del_1p19q_chr_arm_fraction
    chr_arm_fraction = snakemake.params.chr_arm_fraction
    del_chr_arm_cn_limit = snakemake.params.del_chr_arm_cn_limit
    amp_chr_arm_cn_limit = snakemake.params.amp_chr_arm_cn_limit
    normal_cn_lower_limit = snakemake.params.normal_cn_lower_limit
    normal_cn_upper_limit = snakemake.params.normal_cn_upper_limit
    normal_baf_lower_limit = snakemake.params.normal_baf_lower_limit
    normal_baf_upper_limit = snakemake.params.normal_baf_upper_limit
    baseline_fraction_limit = snakemake.params.baseline_fraction_limit
    polyploidy_fraction_limit = snakemake.params.polyploidy_fraction_limit
    max_cnv_fp_size = snakemake.params.max_cnv_fp_size
    TC = float(snakemake.params.tc)
    with open(snakemake.output.tsv_additional_only, "w") as out_additional_only:
        create_tsv_report(
            in_vcfs,
            in_org_vcfs,
            in_del,
            in_amp,
            in_chrom_arm_size,
            in_gatk_cnr,
            amp_cn_limit,
            out_tsv,
            out_additional_only,
            out_tsv_chrom_arms,
            out_vcf_filename,
            del_1p19q_cn,
            del_1p19q_chr_arm_fraction,
            chr_arm_fraction,
            del_chr_arm_cn_limit,
            amp_chr_arm_cn_limit,
            normal_cn_lower_limit,
            normal_cn_upper_limit,
            normal_baf_lower_limit,
            normal_baf_upper_limit,
            baseline_fraction_limit,
            polyploidy_fraction_limit,
            TC,
            max_cnv_fp_size,
        )
