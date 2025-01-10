
__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2024, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"

import pysam
import statistics
import numpy as np
import scipy.stats as stats


def read_segments(input_segments):
    segments = pysam.VariantFile(input_segments)

    segment_dict = {}
    for segment in segments.fetch():
        chrom = segment.chrom
        start_pos = segment.pos
        end_pos = segment.stop
        CN = segment.info["CORR_CN"]
        if chrom.endswith("X") or chrom.endswith("Y"):
            continue
        if chrom not in segment_dict:
            segment_dict[chrom] = []
        segment_dict[chrom].append([start_pos, end_pos, CN, []])
    return segment_dict


def read_germline_vcf(input_germline_vcf, segment_dict, min_germline_af):
    germline_vcf = pysam.VariantFile(input_germline_vcf)

    for record in germline_vcf.fetch():
        chrom = record.chrom
        pos = record.pos
        AF = record.info["AF"][0]
        if AF < min_germline_af or AF > (1 - min_germline_af):
            continue
        if chrom in segment_dict:
            i = 0
            for segment in segment_dict[chrom]:
                if pos >= segment[0] and pos <= segment[1]:
                    segment_dict[chrom][i][3].append(AF)
                i += 1
    return segment_dict


def read_snv_vcf_and_find_max_af(input_snv_vcf, segment_dict, max_somatic_af, gnomAD_AF_limit):
    snv_vcf = pysam.VariantFile(input_snv_vcf)
    snv_list = []

    # VEP annotation
    vep_fields = {}
    for record in snv_vcf.header.records:
        if record.type == "INFO":
            if record['ID'] == "CSQ":
                vep_fields = {v: c for c, v in enumerate(record['Description'].split("Format: ")[1].split('">')[0].split("|"))}

    for record in snv_vcf.fetch():
        chrom = record.chrom
        pos = record.pos
        ref = record.ref
        alt = record.alts[0]
        # Only keep variants called by both callers
        try:
            CALLERS = record.info["CALLERS"]
        except KeyError:
            continue
        if "vardict" not in CALLERS or "gatk_mutect2" not in CALLERS:
            continue
        # Only keep SNVs
        if len(ref) > 1 or len(alt) > 1:
            continue
        # Check if SNVs is in segment with clear CNA based on copy number and then skip the SNV
        CNA = False
        if chrom not in segment_dict:
            continue
        for seg in segment_dict[chrom]:
            if pos >= seg[0] and pos <= seg[1]:
                if seg[2] < 1.6 or seg[2] > 2.4:
                    CNA = True
                    break
        if CNA:
            continue
        # Skip low AF somatic SNVs and germline SNPs
        gnomAD_AF = record.info["CSQ"][0].split("|")[vep_fields["gnomAD_AF"]]
        if gnomAD_AF == "":
            gnomAD_AF = 0
        else:
            gnomAD_AF = float(gnomAD_AF)
        if gnomAD_AF > gnomAD_AF_limit:
            continue
        AF = record.samples.items()[0][1]["AF"][0]
        # Skip low AF somatic SNVs
        if AF > max_somatic_af:
            continue
        # Only keep variants in exons
        consequence = record.info["CSQ"][0].split("|")[vep_fields["Consequence"]].split("&")
        if not (
                "missense_variant" in consequence or
                "synonymous_variant" in consequence or
                "stop_lost" in consequence or
                "stop_gained" in consequence or
                "stop_retained_variant" in consequence
        ):
            continue

        snv_list.append(AF)

    if len(snv_list) > 0:
        return 2 * max(snv_list)
    else:
        return 0


def baf_to_tc(abs_value_seg_median, CN, CN_list, median_noise_level):
    # Remove the medium noise level
    abs_value_seg_median -= median_noise_level
    # Get highest and lowest difference to normal copy number for a segment
    min_CN_diff = 0.00001
    max_CN_diff = 0.00001
    if len(CN_list) > 0:
        min_CN_diff = max(0.00001, 2.0 - min(CN_list))
        max_CN_diff = max(0.00001, max(CN_list) - 2.0)
    # If clear aneuploidy (only based on copy number)
    if len(CN_list) > 0 and (min(CN_list) < 1.85 or max(CN_list) > 2.15):
        # Deletion if copy number of at least 40% of the lowest copy number segment
        if (2.0 - CN) / min_CN_diff > 0.4:
            return [(2 * abs_value_seg_median)/(abs_value_seg_median + 0.5), "Del"]
        # Duplication (not reported) if copy number of at least 40% of the highest copy number segment
        elif (CN - 2.0) / max_CN_diff > 0.4:
            return [(2 * abs_value_seg_median)/(abs_value_seg_median + 0.333), "Dup"]
        # Copy neutral LoH (Reported if no deletion) if copy number looks normal
        elif (2.0 - CN) / min_CN_diff < 0.4 and (CN - 2.0) / max_CN_diff < 0.4:
            return [(2 * abs_value_seg_median)/(abs_value_seg_median + 0.667), "CNLoH"]
        # Unknown, assuming deletion (not reported)
        else:
            return [(2 * abs_value_seg_median)/(abs_value_seg_median + 0.5), "Unknown"]
    # If no clear CNA signal, assume deletion for TC adjustment model for all segments
    else:
        return [(2 * abs_value_seg_median)/(abs_value_seg_median + 0.5), "Del"]


def test_if_signal_in_segment(segment, abs_value_seg_median, median_noise_level, vaf_baseline):
    # Get VAF for max density on the lower AF half of the BAF plot
    kde1 = stats.gaussian_kde(segment)
    kde1.set_bandwidth(bw_method=kde1.factor / 2.)
    x = np.linspace(0, vaf_baseline, 200)
    kde1 = kde1(x)
    max1 = x[kde1.argmax()]

    # Get VAF for max density on the upper AF half of the BAF plot
    kde2 = stats.gaussian_kde(segment)
    kde2.set_bandwidth(bw_method=kde2.factor / 2.)
    x = np.linspace(vaf_baseline, 1, 200)
    kde2 = kde2(x)
    max2 = x[kde2.argmax()]

    # Get VAF for min density between the two maxima
    kde3 = stats.gaussian_kde(segment)
    kde3.set_bandwidth(bw_method=kde3.factor / 2.)
    x = np.linspace(max1, max2, 200)
    kde3 = kde3(x)
    minimum = x[kde3.argmin()]

    # Get the density for the maxima and minimum
    density_max1 = kde1[kde1.argmax()]
    density_max2 = kde2[kde2.argmax()]
    density_min = kde3[kde3.argmin()]
    '''
    Return True (the segment has CNA) if all of these are true:
        * the signal is higher than the noise
        * the minimum VAF is not the same as any maxima
        * the minimum VAF is between 0.4 and 0.6
        * the maxima ratio is between 0.5 and 2
        * the minimum maximum ration is lower than 0.85
    '''
    if abs_value_seg_median > median_noise_level:
        if minimum != max1 and minimum != max2 and minimum > 0.4 and minimum < 0.6:
            if (
                density_max1 / density_max2 < 2 and
                density_max1 / density_max2 > 0.5 and
                density_min / density_max1 < 0.85 and
                density_min / density_max2 < 0.85
            ):
                return True
    return False


def calculate_cnv_tc(segment_dict_AF, min_nr_SNPs_per_segment, vaf_baseline, min_segment_length):
    # Calculate BAF noise levels for segments without CNA (noise_level)
    # Save copy numbers for segment with CNA signal (CN_signal_list)
    CN_signal_list = []
    noise_level = []
    for chrom in segment_dict_AF:
        for segment in segment_dict_AF[chrom]:
            if len(segment[3]) > min_nr_SNPs_per_segment:
                if not test_if_signal_in_segment(segment[3], 1, 0, vaf_baseline):
                    for AF in segment[3]:
                        noise_level.append(abs(AF-vaf_baseline))
                else:
                    CN_signal_list.append(segment[2])
    median_noise_level = statistics.median(noise_level)

    # Iterate through segments to find final CNAs
    tc_dict = {"Del": [], "Dup": [], "CNLoH": [], "Unknown": []}
    for chrom in segment_dict_AF:
        for segment in segment_dict_AF[chrom]:
            data_points_seg = len(segment[3])
            segment_length = segment[1] - segment[0]
            # Only use segments that have sufficient number of SNPs and segment length (in bases)
            if data_points_seg > min_nr_SNPs_per_segment and segment_length > min_segment_length:
                i = 0
                short = True
                # For extra long segments, split them up into 100 datapoint segments to find smaller CNAs
                while i + min_nr_SNPs_per_segment < data_points_seg:
                    seg = segment[3]
                    if data_points_seg > 100:
                        seg = segment[3][i:i+100]
                        short = False
                    abs_value_seg = []
                    for AF in seg:
                        abs_value_seg.append(abs(AF-vaf_baseline))
                    abs_value_seg_median = statistics.median(abs_value_seg)
                    # If signal is found calculate TC based on the median separation in BAF around the BAF-baseline (~50%)
                    if test_if_signal_in_segment(seg, abs_value_seg_median, median_noise_level, vaf_baseline):
                        tc_seg = baf_to_tc(abs_value_seg_median, segment[2], CN_signal_list, median_noise_level)
                        tc_dict[tc_seg[1]].append([tc_seg[0], segment, chrom])
                    i += 50
                    if short:
                        break
    # Report highest TC based on CNAs with deletions. If none are found report highest copy neutral LoH CNA.
    max_tc = 0
    found_del = False
    seg_list = []
    for seg_info in tc_dict["Del"]:
        seg_list.append([seg_info, "Deletion"])
        if seg_info[0] > max_tc:
            max_tc = seg_info[0]
            found_del = True
    for seg_info in tc_dict["CNLoH"]:
        seg_list.append([seg_info, "CNLoH"])
        if found_del and seg_info[0] > max_tc:
            max_tc = seg_info[0]
    return max_tc, seg_list


def write_tc(output_tc, tc_cnv, tc_snv):
    output = open(output_tc, "w")
    output.write("Percentage ctDNA based on CNV data\tPercentage ctDNA based on SNV data\n")
    output.write(f"{tc_cnv*100:.1f}%\t{tc_snv*100:.1f}%\n")
    output.close()
    return f"{tc_cnv*100:.1f}%\t{tc_snv*100:.1f}%\n"


def write_seg_list(output_file_name, seg_list):
    output = open(output_file_name, "w")
    output.write("ctDNA_fraction\tCNV_type\tchromosome\tstart_pos\tend_pos\n")
    for seg in seg_list:
        output.write(f"{seg[0][0]*100:.1f}%\t{seg[1]}\t{seg[0][2]}\t{seg[0][1][0]}\t{seg[0][1][1]}\n")


if __name__ == "__main__":
    input_segments = snakemake.input.segments
    input_germline_vcf = snakemake.input.germline_vcf
    input_vcf = snakemake.input.vcf
    output_ctDNA_fraction = snakemake.output.ctDNA_fraction
    output_ctDNA_info = snakemake.output.ctDNA_info

    min_germline_af = float(snakemake.params.min_germline_af)
    max_somatic_af = float(snakemake.params.max_somatic_af)
    min_nr_SNPs_per_segment = int(snakemake.params.min_nr_SNPs_per_segment)
    min_segment_length = int(snakemake.params.min_segment_length)
    gnomAD_AF_limit = float(snakemake.params.gnomAD_AF_limit)
    vaf_baseline = float(snakemake.params.vaf_baseline)

    # Read CNV segments
    segment_dict = read_segments(input_segments)
    # Read germline SNPs
    segment_dict_AF = read_germline_vcf(input_germline_vcf, segment_dict, min_germline_af)
    tc_cnv, seg_list = calculate_cnv_tc(segment_dict_AF, min_nr_SNPs_per_segment, vaf_baseline, min_segment_length)
    # Read SNVs and report TC based on max VAF of somatic SNV.
    tc_snv = read_snv_vcf_and_find_max_af(input_vcf, segment_dict, max_somatic_af, gnomAD_AF_limit)
    write_tc(output_ctDNA_fraction, tc_cnv, tc_snv)
    # Write additional info regarding which chromosomes have deletions
    write_seg_list(output_ctDNA_info, seg_list)
