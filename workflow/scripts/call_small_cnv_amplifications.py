
import statistics
import logging

log = logging.getLogger()


def read_cnv_data(cnv_file_name, sample_name, region):
    probe_data = []
    probe_positions = []
    gene_probe_index = []
    with open(cnv_file_name) as cnv_file:
        header = True
        i = 0
        for line in cnv_file:
            if header:
                if line.find("CONTIG") != -1:
                    header = False
                continue
            columns = line.strip().split("\t")
            chrom = columns[0]
            start = int(columns[1])
            end = int(columns[2])
            log2_copy_ratio = float(columns[3])
            if chrom == region[0] and ((start >= region[1] and start <= region[2]) or (end >= region[1] and end <= region[2])):
                probe_data.append(log2_copy_ratio)
                probe_positions.append([chrom, start, end])
                if (
                    chrom == region[0] and
                    ((start >= region[3] and start <= region[4]) or (end >= region[3] and end <= region[4]))
                ):
                    gene_probe_index.append(i)
                i += 1
    return probe_data, gene_probe_index, probe_positions


def find_max_probe_diff(probe_data, window_size):
    probe_len = len(probe_data)
    i = 0
    j = i + 1 + window_size
    probe_windows1 = []
    probe_windows2 = []
    while j + window_size + 1 < probe_len:
        probe_windows1.append([sum(probe_data[i:i + window_size]) / window_size, i, i + window_size])
        probe_windows2.append([sum(probe_data[j:j + window_size]) / window_size, j, j + window_size])
        i += 1
        j += 1
    k = 0
    max_probe_window_diff = 0.0
    min_probe_window_diff = 0.0
    k_max_index = 0
    k_min_index = 0
    for probe_window in probe_windows1:
        probe_window_diff = probe_window[0] - probe_windows2[k][0]
        if probe_window_diff > max_probe_window_diff:
            max_probe_window_diff = probe_window_diff
            k_max_index = k
        if probe_window_diff < min_probe_window_diff:
            min_probe_window_diff = probe_window_diff
            k_min_index = k
        k += 1
    return k_max_index, k_min_index


def filter_amplifications(
        max_probe_diff_index, min_probe_diff_index, probe_data, gene_probe_index, probe_positions, region,
        amplifications, region_max_size, window_size, min_log_odds_diff, min_nr_stdev_diff
):
    # Filter deletions
    if max_probe_diff_index <= min_probe_diff_index:
        return "Deletion"
    # Filter large amplifications (captured by other tools)
    if max_probe_diff_index - min_probe_diff_index >= region_max_size:
        return "Too_large"
    # Filter amplifications not in the actual gene of interest
    in_gene = False
    # Adding one extra to allow the breakpoint to not be included
    k = min_probe_diff_index + window_size + 1
    while k < max_probe_diff_index + window_size:
        if k in gene_probe_index:
            in_gene = True
            break
        k += 1
    if not in_gene:
        return "Not_in_gene"
    # Filter too short amplifications
    high_probes_i_start = min_probe_diff_index + window_size + 1
    high_probes_i_stop = max_probe_diff_index + window_size
    high_probes = probe_data[high_probes_i_start:high_probes_i_stop]
    if len(high_probes) < window_size:
        return "Too_small"
    # Calculate low probes window averages
    low_probes_window_averages = []
    k = 0
    while k + window_size < min_probe_diff_index:
        low_probes_window_averages.append(sum(probe_data[k:k+window_size]) / window_size)
        k += 1
    k = max_probe_diff_index + len(high_probes) + 1
    while k + window_size < len(probe_data):
        low_probes_window_averages.append(sum(probe_data[k:k+window_size]) / window_size)
        k += 1
    # Filter amplifications with too few data points outside to calculate median and standard deviations
    if len(low_probes_window_averages) < 4:
        return "Too_few_outside"
    # Calculate high and low medians and stdev for high probes window averages
    median_high = statistics.median(high_probes)
    median_low = statistics.median(low_probes_window_averages)
    stdev_low = statistics.stdev(low_probes_window_averages)
    if stdev_low == 0:
        stdev_low = 0.01
    # Amplifications must have positive median log odds ration
    if median_high < 0:
        return "Not_amplification"
    # Filter based on difference in log odds ration between high and low median (min_log_odds_diff)
    if median_high - median_low < min_log_odds_diff:
        return "Low_abs_diff"
    # Filter based on the number of standard deviation between the high and low median (min_nr_stdev_diff)
    if median_high < median_low + stdev_low * min_nr_stdev_diff:
        return "Low_nr_std_diff"
    median_diff = median_high - median_low
    nr_std_diff = abs(median_diff) / stdev_low
    start_pos = probe_positions[high_probes_i_start][1]
    end_pos = probe_positions[high_probes_i_stop - 1][2]
    amplifications.write(f"{region[5]}\t{region[0]}\t{start_pos}\t{end_pos}")
    amplifications.write(f"\t{median_diff}\t{median_high}\t{median_low}\t{len(high_probes)}\t{nr_std_diff}\n")
    return "Unfiltered"


def read_regions_data(regions_file):
    next(regions_file)
    regions = []
    for line in regions_file:
        columns = line.strip().split("\t")
        gene_name = columns[0]
        chrom = columns[1]
        start_cnv_region = int(columns[2])
        stop_cnv_region = int(columns[3])
        gene_start = int(columns[4])
        gene_stop = int(columns[5])
        regions.append([chrom, start_cnv_region, stop_cnv_region, gene_start, gene_stop, gene_name])
    return regions


def call_small_cnv_amplifications(
    cnv_file_name, regions_file, amplifications, window_size, region_max_size, min_nr_stdev_diff, min_log_odds_diff,
):
    regions = read_regions_data(regions_file)
    sample_name = cnv_file_name.split("/")[1].split("_")[0]
    amplifications.write("Gene(s)\tChromosome\tGene_start\tGene_end\tLog2_ratio_diff\tMedian_L2R_amplification\t")
    amplifications.write("Median_L2R_surrounding\tNumber_of_data_points\tNumber_of_stdev\n")
    for region in regions:
        probe_data, gene_probe_index, probe_positions = read_cnv_data(cnv_file_name, sample_name, region)
        # Warning about to small region
        # 1 window for Amplifications and one window plus 3 for high probes (for statistics) plus 2 spaces between windows)
        if len(probe_data) < window_size * 2 + 3 + 2:
            log.info(f"Too few data points for region: {region}")
        max_probe_diff_index, min_probe_diff_index = find_max_probe_diff(probe_data, window_size)
        filter = filter_amplifications(
            max_probe_diff_index, min_probe_diff_index, probe_data, gene_probe_index, probe_positions, region,
            amplifications, region_max_size, window_size, min_log_odds_diff, min_nr_stdev_diff,
        )
    return filter


if __name__ == "__main__":
    log = snakemake.log_fmt_shell(stdout=False, stderr=True)

    call_small_cnv_amplifications(
        snakemake.input.cnv_data,
        open(snakemake.input.regions_file),
        open(snakemake.output.amplifications, "w"),
        snakemake.params.window_size,
        snakemake.params.region_max_size,
        snakemake.params.min_nr_stdev_diff,
        snakemake.params.min_log_odds_diff,
    )
