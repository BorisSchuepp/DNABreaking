# Filtering algorithm for the analysis of the breaking pattern
# of (nicked) dsDNA from NGS data
# Author: Boris N. Schuepp
# Version: 1.1

import statistics
import time
from datetime import datetime
import ConfigurationHandeler
from PlottingNGS import *


# Translate Phred quality characters to numerical values. Provide the quality characters as string,
# returns list of integer quality values
def translate_quality_values(quality_string_in):
    quality_dict = {"!": 0, "\"": 1, "#": 2, "$": 3, "%": 4, "&": 5, "'": 6, "(": 7, ")": 8, "*": 9, "+": 10, ",": 11,
                    "-": 12, ".": 13, "/": 14, "0": 15, "1": 16, "2": 17, "3": 18, "4": 19, "5": 20, "6": 21, "7": 22,
                    "8": 23, "9": 24, ":": 25, ";": 26, "<": 27, "=": 28, ">": 29, "?": 30, "@": 31, "A": 32, "B": 33,
                    "C": 34, "D": 35, "E": 36, "F": 37, "G": 38, "H": 39, "I": 40}

    quality_numbers_list_out = []
    for quality_char in quality_string_in:
        quality_numbers_list_out.append(quality_dict[quality_char])

    return quality_numbers_list_out


# Calculate the quality average for every base distance from the 5'-end, provide a list of numerical
# quality values, returns a dictionary[distance from 5'-end] = average quality value
def calculate_quality_average(quality_values_in):
    longest_read = max([len(i) for i in quality_values_in])
    quality_dict_out = {}
    for i in range(0, longest_read):
        current_count = 0
        current_quality_sum = 0
        for quality_value in quality_values_in:
            if len(quality_value) > i:
                current_quality_sum += quality_value[i]
                current_count += 1
        quality_dict_out[i] = current_quality_sum / current_count

    return quality_dict_out


# Read an input file in fastq format and return a dictionary with keys = sequence ids
# and values = [read sequence, quality values]
def read_input_file(input_file_name_in):
    try:
        input_file = open(input_file_name_in, "r")
    except FileNotFoundError:
        print(f"ERROR: File {input_file_name_in} not found! Please check the file and the configuration!. Aborting.")
        exit(3)

    # Fastq file format: Sequence id, sequence, not relevant field, quality data
    raw_data_out = {}
    running_count = 0
    current_read = []
    for line in input_file:
        line = line.strip("\n")
        if running_count != 2:
            current_read.append(line)
        if running_count == 3:
            index = current_read[0].split()[0]
            if index in raw_data_out.keys():
                print(f"WARNING: Multiple entries for same index {index}!")
            raw_data_out[index] = [current_read[1], current_read[2]]
            current_read = []
        running_count = (running_count + 1) % 4

    input_file.close()
    return raw_data_out


# Match single reads to create paired reads using their sequence ids, print warnings if any single
# read stays unmatched.
def match_paired_reads(reads_1_in, reads_2_in):
    paired_reads_out = []
    for index in reads_1_in.keys():
        if index in reads_2_in.keys():
            paired_reads_out.append(
                (reads_1_in[index], reads_2_in[index]))
        else:
            print(f"WARNING: Index {index} wasn't found in both files!")

    for index in reads_2_in.keys():
        if index not in reads_1_in.keys():
            print(f"WARNING: Index {index} wasn't found in both files!")

    return paired_reads_out


# Provided a paired read and a reference sequence as well as filtering parameters both reads will be matched
# against the reference. Two list are returned one is the list of all assigments (even those that have not been
# sucessful) and the other is the list of only valid paired reads.
def match_to_reference(paired_reads_in, reference_in, no_mismatch_tolerance_distance_in, maxium_further_mismatches_in):
    strand_a_53 = reference_in
    strand_b_53 = get_complementary_strand(reference_in)
    all_matching_results = []
    filtered_matching_results = []

    for paired_read in paired_reads_in:
        strand_read_1, matching_position_read_1, number_of_mismatches_read_1 = match_to_reference_single_read(
            paired_read[0], strand_a_53, strand_b_53, no_mismatch_tolerance_distance_in, maxium_further_mismatches_in)
        strand_read_2, matching_position_read_2, number_of_mismatches_read_2 = match_to_reference_single_read(
            paired_read[1], strand_a_53, strand_b_53, no_mismatch_tolerance_distance_in, maxium_further_mismatches_in)
        # Collect a list of all filtering results (even those that will be ignored)
        all_matching_results.append((
            strand_read_1, matching_position_read_1, number_of_mismatches_read_1,
            strand_read_2, matching_position_read_2, number_of_mismatches_read_2))
        # Both reads are found and are on different strands of the template dsDNS
        if strand_read_1 >= 0 and strand_read_2 >= 0 and strand_read_1 + strand_read_2 == 1:
            filtered_matching_results.append((
                strand_read_1, matching_position_read_1, number_of_mismatches_read_1,
                strand_read_2, matching_position_read_2, number_of_mismatches_read_2))

    return all_matching_results, filtered_matching_results


# Match a single read to the template strand or its complement using a region (no_tolerance_distance) that has
# to exactly fit some part of the template (or complement). After the position has been determined count how
# many of the remaining bases don't match the expected base from the template (or its complement).
# Return three values: strand on which the read was found (1 = strandA , 2 = strandB , -1 = wasn't found),
# Base index on template (strandA) where the read starts (-1 if read was not found),
# Number of mismatches found in the remaining part of the read (-1 if read was not found).
def match_to_reference_single_read(read_in, reference_strand_a, reference_strand_b,
                                   no_mismatch_tolerance_distance_in, maxium_further_mismatches_in):
    initial_bases = read_in[0][:no_mismatch_tolerance_distance_in]
    pos_strand_a = reference_strand_a.find(initial_bases)
    pos_strand_b = reference_strand_b.find(initial_bases)
    # The inital bases are not found on the reference or its complement
    if pos_strand_a == -1 and pos_strand_b == -1:
        return -1, -1, -1

    if pos_strand_a != -1:
        strand_found = reference_strand_a
        pos_found = pos_strand_a
        found_index = 0
        result_pos = pos_strand_a + 1
    else:
        strand_found = reference_strand_b
        pos_found = pos_strand_b
        found_index = 1
        result_pos = len(reference_strand_b) - pos_strand_b

    remaining_sequence = read_in[0][no_mismatch_tolerance_distance_in:]
    reference_sequence_stretch = strand_found[pos_found + no_mismatch_tolerance_distance_in:
                                              pos_found + no_mismatch_tolerance_distance_in + len(remaining_sequence)]

    mismatch_indices = get_mismatches(remaining_sequence, reference_sequence_stretch)

    if len(mismatch_indices) > maxium_further_mismatches_in:
        return -2, -1, len(mismatch_indices)
    else:
        return found_index, result_pos, len(mismatch_indices)

# Tries to match a read on a single reference strand. Firstly, the inital bases are matched without tolerance for mismatches
# Afterwards the remainder is matched with a number of mismatches that are tolerated. 
def match_to_reference_single_read_single_sequence(read_in, reference_strand,
                                   no_mismatch_tolerance_distance_in, maxium_further_mismatches_in):


    # Check if the initial bases are found on the reference 
    initial_bases = read_in[0][:no_mismatch_tolerance_distance_in]
    pos_strand = reference_strand.find(initial_bases)
    if pos_strand == -1:
        return -1, -1
    else:
        pos_found = pos_strand
        result_pos = pos_strand + 1

    # Count the missmatches on the remainder 
    remaining_sequence = read_in[0][no_mismatch_tolerance_distance_in:]
    reference_sequence_stretch = reference_strand[pos_found + no_mismatch_tolerance_distance_in:
                                              pos_found + no_mismatch_tolerance_distance_in + len(remaining_sequence)]

    mismatch_indices = get_mismatches(remaining_sequence, reference_sequence_stretch)

    if len(mismatch_indices) > maxium_further_mismatches_in:
        return -1, len(mismatch_indices)
    else:
        return result_pos, len(mismatch_indices)


# Provide the remaing part of the read and the expected sequence as determined from the inital positional match
# returns a list of indices in which the read does not match the expected sequence
def get_mismatches(remaining_read_in, reference_sequence_stretch_in):
    index = 0
    mismatch_indices = []
    for base_A, base_B in zip(remaining_read_in, reference_sequence_stretch_in):
        if base_A != base_B:
            mismatch_indices.append(index)
        index += 1
    return mismatch_indices


# Returns the first and the last base index for each read based on the 5'-end positions measured on
# the upper strand (the one provided in 5'->3' direction in the configuration). Due to the orientation
# of the upper strand the first base index is the 5'-end on the upper strand while the last base index
# is the 5'-end on the lower strand.
def fragment_start_end(filtered_paired_reads_in):
    fragments_start_end = []
    for paired_read in filtered_paired_reads_in:
        start_a = paired_read[1]
        start_b = paired_read[4]
        left_start = min(start_a, start_b)
        right_start = max(start_a, start_b)
        fragments_start_end.append((left_start, right_start))
    return fragments_start_end


# Returns the complementary DNA strand for the provided input strand
def get_complementary_strand(strand_in):
    complement_35 = ""
    for base in strand_in:
        if base == "A":
            complement_35 += "T"
        if base == "T":
            complement_35 += "A"
        if base == "G":
            complement_35 += "C"
        if base == "C":
            complement_35 += "G"
    complement_53 = complement_35[::-1]
    return complement_53


# Used to bin the results from the filtering into some predefined categories, returns a list of
# bin descriptors and fractional bin sizes in terms of fraction of total single reads
def classify_filtering(mismatch_results_in):
    mismatch_keys = ["0 mismatches", "1 mismatch", "2 mismatches", "3 to 10 mismatches",
                     "5'-end not matching reference", "More than 10 mismatches"]

    classified_results_dict = {key: 0 for key in mismatch_keys}
    for result in mismatch_results_in:
        res = [result[0], result[3]]
        mis = [result[2], result[5]]
        for res_i, mis_i in zip(res, mis):
            if res_i == -2:
                classified_results_dict["More than 10 mismatches"] += 1
            if res_i == -1:
                classified_results_dict["5'-end not matching reference"] += 1
            if res_i == 0 or res_i == 1:
                if mis_i == 0:
                    classified_results_dict["0 mismatches"] += 1
                if mis_i == 1:
                    classified_results_dict["1 mismatch"] += 1
                if mis_i == 2:
                    classified_results_dict["2 mismatches"] += 1
                if 3 <= mis_i <= 10:
                    classified_results_dict["3 to 10 mismatches"] += 1

    classified_results_list = [classified_results_dict[key] for key in mismatch_keys]
    classified_results_dict = {i: classified_results_dict[i] / sum(classified_results_list)
                               for i in classified_results_dict.keys()}
    return classified_results_dict


# Use reassembly paired reads to obtain a framgent length distribution. Provide a list of (i,j) position, returns
# a dictionary[fragment length] = fractional of reads with this length
def get_length_dict(fragment_ends_in):
    length_dict = {}
    lengths_all_out = []
    for fragment in fragment_ends_in:
        current_length = fragment[1] - fragment[0] + 1
        lengths_all_out.append(current_length)
        if current_length in length_dict.keys():
            length_dict[current_length] += 1
        else:
            length_dict[current_length] = 1
    mean_length_out = round(statistics.mean(lengths_all_out), 2)
    std_length_out = round(statistics.stdev(lengths_all_out), 2)

    total_count = sum(list(length_dict.values()))
    length_dict = {i: length_dict[i] / total_count for i in sorted(list(length_dict.keys()))}
    return length_dict, mean_length_out, std_length_out


# Calculate base counts (see SI) from a list of (i,j) values. Returns a dictionary[base index] = fractional base count
def get_base_counts(fragment_ends_in):
    base_count_dict_out = {i: 0 for i in range(1, max([j[1] for j in fragment_ends_in]) + 1)}
    total_fragments = len(fragment_ends_in)
    for fragment in fragment_ends_in:
        start_idx = fragment[0]
        end_idx = fragment[1]
        for i in range(start_idx, end_idx + 1):
            base_count_dict_out[i] += 1

    base_count_dict_out = {i: base_count_dict_out[i] / total_fragments for i in base_count_dict_out.keys()}
    return base_count_dict_out


# Determines the i and j value (see paper) for each read, one (j) denotes the 5'-end position on the upper strand
# (provided in the configuration) the other (i) describes the position of the 5'-end on the complementary strand.
def get_break_counts(fragment_ends_in, window=None, normalize=False):
    break_j_dict = {}
    break_i_dict = {}
    for fragment in fragment_ends_in:
        current_j = fragment[0]
        current_i = fragment[1]
        if current_j not in break_j_dict.keys():
            break_j_dict[current_j] = 1
        else:
            break_j_dict[current_j] += 1
        if current_i not in break_i_dict.keys():
            break_i_dict[current_i] = 1
        else:
            break_i_dict[current_i] += 1

    keys_temp = sorted(list(break_j_dict.keys()))

    # Add zeros to ensure same x axes
    for n in range(1, max(keys_temp) + 1):
        if n not in break_j_dict.keys():
            break_j_dict[n] = 0

    # Optionally restrict to window
    if window:
        break_j_dict_sorted = {n: break_j_dict[n] for n in sorted(list(break_j_dict.keys()))
                               if window[0] < n <= window[1]}
    else:
        break_j_dict_sorted = {n: break_j_dict[n] for n in sorted(list(break_j_dict.keys()))}

    keys_temp = sorted(list(break_i_dict.keys()))
    # Add zeros to ensure same x axes
    for n in range(1, max(keys_temp) + 1):
        if n not in break_i_dict.keys():
            break_i_dict[n] = 0

    # Optionally restrict to window
    if window:
        break_i_dict_sorted = {n: break_i_dict[n] for n in sorted(list(break_i_dict.keys()))
                               if window[0] < n <= window[1]}
    else:
        break_i_dict_sorted = {n: break_i_dict[n] for n in sorted(list(break_i_dict.keys()))}

    # Optionally normalize the breaking distribution
    if normalize:
        sum_j = sum(list(break_j_dict_sorted.values()))
        sum_i = sum(list(break_i_dict_sorted.values()))
        break_j_dict_sorted_normed = {n: break_j_dict_sorted[n] / sum_j for n in break_j_dict_sorted.keys()}
        break_i_dict_sorted_normed = {n: break_i_dict_sorted[n] / sum_i for n in break_i_dict_sorted.keys()}
    else:
        break_i_dict_sorted_normed = break_i_dict_sorted
        break_j_dict_sorted_normed = break_j_dict_sorted

    return break_i_dict_sorted_normed, break_j_dict_sorted_normed


# Calculate the average and standard deviation of break counts for a group of samples. Returns a
# dictionary of dictionaries, contain averages as a dictionary of base index : average break count and
# standard deviations in the same fashion
def average_break_counts(break_counts_group_in):
    x_values_all = [list(i.keys()) for i in break_counts_group_in]
    combined_x_values = []
    for x_values in x_values_all:
        combined_x_values += x_values
    combined_x_values = list(set(combined_x_values))
    averaged_break_counts = {"Average": {}, "Standard Deviation": {}}
    for x_value in combined_x_values:
        values = []
        for break_counts in break_counts_group_in:
            if x_value in break_counts.keys():
                values.append(break_counts[x_value])
        averaged_break_counts["Average"][x_value] = statistics.mean(values)
        if len(values) > 1:
            averaged_break_counts["Standard Deviation"][x_value] = statistics.stdev(values)
        else:
            averaged_break_counts["Standard Deviation"][x_value] = 0

    return averaged_break_counts


def match_reads_single(reads_in, sequences_in, window_in, centroids_in, sample_name_in, output_dir_in, norm=True):
    top_sequence = sequences_in[0]
    bot_sequence = sequences_in[1]
    top_indices = {}
    bot_indices = {}
    count_invalid = 0
    count_invalid_ex = 0

    for data_set in reads_in:
        for read_idx in data_set.keys():
            read_seq = data_set[read_idx]
            res_1 = match_to_reference_single_read_single_sequence(read_seq, top_sequence, no_mismatch_tolerance_distance, maxium_further_missmatches)
            res_2 = match_to_reference_single_read_single_sequence(read_seq, bot_sequence, no_mismatch_tolerance_distance, maxium_further_missmatches)
            
            if res_1[0] == -1 and res_2[0] == -1:
                count_invalid += 1
                if top_sequence != get_complementary_strand(bot_sequence):
                    # Check if taking the complementary sequence reduces the error. This is only relevant/different if top and bot sequence are not complementary
                    res_3 = match_to_reference_single_read_single_sequence(read_seq, get_complementary_strand(top_sequence), no_mismatch_tolerance_distance, maxium_further_missmatches)
                    res_4 = match_to_reference_single_read_single_sequence(read_seq, get_complementary_strand(bot_sequence), no_mismatch_tolerance_distance, maxium_further_missmatches)
                    if res_3[0] == -1 and res_4[0] == -1:
                        count_invalid_ex += 1
            
            if res_1[0] != -1 and res_2[0] != -1:
                print("WARNING: Found a single read on both strands!")
            
            if res_1[0] != -1:
                if res_1[0] not in top_indices.keys():
                    top_indices[res_1[0]] = 1
                else:
                    top_indices[res_1[0]] += 1

            if res_2[0] != -1:
                if res_2[0] not in bot_indices.keys():
                    bot_indices[res_2[0]] = 1
                else:
                    bot_indices[res_2[0]] += 1

    x_vals = sorted([key-centroids_in[0] for key in top_indices.keys()
                     if centroids_in[0]-window_in <= key <= centroids_in[0]+window_in])
    y_vals = [top_indices[i+centroids_in[0]] for i in x_vals]
    
    if norm:
        y_vals = [i/sum(y_vals) for i in y_vals]


    res_dict = {i: y_vals[x_vals.index(i)] for i in x_vals}
    res_dict_new = {f"{sample_name_in}#TOP": res_dict}

    x_vals = sorted([key-centroids_in[1] for key in bot_indices.keys()
                     if centroids_in[1]-window_in <= key <= centroids_in[1]+window_in])
    y_vals = [bot_indices[i+centroids_in[1]] for i in x_vals]

    if norm:
        y_vals = [i/sum(y_vals) for i in y_vals]

    res_dict = {i: y_vals[x_vals.index(i)] for i in x_vals}
    res_dict_new[f"{sample_name_in}#BOT"] = res_dict


    if norm:
        tag = "Norm_"
    else:
        tag = "Full_"

    DataUtilities.print_breaking(res_dict_new, True, f"{sample_name_in}_Breaking_Dist_{tag}Shifted_Single.csv", output_dir_in)
    print(f"{sample_name_in}: Invalid, no complement (single reads)",
          100*count_invalid/sum([len(data) for data in reads_in]))
    if top_sequence != get_complementary_strand(bot_sequence):
        print(f"{sample_name_in}: Invalid, with complement (single reads)",
            100*count_invalid_ex/sum([len(data) for data in reads_in]))

    return top_indices, bot_indices


# Global properties
config_file_path = "../RunConfigurations/Config_ss_Batch5.txt"
create_graphics = True
single_read_analysis = True 
no_mismatch_tolerance_distance = 20
maxium_further_missmatches = 10
gaussian_fit_window = 50

# Read configuration

configuration_good, configuration_data = ConfigurationHandeler.read_config(config_file_path)

if not configuration_good:
    print("Configuration faulty. Aborting.")
    exit(2)

(sample_file_names, sample_names, sample_sequences, sample_sequences_second, sample_nicks,
 data_dir, graphics_dir, data_export_dir) = configuration_data

log_file = open(f"{data_export_dir}log.txt", "w")
formatted_time_start = datetime.now().strftime("%Y-%m-%d, %H:%M:%S")
start_time = time.time()
out_str = (f"Starting analysis on {formatted_time_start}, using no_mismatch_tolerance_distance "
           f"{no_mismatch_tolerance_distance}, maxium_further_missmatches = {maxium_further_missmatches} "
           f"and gaussian_fit_window = {gaussian_fit_window}.")
print(out_str)
log_file.write(out_str + "\n")

sample_groups = {sample_name: sample_name.split("_")[0] for sample_name in sample_names}
paired_read_sample_paths = [[f"{data_dir}{a[0]}", f"{data_dir}{a[1]}"] for a in sample_file_names]

all_sigma = {}
all_mu = {}
all_mu_shifted = {}


# Perform various data analysis steps on each sample
for (paired_read_sample_path, sample_name,
     sample_sequence, secondary_sequence, sample_nick) in zip(paired_read_sample_paths,
                                                              sample_names, sample_sequences, sample_sequences_second, sample_nicks):
    # Read raw data
    data_read_1 = read_input_file(paired_read_sample_path[0])
    data_read_2 = read_input_file(paired_read_sample_path[1])
    out_str = f"{sample_name}: Number of raw single reads are {len(data_read_1)} (read 1) "\
              f"and {len(data_read_2)} (read 2)."
    print(out_str)
    log_file.write(out_str + "\n")

    # Quality analysis

    quality_values = [data_read_1[i][1] for i in data_read_1.keys()]
    quality_values += [data_read_2[i][1] for i in data_read_2.keys()]
    quality_values_translated = [translate_quality_values(i) for i in quality_values]
    quality_values_averages = calculate_quality_average(quality_values_translated)
    quality_values_averages_dict = {sample_name: quality_values_averages}
    DataUtilities.print_quality_results(quality_values_averages_dict, f"{sample_name}_Quality.csv", data_export_dir)

    # Single read analysis
    if single_read_analysis:
        both_sequences = [sample_sequence, secondary_sequence]

        match_reads_single([data_read_1, data_read_2], both_sequences,
                           gaussian_fit_window, [sample_nick[1], sample_nick[3]], sample_name, data_export_dir)
        match_reads_single([data_read_1, data_read_2], both_sequences,
                           gaussian_fit_window, [sample_nick[1], sample_nick[3]], sample_name, data_export_dir, False)
    
    # Combine paired reads from sequence index

    data_paired = match_paired_reads(data_read_1, data_read_2)
    out_str = f"{sample_name}: Number of raw paired reads found is {len(data_paired)}."
    print(out_str)
    log_file.write(out_str + "\n")

    # Match paired reads against reference

    strand_top_53 = sample_sequence
    paired_reads_matched, paired_reads_matched_filtered = (
        match_to_reference(data_paired, strand_top_53, no_mismatch_tolerance_distance, maxium_further_missmatches))

    out_str = f"{sample_name}: Number of valid paired reads is {len(paired_reads_matched_filtered)}, which is "\
              f"{round(len(paired_reads_matched_filtered) / len(paired_reads_matched) * 100, 1)} % of total."
    print(out_str)
    log_file.write(out_str + "\n")

    # Obtain binned filtering results and save them to file

    filtering_results = {sample_name: classify_filtering(paired_reads_matched)}
    DataUtilities.print_filtering_results(filtering_results,
                                          f"{sample_name}_Filtering.csv", data_export_dir)

    # Calculate the 5'-end positions (i and j in the paper)
    i_j_positions = fragment_start_end(paired_reads_matched_filtered)

    # Length analysis
    fragment_lengths, lengths_mean, lengths_std = get_length_dict(i_j_positions)
    DataUtilities.print_lengths(fragment_lengths,
                                f"{sample_name}_Lengths.csv", data_export_dir)

    out_str = f"{sample_name}: Lengths mean: {lengths_mean}, Standard deviation: {lengths_std}."
    print(out_str)
    log_file.write(out_str + "\n")

    # Base count analysis
    base_count_results = {sample_name: get_base_counts(i_j_positions)}
    DataUtilities.print_basecounts(base_count_results,
                                   f"{sample_name}_Basecounts.csv", data_export_dir)

    # Breaking analysis (i = lower strand, j = upper strand with upper = strand from config)
    # Obtain for i + j the following: full breaking histogram, breaking histogram restriced + normed,
    # full breaking histogram shifted, breaking histogram restricted + normed + shifed

    nick_strand = sample_nick[0]
    nick_pos = sample_nick[1]

    i_break_hist_full, j_break_hist_full = get_break_counts(i_j_positions)
    i_break_hist_restr_norm, j_break_hist_restr_norm = (
        get_break_counts(i_j_positions, [nick_pos - gaussian_fit_window, nick_pos + gaussian_fit_window], True))

    # Create shifted datasets
    j_break_hist_full_shifted = {k - nick_pos: j_break_hist_full[k] for k in j_break_hist_full.keys()}
    i_break_hist_full_shifted = {k - nick_pos + 1: i_break_hist_full[k] for k in i_break_hist_full.keys()}
    j_break_hist_restr_norm_shifted = {k - nick_pos: j_break_hist_restr_norm[k]
                                       for k in j_break_hist_restr_norm.keys()}
    i_break_hist_restr_norm_shifted = {k - nick_pos + 1: i_break_hist_restr_norm[k]
                                       for k in i_break_hist_restr_norm.keys()}

    # Prepare for printout
    breaking_dists_full = {f"{sample_name}#TOP": j_break_hist_full,
                           f"{sample_name}#BOT": i_break_hist_full}
    breaking_dists_full_shifted = {f"{sample_name}#TOP": j_break_hist_full_shifted,
                                   f"{sample_name}#BOT": i_break_hist_full_shifted}
    breaking_dists_rests_norm = {f"{sample_name}#TOP": j_break_hist_restr_norm,
                                 f"{sample_name}#BOT": i_break_hist_restr_norm}
    breaking_dists_rests_norm_shifted = {f"{sample_name}#TOP": j_break_hist_restr_norm_shifted,
                                         f"{sample_name}#BOT": i_break_hist_restr_norm_shifted}

    # Save data to file
    DataUtilities.print_breaking(breaking_dists_full, False,
                                 f"{sample_name}_Breaking_Dist_Full.csv", data_export_dir)
    DataUtilities.print_breaking(breaking_dists_full_shifted, True,
                                 f"{sample_name}_Breaking_Dist_Full_Shifted.csv", data_export_dir)
    DataUtilities.print_breaking(breaking_dists_rests_norm, False,
                                 f"{sample_name}_Breaking_Dist_Norm.csv", data_export_dir)
    DataUtilities.print_breaking(breaking_dists_rests_norm_shifted, True,
                                 f"{sample_name}_Breaking_Dist_Norm_Shifted.csv", data_export_dir)

    # Perform Gaussian Fit (only on the non-nicked strand)

    current_mu = 0
    current_sigma = 0
    if nick_strand == 1:
        current_mu, current_sigma = Fitting.fit_gauss(j_break_hist_restr_norm)
        all_sigma[sample_name] = current_sigma
        all_mu[sample_name] = current_mu
        all_mu_shifted[sample_name] = current_mu - nick_pos
    elif nick_strand == 0:
        current_mu, current_sigma = Fitting.fit_gauss(i_break_hist_restr_norm)
        all_sigma[sample_name] = current_sigma
        all_mu[sample_name] = current_mu
        all_mu_shifted[sample_name] = current_mu - nick_pos + 1

    out_str = f"{sample_name}: Mu: {round(current_mu, 2)}, Sigma: {round(current_sigma, 2)}"
    print(out_str)
    log_file.write(out_str + "\n")

# Perform analysis schemes on groups of samples (averages etc.)
for sample_group in set(list(sample_groups.values())):
    group_sample_names_norm = []
    group_sample_names_norm_shifted = []
    group_sample_names_norm_shifted_single = []
    for sample in sample_groups.keys():
        if sample_groups[sample] == sample_group:
            group_sample_names_norm.append(f"{sample}_Breaking_Dist_Norm.csv")
            group_sample_names_norm_shifted.append(f"{sample}_Breaking_Dist_Norm_Shifted.csv")
            group_sample_names_norm_shifted_single.append(f"{sample}_Breaking_Dist_Norm_Shifted_Single.csv")

    # Calculate average break counts with standard deviation
    data_breaking_group_norm_dict = DataUtilities.read_breaking_results(group_sample_names_norm, data_export_dir)
    data_breaking_group_norm_shifted_dict = DataUtilities.read_breaking_results(
        group_sample_names_norm_shifted, data_export_dir)
    data_breaking_group_norm_shifted_single_dict = DataUtilities.read_breaking_results(
        group_sample_names_norm_shifted_single, data_export_dir)

    data_breaking_group_norm_averaged = {}
    data_breaking_group_norm_shifted_averaged = {}
    data_breaking_group_norm_shifted_single_averaged = {}

    for strand in ["TOP", "BOT"]:
        data_breaking_group_norm = [data_breaking_group_norm_dict[i] for i in
                                    data_breaking_group_norm_dict.keys() if i.find(f"{strand}") != -1]
        data_breaking_group_norm_shifted = [data_breaking_group_norm_shifted_dict[i] for i in
                                            data_breaking_group_norm_shifted_dict.keys() if i.find(f"{strand}") != -1]
        data_breaking_group_norm_shifted_single = [data_breaking_group_norm_shifted_single_dict[i] for i in
                                            data_breaking_group_norm_shifted_single_dict.keys() if i.find(f"{strand}") != -1]

        data_breaking_average_norm = average_break_counts(data_breaking_group_norm)
        data_breaking_group_norm_shifted = average_break_counts(data_breaking_group_norm_shifted)
        data_breaking_group_norm_shifted_single = average_break_counts(data_breaking_group_norm_shifted_single)

        data_breaking_group_norm_averaged[
            f"{sample_group}#{strand}#AVERAGE"] = data_breaking_average_norm["Average"]
        data_breaking_group_norm_averaged[
            f"{sample_group}#{strand}#STD"] = data_breaking_average_norm["Standard Deviation"]
        data_breaking_group_norm_shifted_averaged[
            f"{sample_group}#{strand}#AVERAGE"] = data_breaking_group_norm_shifted["Average"]
        data_breaking_group_norm_shifted_averaged[
            f"{sample_group}#{strand}#STD"] = data_breaking_group_norm_shifted["Standard Deviation"]
        data_breaking_group_norm_shifted_single_averaged[
            f"{sample_group}#{strand}#AVERAGE"] = data_breaking_group_norm_shifted_single["Average"]
        data_breaking_group_norm_shifted_single_averaged[
            f"{sample_group}#{strand}#STD"] = data_breaking_group_norm_shifted_single["Standard Deviation"]

    DataUtilities.print_breaking(data_breaking_group_norm_averaged, False,
                                 f"{sample_group}_Breaking_Dist_Norm_Average.csv", data_export_dir)
    DataUtilities.print_breaking(data_breaking_group_norm_shifted_averaged, True,
                                 f"{sample_group}_Breaking_Dist_Norm_Shifted_Average.csv", data_export_dir)
    DataUtilities.print_breaking(data_breaking_group_norm_shifted_single_averaged, True,
                                 f"{sample_group}_Breaking_Dist_Norm_Shifted_Single_Average.csv", data_export_dir)


DataUtilities.print_fit_parameters(all_mu, sample_groups, "Gaussian_Mu_All.csv", data_export_dir)
DataUtilities.print_fit_parameters(all_sigma, sample_groups, "Gaussian_Sigma_All.csv", data_export_dir)
DataUtilities.print_fit_parameters(all_mu_shifted, sample_groups, "Gaussian_Mu_Shifted_All.csv", data_export_dir)

formatted_time_end = datetime.now().strftime("%Y-%m-%d, %H:%M:%S")
end_time = time.time()
total_time = end_time - start_time
minutes = int(total_time // 60)
seconds = int(total_time % 60)
formatted_total_time = f"{minutes:02d}:{seconds:02d}"

out_str = (f"Finished analysis at {formatted_time_end}. Might take a moment to generate graphics.\n"
           f"Total time: {formatted_total_time}")
log_file.write(out_str)
print(out_str)
log_file.close()


# Graphics creation

if create_graphics:
    plotter = PlottingNGS()

    # Single sample graphics (length distribution, breaking distribution for each strand)
    for sample_name, sample_nick in zip(sample_names, sample_nicks):
        centroid = sample_nick[1]
        nick_strand = sample_nick[0]
        if nick_strand == 0:
            gaussian_top = False
        else:
            gaussian_top = True
        
        plotter.plot_single_breaking([f"{sample_name}_Breaking_Dist_Norm.csv"], data_export_dir,
                                     True, gaussian_top, False, False, centroid, "blue",
                                     f"{graphics_dir}{sample_name}_Breaking_Dist_Norm_Top.pdf")
        plotter.plot_single_breaking([f"{sample_name}_Breaking_Dist_Norm.csv"], data_export_dir,
                                     False, not gaussian_top, False, False, centroid, "blue",
                                     f"{graphics_dir}{sample_name}_Breaking_Dist_Norm_Bot.pdf")
        plotter.plot_lengths(f"{sample_name}_Lengths.csv", data_export_dir,
                             f"{graphics_dir}{sample_name}_Lengths.pdf")


    # Grouped sample graphics (average breaking for each strand, compariosn of breaking for each strand,
    # basecouts and filtering results
    for sample_group in set(list(sample_groups.values())):
        sample_names_group = [i for i in sample_groups.keys() if sample_groups[i] == sample_group]
        if len(sample_names_group) != 3:
            continue
        # Plot quality values

        quality_names = [f"{sample_name}_Quality.csv" for sample_name in sample_names_group]
        plotter.plot_quality_values(quality_names, data_export_dir,
                                f"{graphics_dir}{sample_group}_Quality.pdf")

        # Plot basecounts

        base_count_names = [f"{sample_name}_Basecounts.csv" for sample_name in sample_names_group]
        plotter.plot_basecounts(base_count_names, data_export_dir,
                                f"{graphics_dir}{sample_group}_Basecounts.pdf")

        # Plot filtering

        filtering_names = [f"{sample_name}_Filtering.csv" for sample_name in sample_names_group]
        plotter.plot_filtering_results(filtering_names, data_export_dir,
                                f"{graphics_dir}{sample_group}_Filtering.pdf")

        sample_nick_final = (0, 0)
        for sample_name, sample_nick in zip(sample_names, sample_nicks):
            if sample_name in sample_names_group:
                sample_nick_final = sample_nick
                break

        centroid = sample_nick_final[1]
        nick_strand = sample_nick_final[0]
        if nick_strand == 0:
            gaussian_top = False
        else:
            gaussian_top = True

        # Plot average breaking distribution

        plotter.plot_single_breaking([f"{sample_group}_Breaking_Dist_Norm_Shifted_Average.csv"], data_export_dir,
                                True, gaussian_top, True, True, centroid, "blue",
                                     f"{graphics_dir}{sample_group}_Break_Dist_Norm_Shifted_Average_Top.pdf")
        print("Printed Top Average", sample_group)
        plotter.plot_single_breaking([f"{sample_group}_Breaking_Dist_Norm_Shifted_Average.csv"], data_export_dir,
                                False, not gaussian_top, True, True, centroid, "blue",
                                     f"{graphics_dir}{sample_group}_Break_Dist_Norm_Shifted_Average_Bot.pdf")
        print("Printed Bot Average", sample_group)

        breaking_dist_names = [f"{sample_name}_Breaking_Dist_Norm.csv" for sample_name in sample_names_group]

        # Plot comparision of breaking distributions
        plotter.plot_three_breaking(breaking_dist_names, data_export_dir, True, gaussian_top, centroid,
                                    f"{graphics_dir}{sample_group}_Break_Dist_Norm_Comparison_Top.pdf")
        plotter.plot_three_breaking(breaking_dist_names, data_export_dir, False, not gaussian_top, centroid,
                                    f"{graphics_dir}{sample_group}_Break_Dist_Norm_Comparison_Bot.pdf")

    # Whole dataset graphics (mu/sigma scatter plot)

    plotter.plot_fit_paramameter_scatter_single(["Gaussian_Mu_Shifted_All.csv"], data_export_dir, "Mu",
                                                f"{graphics_dir}Mu_Shifted_Scatter_Plot.pdf")

    plotter.plot_fit_paramameter_scatter_single(["Gaussian_Sigma_All.csv"], data_export_dir, "Sigma",
                                                f"{graphics_dir}Sigma_Scatter_Plot.pdf")
