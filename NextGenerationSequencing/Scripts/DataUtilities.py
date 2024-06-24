# Utility file for handling and combining datasets
# Author: Boris N. Schüpp

import statistics
import numpy as np


# Read fit parameter from .csv file and return a dictionary of dictionaries dict[sample_group][sample_index] = data
def read_fit_parameters(file_names_in, directory_in):
    return read_multiple_dictionaries(file_names_in, directory_in, "Sample number", "string")


# Read filtering results from .csv file and return a dictionary of dictionaries dict[sample][filtering cat] = data
def read_filtering_results(file_names_in, directory_in):
    return read_multiple_dictionaries(file_names_in, directory_in, "Filtering Category", "string")


# Read basecount results from .csv file and return a dictionary of dictionaries dict[sample][base index] = data
def read_basecount_results(file_names_in, directory_in):
    return read_multiple_dictionaries(file_names_in, directory_in, "Fractional Basecounts", "int")


# Read breaking results from .csv file and return a dictionary of dictionaries dict[sample][base index] = data
def read_breaking_results(file_names_in, directory_in):
    return read_multiple_dictionaries(file_names_in, directory_in, "Base index", "float")


# Read breaking results from .csv file and return a dictionary of dictionaries dict[sample][distance from 5'-end] = data
def read_quality_values(file_names_in, directory_in):
    return read_multiple_dictionaries(file_names_in, directory_in, "Distance from 5'-end", "int")


# General purpose funciton to obtain a dictionary of dictionaries from one or multiple csv files, provide the start of
# the header line as first_line_str and the data type of the x-values (string, int, float). Needs same type of .csv's
# (shared x-axis and same header). Returns dictionary of dictionaries result[sample][xvalue] = data
def read_multiple_dictionaries(file_names_in, directory_in, first_line_str, x_type):
    results_dict = {}
    for file_name in file_names_in:
        sample_names = []
        current_file = open(f"{directory_in}{file_name}", "r")
        for line in current_file:
            line = line.strip("\n")
            # Identify header row
            if line.startswith(first_line_str):
                sample_names = line.split(",")[1:]
                for sample_name in sample_names:
                    results_dict[sample_name] = {}
            else:
                data = line.split(",")[1:]
                x_value = None
                if x_type == "string":
                    x_value = line.split(",")[0]
                elif x_type == "int":
                    x_value = int(line.split(",")[0])
                elif x_type == "float":
                    x_value = float(line.split(",")[0])

                for sample, data in zip(sample_names, data):
                    results_dict[sample][x_value] = float(data)
    return results_dict


# Read fractional fragment lengths results from CSV File and return a dictionary length -> fractional of reads
def read_lengths(file_name_in, directory_in):
    lengths_dict = {}
    file_in = open(f"{directory_in}{file_name_in}", "r")
    for line in file_in:
        # Skip header
        if line.startswith("Fragment"):
            continue

        line = line.strip("\n")
        lengths_dict[int(line.split(",")[0])] = float(line.split(",")[1])
    return lengths_dict


# Read a p matrix from a provided input file and directory, returns the p matrix and the row/column names
def read_p_matrix(file_name_in, directory_in):
    with open(f"{directory_in}{file_name_in}") as input_file:
        names = []
        out_matrix = []
        row_count = 0
        for line in input_file:
            line = line.strip("\n")
            if line.startswith("Compared Sample"):
                names = line.split(",")[1:]
                out_matrix = np.zeros((len(names), len(names)))
            else:
                data = [float(i) for i in line.split(",")[1:]]
                for idx, datum in enumerate(data):
                    out_matrix[row_count][idx] = datum
                row_count += 1
        return out_matrix, names


# Print filtering results dict to a desired output file
def print_filtering_results(filtering_results_dict_in, output_file_name_in, output_dir_in):
    print_multiple_dictionaries(filtering_results_dict_in, "Filtering Category",
                                output_file_name_in, output_dir_in)


# Print quality results dict to a desired output file
def print_quality_results(quality_results_dict_in, output_file_name_in, output_dir_in):
    print_multiple_dictionaries(quality_results_dict_in, "Distance from 5'-end",
                                output_file_name_in, output_dir_in)


# Print base count dictionaries to a desired output file
def print_basecounts(basecounts_in, output_file_name_in, output_dir_in):
    print_multiple_dictionaries(basecounts_in, "Fractional Basecounts",
                                output_file_name_in, output_dir_in)


# Print breaking results dict to a desired output file
def print_breaking(breaking_results_dict_in, shifted, output_file_name_in, output_dir_in):
    if shifted:
        print_multiple_dictionaries(breaking_results_dict_in, "Base index shifted",
                                    output_file_name_in, output_dir_in)
    else:
        print_multiple_dictionaries(breaking_results_dict_in, "Base index",
                                    output_file_name_in, output_dir_in)


# Print fragment length dictionaries results dict to a desired output file
def print_lengths(lengths_in, output_file_name_in, output_dir_in):
    with open(f"{output_dir_in}{output_file_name_in}", "w") as output_file:
        fieldnames = "Fragment length, Fraction of reads\n"
        output_file.write(fieldnames)
        for length in lengths_in.keys():
            out_str = f"{str(length)},{str(lengths_in[length])}\n"
            output_file.write(out_str)


# Print fit parameter dict to a desired output file
def print_fit_parameters(fit_parameters_in, groups_dict_in, output_file_name_in, output_dir_in):
    grouped_parameters = {}
    for key in fit_parameters_in.keys():
        if groups_dict_in[key] not in grouped_parameters.keys():
            grouped_parameters[groups_dict_in[key]] = {}
        grouped_parameters[groups_dict_in[key]][key.split("_")[1]] = fit_parameters_in[key]

    for key in grouped_parameters.keys():
        grouped_parameters[key]["Average"] = statistics.mean([grouped_parameters[key][i]
                                                              for i in grouped_parameters[key].keys()])

    print_multiple_dictionaries(grouped_parameters, "Sample number", output_file_name_in, output_dir_in)


# General purpose funcition to print a dictionary of dictionaies to a desired output csv file. Provide the name
# of the x_axis (x_name) which needs to be shared across the dictionary of dictionaries.
def print_multiple_dictionaries(data_in, x_name, output_file_name_in, output_dir_in):
    filtering_results_sample_names = list(data_in.keys())

    all_keys = []
    for temp_dict in data_in.values():
        all_keys += list(temp_dict.keys())
    x_values = sorted(list(set(all_keys)))

    out_file_path = f"{output_dir_in}{output_file_name_in}"
    with open(out_file_path, "w") as filtering_out_file:
        fieldnames = f"{x_name}"
        for sample_name in filtering_results_sample_names:
            fieldnames += f",{sample_name}"
        fieldnames += "\n"
        filtering_out_file.write(fieldnames)
        for x_value in x_values:
            data_string = str(x_value)
            for sample in filtering_results_sample_names:
                if x_value in data_in[sample].keys():
                    data_string += f",{data_in[sample][x_value]}"
                else:
                    data_string += f",0"
            data_string += "\n"
            filtering_out_file.write(data_string)


# Print a p matrix to an output file and directory. Row/Column names need to be provided
def print_p_matrix(p_matrix_in, names_in, output_file_name_in, output_dir_in):
    with open(f"{output_dir_in}{output_file_name_in}", "w") as out_file:
        header_str = "Compared Sample"
        for name in names_in:
            header_str += "," + name
        out_file.write(header_str + "\n")

        for index, name in enumerate(names_in):
            out_str = name
            current_row = p_matrix_in[index]
            for value in current_row:
                out_str += "," + str(value)
            out_file.write(out_str + "\n")


if __name__ == '__main__':
    sample = "ATMi4"
    fig_num = "57"
    data = read_breaking_results([f"{sample}_1_Breaking_Dist_Norm.csv", f"{sample}_2_Breaking_Dist_Norm.csv", f"{sample}_3_Breaking_Dist_Norm.csv"], "../ProcessedData/704_AT_GC_400_1500_ATMi/")
    print_breaking(data, False, f"FigS{fig_num}.csv", f"C:/Users/Boris/Desktop/Data/Supplement/NEW/FigS{fig_num}/")
