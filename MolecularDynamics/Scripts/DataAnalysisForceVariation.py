# Script to perform analysis and visualisation on processed trajectory distance data for force variation simulations
# Author: Boris N. Schüpp

from PlottingMD import *
import statistics
import scipy.stats as stats


def read_base_distances(base_distance_input_file_name_in):
    base_distance_input_file = open(base_distance_input_file_name_in)
    output_dict = {}
    for line in base_distance_input_file:
        line_data = line.strip("\n").split(",")
        force_amount = line_data[1]
        nick_type = int(line_data[0])
        rerun_num = int(line_data[2])
        distances = [float(i) for i in line_data[3:]]
        output_dict[(force_amount, nick_type, rerun_num)] = distances
    return output_dict


def read_end_to_end_distances(end_to_end_distance_input_file_name_in):
    end_to_end_distance_input_file = open(end_to_end_distance_input_file_name_in)
    output_dict = {}
    for line in end_to_end_distance_input_file:
        line_data = line.strip("\n").split(",")
        force_amount = line_data[1]
        nick_type = int(line_data[0])
        rerun_num = int(line_data[2])
        distances = [float(i) for i in line_data[3:]]
        if len(distances) < 5000:
            print("End-to-end distances: Removed", force_amount, nick_type, rerun_num, "due to no complete trajecory.")
            continue
        output_dict[(force_amount, nick_type, rerun_num)] = distances
    return output_dict


def read_backbone_bond_distances(backbone_bond_distance_input_file_name_in):
    backbone_bond_distance_input_file = open(backbone_bond_distance_input_file_name_in)
    output_dict = {}
    for line in backbone_bond_distance_input_file:
        line_data = line.strip("\n").split(",")
        force_amount = line_data[1]
        nick_type = int(line_data[0])
        rerun_num = int(line_data[2])
        base_index = int(line_data[3])
        bond_type = line_data[4]
        distance = float(line_data[5])
        if (force_amount, nick_type, rerun_num, bond_type) not in output_dict.keys():
            output_dict[(force_amount, nick_type, rerun_num, bond_type)] = ([base_index], [distance])
        else:
            output_dict[(force_amount, nick_type, rerun_num, bond_type)][0].append(base_index)
            output_dict[(force_amount, nick_type, rerun_num, bond_type)][1].append(distance)

    return output_dict


def average_arrays(arrays_in):
    average_array = [0 for _ in range(0, max([len(i) for i in arrays_in]))]
    counts_array = [0 for _ in range(0, max([len(i) for i in arrays_in]))]
    for i in range(0, max([len(i) for i in arrays_in])):
        for array in arrays_in:
            if len(array) > i:
                average_array[i] += array[i]
                counts_array[i] += 1
    output_array = [i / j for i, j in zip(average_array, counts_array)]
    return output_array


def calculate_averages(end_to_end_distances_in):
    sequences_to_average = set([key[0] for key in end_to_end_distances_in])

    # Average nicked/non nicked sequence wise
    for sequence in sequences_to_average:
        values_nicked = [end_to_end_distances_in[key] for key in end_to_end_distances_in.keys()
                         if key[0] == sequence and key[1] in [1, 2] and key[2] in [1, 2, 3]]
        values_non_nicked = [end_to_end_distances_in[key] for key in end_to_end_distances_in.keys()
                             if key[0] == sequence and key[1] == 0 and key[2] in [1, 2, 3]]

        average_nicked = average_arrays(values_nicked)
        average_non_nicked = average_arrays(values_non_nicked)
        end_to_end_distances_in[(sequence, 3, -1)] = average_nicked
        end_to_end_distances_in[(sequence, 0, -1)] = average_non_nicked

    return end_to_end_distances_in


def export_data(data_in, x_data_in, data_output_dir_in, data_output_file_in, header_line_in):
    out_file = open(f"{data_output_dir_in}{data_output_file_in}", "w")
    nick_translate = {0: "no nick", 1: "nicked", 2: "nicked", 3: "nicked"}
    name_translate = {"05": "0.5 nN", "10": "1.0 nN", "15": "1.5 nN", "20": "20 nN", "25": "2.5 nN", "30": "3.0 nN"}
    header = header_line_in
    for key in data_in.keys():
        name = ""
        if key[2] == -1:
            name += "Average "
            name += nick_translate[key[1]]
            name += " " + name_translate[key[0]]
        else:
            name += name_translate[key[0]] + " "
            name += nick_translate[key[1]]
            name += f" run {key[2]}"
        header += f",{name}"
    header += "\n"
    out_file.write(header)

    for x_idx, x_val in enumerate(x_data_in):
        line_out = str(x_val)
        for key in data_in.keys():
            line_out += "," + str(data_in[key][x_idx])
        line_out += "\n"
        out_file.write(line_out)
    out_file.close()


def export_data_bond_type(data_in, x_data_in, data_output_dir_in, data_output_file_in, header_line_in, sequences_in):
    out_file = open(f"{data_output_dir_in}{data_output_file_in}", "w")
    sequences = sequences_in
    nick_translate = {3: ""}
    name_translate = {"05": "0.5 nN", "10": "1.0 nN", "15": "1.5 nN", "20": "20 nN", "25": "2.5 nN", "30": "3.0 nN"}
    bond_names = ["P-O5'", "O5'-C5'", "C5'-C4'", "C4'-C3'", "C3'-O3'", "O3'-P"]
    header = header_line_in
    for sequence in sequences:
        for bond in bond_names:
            name_new = name_translate[sequence]
            name_new += " " + bond
            header += f",{name_new}"
    header += "\n"
    out_file.write(header)

    for x_idx, x_val in enumerate(x_data_in):
        line_out = str(x_val)
        for sequence in sequences:
            for bond in bond_names:
                if x_idx + 1 not in data_in[(sequence, 1, -1, bond)][0]:
                    line_out += ","
                else:
                    line_out += "," + str(data_in[(sequence, 1, -1, bond)][1]
                                          [data_in[(sequence, 1, -1, bond)][0].index(x_idx + 1)])
        line_out += "\n"
        out_file.write(line_out)
    out_file.close()


def find_longest_consecutive_fraying(binary_list):
    current_strech = []
    all_streches = []
    for base_idx, base_binary in enumerate(binary_list):
        if base_binary == 1:
            current_strech.append(base_idx + 1)
        else:
            all_streches.append(current_strech)
            current_strech = []
    if len(current_strech) > 0:
        all_streches.append(current_strech)
    for strech in all_streches:
        if 50 not in strech and 51 not in strech:
            continue
        return strech[0], strech[-1], strech[-1] - strech[0]


def calculate_range(input_dict_in, sequences):
    lefts = []
    rights = []
    widths = []

    for sequence in sequences:
        for key in input_dict_in:

            if key[0] != sequence or key[1] not in [1, 2] or key[2] not in [1, 2, 3]:
                continue
            binary = []
            for base in input_dict_in[key]:
                # Filter for high base distances
                if base > 1.0:
                    binary.append(1)
                else:
                    binary.append(0)
            left, right, width = find_longest_consecutive_fraying(binary)
            lefts.append(left)
            rights.append(right)
            widths.append(width)

    print(sequences, "Fraying start: Average", statistics.mean(lefts), "SD:", str(statistics.stdev(lefts)))
    print(sequences, "Fraying end: Average", statistics.mean(rights), "SD:", str(statistics.stdev(rights)))
    print(sequences, "Fraying width: Average", statistics.mean(widths), "SD:", str(statistics.stdev(widths)))

    return lefts, rights, widths


def export_fraying_data(in_05, in_10, in_15, in_20, in_25, in_30, data_output_dir_in, data_output_file_in):
    out_file = open(f"{data_output_dir_in}{data_output_file_in}", "w")

    header = "Simulation number,0.5 nN,1.0 nN,1.5 nN,2.0 nN,2.5 nN,3.0 nN\n"
    out_file.write(header)
    for sim_num, value_cur in enumerate(in_05):
        out_str = str(sim_num + 1) + ","
        out_str += str(value_cur) + ","
        out_str += str(in_10[sim_num]) + ","
        out_str += str(in_15[sim_num]) + ","
        out_str += str(in_20[sim_num]) + ","
        out_str += str(in_25[sim_num]) + ","
        out_str += str(in_30[sim_num])
        out_str += "\n"
        out_file.write(out_str)
    averages_str = f"Average,{statistics.mean(in_05)},{statistics.mean(in_10)},{statistics.mean(in_15)},{statistics.mean(in_20)}, {statistics.mean(in_25)}, {statistics.mean(in_30)}"
    out_file.write(averages_str)
    out_file.close()


def calculate_force_average(force_dict_in):
    # Average for 704 with and without nick

    bond_names = ["P-O5'", "O5'-C5'", "C5'-C4'", "C4'-C3'", "C3'-O3'", "O3'-P"]
    sequences = ["05", "10", "15", "20", "25", "30"]

    for bond in bond_names:
        for sequence in sequences:
            grouped_keys_nick_type_0 = [key for key in force_dict_in.keys() if key[0] == sequence and
                                        key[1] == 0 and key[2] in [1, 2, 3] and key[3] == bond]
            grouped_keys_nick_type_1 = [key for key in force_dict_in.keys() if key[0] == sequence and
                                        key[1] == 1 and key[2] in [1, 2, 3] and key[3] == bond]

            base_idx_no_nick = force_dict_in[grouped_keys_nick_type_0[0]][0]
            grouped_forces_nick_type_0 = [force_dict_in[key][1] for key in grouped_keys_nick_type_0]
            grouped_forces_nick_type_1 = [force_dict_in[key][1] for key in grouped_keys_nick_type_1]

            average_forces_nick_type_0 = average_arrays(grouped_forces_nick_type_0)
            average_forces_nick_type_1 = average_arrays(grouped_forces_nick_type_1)

            force_dict_in[(sequence, 0, -1, bond)] = (base_idx_no_nick, average_forces_nick_type_0)
            force_dict_in[(sequence, 1, -1, bond)] = (base_idx_no_nick, average_forces_nick_type_1)

    return force_dict_in


def convert_to_forces(distance_dict_in):
    bond_names = ["P-O5'", "O5'-C5'", "C5'-C4'", "C4'-C3'", "C3'-O3'", "O3'-P"]
    equ_dists = [0.161, 0.141, 0.1526, 0.1562, 0.141, 0.161]
    force_constants = [192464, 267776, 259408, 259408, 267776, 192464]
    force_conversion = 10 ** (-2) / 6.022140
    force_constants = [a * force_conversion for a in force_constants]
    force_constants_dict = {a: b for a, b in zip(bond_names, force_constants)}
    equ_dists_dict = {a: b for a, b in zip(bond_names, equ_dists)}
    force_dict_out = {}
    for key in distance_dict_in.keys():
        out_x = distance_dict_in[key][0]
        out_y = [abs(i - equ_dists_dict[key[3]]) * force_constants_dict[key[3]] for i in distance_dict_in[key][1]]
        force_dict_out[key] = (out_x, out_y)
    return force_dict_out


# Dictionary naming (NICK-TYPE, FORCE, RUN-NUM),
# FORCE: 05 = 0.5nN, etc.
# NICK-TYPE: 0 = no nick, 1 = nick type C
# RUN-NUM 1,2,3 = run 1,2,3, -1 = Average


plotter = PlottingMD()

# Read in Data
forces = ["05", "10", "15", "20", "25", "30"]
data_export_directory = "../ProcessedData/ForceVariation/"
base_distance_input_file_name = f"{data_export_directory}BaseDistance_all.csv"
end_to_end_distance_input_file_name = f"{data_export_directory}EndToEndDistance_fractional_all.csv"
backbone_bond_distance_input_file_name = f"{data_export_directory}BackboneBondDistance_all.csv"
graphics_directory = "../Graphics/ForceVariation/"
base_distance_data_dict = read_base_distances(base_distance_input_file_name)
end_to_end_distance_data_dict = read_end_to_end_distances(end_to_end_distance_input_file_name)
backbone_bond_distance_data_dict = read_backbone_bond_distances(backbone_bond_distance_input_file_name)

# End-to-End distance analysis

time_values = [0.02 * i for i in range(0, 5001)]
time_based_means = {}
end_to_end_distance_data_dict = calculate_averages(end_to_end_distance_data_dict)

# Plot seperately for each sequence
plotter.plot_end_to_end_distance_all_force_variation(end_to_end_distance_data_dict, time_values,
                                                     f"{graphics_directory}End_to_end_distance_all.pdf")

# Export Figure data

exported_end_to_end_distances = {key: end_to_end_distance_data_dict[key] for key
                                 in end_to_end_distance_data_dict.keys() if key[0] in forces and
                                 key[1] in [0, 1, 3] and key[2] in [1, 2, 3, -1]}

export_data(exported_end_to_end_distances, time_values, data_export_directory,
            "End_to_end_distance_compare_Nick_No_Nick.csv", "Simulation time [ns]")

# Base distance analysis

plot_base_pair_distances_all = calculate_averages(base_distance_data_dict)

# Plot seperately for each sequence
plotter.plot_base_pair_distances_all_force_variation(base_distance_data_dict,
                                                     f"{graphics_directory}Base_Pair_Distance_all.pdf")
exported_base_pair_distances = {key: base_distance_data_dict[key] for key
                                in base_distance_data_dict.keys() if key[0] in forces and
                                key[1] in [0, 1, 3] and key[2] in [1, 2, 3, -1]}
export_data(exported_base_pair_distances, [i for i in range(-49, 51)],
            data_export_directory, "Base_Pair_Distance_Comparison_Nick_No_Nick.csv",
            "Relative base index")

# Calculate fraying start, end, width
left_05, right_05, width_05 = calculate_range(base_distance_data_dict, ["05"])
left_10, right_10, width_10 = calculate_range(base_distance_data_dict, ["10"])
left_15, right_15, width_15 = calculate_range(base_distance_data_dict, ["15"])
left_20, right_20, width_20 = calculate_range(base_distance_data_dict, ["20"])
left_25, right_25, width_25 = calculate_range(base_distance_data_dict, ["25"])
left_30, right_30, width_30 = calculate_range(base_distance_data_dict, ["30"])
values = [width_05, width_10, width_15, width_20, width_25, width_30]

# Statitics

names = ["0.5 nN", "1 nN", "1.5 nN", "2.0 nN", "2.5 nN", "3 nN"]
pvals = []
for value, name in zip(values, names):
    pval_name = []
    pos = names.index(name)
    for value2, name2 in zip(values, names):
        pos2 = names.index(name2)
        pval = stats.f_oneway(value, value2).pvalue
        if pval > 0.05:
            res = "n.s."
        elif pval > 0.01:
            res = "*"
        elif pval > 0.001:
            res = "**"
        else:
            res = "***"
        pval_name.append(res)
    pvals.append(pval_name)

print("      ", names)
for idx, name in enumerate(names):
    print(name, pvals[idx])

# Plot fraying data
plotter.plot_fraying_force_variation(left_05, left_10, left_15, left_20, left_25, left_30, "Left",
                                     f"{graphics_directory}Fraying_Left_Comparison.pdf")
plotter.plot_fraying_force_variation(right_05, right_10, right_15, right_20, right_25, right_30, "Right",
                                     f"{graphics_directory}Fraying_Right_Comparison.pdf")
plotter.plot_fraying_force_variation(width_05, width_10, width_15, width_20, width_25, width_30, "Width",
                                     f"{graphics_directory}Fraying_Width_Comparison.pdf")

export_fraying_data(left_05, left_10, left_15, left_20, left_25, left_30, data_export_directory, "Fraying_Start_Indices.csv")
export_fraying_data(right_05, right_10, right_15, right_20, right_25, right_30, data_export_directory, "Fraying_End_Indices.csv")
export_fraying_data(width_05, width_10, width_15, width_20, width_25, width_30, data_export_directory, "Fraying_Width.csv")

# Bond forces
force_dict = convert_to_forces(backbone_bond_distance_data_dict)
force_dict = calculate_force_average(force_dict)

# Compare bond types
plotter.plot_force_bond_type_all_force_variation(force_dict, f"{graphics_directory}Forces_Comparison_Bond_Type.pdf")
export_data_bond_type(force_dict, [i for i in range(-49, 51)], data_export_directory,
            "Forces_Comparison_Bond_Type.csv", "Relative base index",
                      forces)