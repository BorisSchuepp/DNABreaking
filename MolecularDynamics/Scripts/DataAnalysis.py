# Script to perform analysis and visualisation on processed trajectory distance data
# Author: Boris N. Schüpp

from PlottingMD import *
import statistics


def read_base_distances(base_distance_input_file_name_in):
    base_distance_input_file = open(base_distance_input_file_name_in)
    output_dict = {}
    for line in base_distance_input_file:
        line_data = line.strip("\n").split(",")
        sequence = line_data[0]
        nick_type = int(line_data[1])
        rerun_num = int(line_data[2])
        distances = [float(i) for i in line_data[3:]]
        if (sequence == "GAT" and nick_type == 2) or (sequence == "GCT" and nick_type == 2 and rerun_num == 0):
            print(f"Base distances: Removed {sequence} {nick_type} {rerun_num} due to no complete trajectory or "
                  f"error in pull group.", )
            continue
        output_dict[(sequence, nick_type, rerun_num)] = distances
    return output_dict


def read_end_to_end_distances(end_to_end_distance_input_file_name_in):
    end_to_end_distance_input_file = open(end_to_end_distance_input_file_name_in)
    output_dict = {}
    for line in end_to_end_distance_input_file:
        line_data = line.strip("\n").split(",")
        sequence = line_data[0]
        nick_type = int(line_data[1])
        rerun_num = int(line_data[2])
        distances = [float(i) for i in line_data[3:]]
        if len(distances) < 5000:
            print("End-to-end distances: Removed", sequence, nick_type, rerun_num, "due to no complete trajecory.")
            continue
        output_dict[(sequence, nick_type, rerun_num)] = distances
    return output_dict


def read_backbone_bond_distances(backbone_bond_distance_input_file_name_in):
    backbone_bond_distance_input_file = open(backbone_bond_distance_input_file_name_in)
    output_dict = {}
    for line in backbone_bond_distance_input_file:
        line_data = line.strip("\n").split(",")
        sequence = line_data[0]
        nick_type = int(line_data[1])
        rerun_num = int(line_data[2])
        base_index = int(line_data[3])
        bond_type = line_data[4]
        distance = float(line_data[5])
        if (sequence, nick_type, rerun_num, bond_type) not in output_dict.keys():
            output_dict[(sequence, nick_type, rerun_num, bond_type)] = ([base_index], [distance])
        else:
            output_dict[(sequence, nick_type, rerun_num, bond_type)][0].append(base_index)
            output_dict[(sequence, nick_type, rerun_num, bond_type)][1].append(distance)

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
                         if key[0] == sequence and key[1] in [1, 2] and key[2] in [0, 1, 2]]
        values_non_nicked = [end_to_end_distances_in[key] for key in end_to_end_distances_in.keys()
                             if key[0] == sequence and key[1] == 0 and key[2] in [0, 1, 2]]

        average_nicked = average_arrays(values_nicked)
        average_non_nicked = average_arrays(values_non_nicked)
        end_to_end_distances_in[(sequence, 3, -1)] = average_nicked
        end_to_end_distances_in[(sequence, 0, -1)] = average_non_nicked

    # Average nicked/non-nicked 704*

    averages_nicked_704 = average_arrays([end_to_end_distances_in[key] for key in end_to_end_distances_in.keys()
                                          if key[0] not in ["GC", "AT"] and key[1] in [1, 2] and key[2] in [0, 1, 2]])
    averages_non_nicked_704 = average_arrays([end_to_end_distances_in[key] for key in end_to_end_distances_in.keys()
                                              if key[0] not in ["GC", "AT"] and key[1] == 0 and key[2] in [0, 1, 2]])
    end_to_end_distances_in[("704", 3, -1)] = averages_nicked_704
    end_to_end_distances_in[("704", 0, -1)] = averages_non_nicked_704

    # Average all non-nicked
    averages_non_nicked_all = average_arrays([end_to_end_distances_in[key] for key in end_to_end_distances_in.keys()
                                              if key[1] == 0 and key[2] in [0, 1, 2]])

    end_to_end_distances_in[("ALL", 0, -1)] = averages_non_nicked_all

    return end_to_end_distances_in


def export_data(data_in, x_data_in, data_output_dir_in, data_output_file_in, header_line_in):
    out_file = open(f"{data_output_dir_in}{data_output_file_in}", "w")
    nick_translate = {0: "no nick", 1: "nicked", 2: "nicked", 3: "nicked"}

    header = header_line_in
    for key in data_in.keys():
        name = ""
        if key[0] == "704" or key[2] == -1:
            name += "Average "
            name += nick_translate[key[1]]
            name += " " + key[0]
        else:
            name += key[0] + " "
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


def export_data_nick_type(data_in, x_data_in, data_output_dir_in, data_output_file_in, header_line_in, sequences_in):
    out_file = open(f"{data_output_dir_in}{data_output_file_in}", "w")
    sequences = sequences_in
    nick_translate = {1: "nick type B", 2: "nick type C"}

    header = header_line_in
    for sequence in sequences:
        for nick_type in [1, 2]:
            for run in [0, 1, 2, -1]:
                name = ""
                if run == -1:
                    name += "Average "
                    name += nick_translate[nick_type]
                    name += " " + sequence
                else:
                    name += sequence + " "
                    name += nick_translate[nick_type]
                    name += f" run {run}"
                header += f",{name}"
    header += "\n"
    out_file.write(header)

    for x_idx, x_val in enumerate(x_data_in):
        line_out = str(x_val)
        for sequence in sequences:
            for nick_type in [1, 2]:
                for run in [0, 1, 2, -1]:
                    line_out += "," + str(data_in[(sequence, nick_type, run, "C3'-O3'")][x_idx])
        line_out += "\n"
        out_file.write(line_out)
    out_file.close()


def export_data_bond_type(data_in, x_data_in, data_output_dir_in, data_output_file_in, header_line_in, sequences_in):
    out_file = open(f"{data_output_dir_in}{data_output_file_in}", "w")
    sequences = sequences_in
    nick_translate = {3: "both nick types"}
    bond_names = ["P-O5'", "O5'-C5'", "C5'-C4'", "C4'-C3'", "C3'-O3'", "O3'-P"]

    header = header_line_in
    for sequence in sequences:
        for bond in bond_names:
            name = "Average both nick types"
            name += " " + sequence
            name += " " + bond
            header += f",{name}"
    header += "\n"
    out_file.write(header)

    for x_idx, x_val in enumerate(x_data_in):
        line_out = str(x_val)
        for sequence in sequences:
            for bond in bond_names:
                if x_idx+1 not in data_in[(sequence, 3, -1, bond)][0]:
                    line_out += ","
                else:
                    line_out += "," + str(data_in[(sequence, 3, -1, bond)][1]
                                          [data_in[(sequence, 3, -1, bond)][0].index(x_idx+1)])
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
            if key[0] != sequence or key[1] not in [1, 2] or key[2] not in [0, 1, 2]:
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


def export_fraying_data(in_gc, in_704, in_at, data_output_dir_in, data_output_file_in):
    out_file = open(f"{data_output_dir_in}{data_output_file_in}", "w")

    header = "Simulation number, GC*, 704*-like, AT*\n"
    out_file.write(header)
    for sim_num, value in enumerate(in_704):
        out_str = str(sim_num + 1) + ","
        if sim_num < len(in_gc):
            out_str += str(in_gc[sim_num])
        out_str += "," + str(value) + ","
        if sim_num < len(in_at):
            out_str += str(in_at[sim_num])
        out_str += "\n"
        out_file.write(out_str)
    averages_str = f"Average,{statistics.mean(in_gc)},{statistics.mean(in_704)},{statistics.mean(in_at)}"
    out_file.write(averages_str)
    out_file.close()


def calculate_force_average(force_dict_in):
    # Average for 704 with and without nick

    bond_names = ["P-O5'", "O5'-C5'", "C5'-C4'", "C4'-C3'", "C3'-O3'", "O3'-P"]
    sequences = ["AAT", "ACT", "AGT", "GAT", "GGT", "GCT", "GTT"]
    sequences_all = ["AAT", "ACT", "AGT", "GAT", "GGT", "GCT", "GTT", "GC", "AT"]

    # Averages for 704-like sequences
    for bond in bond_names:
        grouped_keys_no_nick = [key for key in force_dict_in.keys() if key[0] in sequences and
                                key[1] == 0 and key[2] in [0, 1, 2] and key[3] == bond]
        grouped_forces_no_nick = [force_dict_in[key][1] for key in grouped_keys_no_nick]
        base_idx_no_nick = force_dict_in[grouped_keys_no_nick[0]][0]
        average_force_no_nick = average_arrays(grouped_forces_no_nick)
        force_dict_in[("704", 0, -1, bond)] = (base_idx_no_nick, average_force_no_nick)

        grouped_keys_nick = [key for key in force_dict_in.keys() if key[0] in sequences and
                             key[1] in [1, 2] and key[2] in [0, 1, 2] and key[3] == bond]
        grouped_forces_nick = [force_dict_in[key][1] for key in grouped_keys_nick]
        base_idx_nick = force_dict_in[grouped_keys_nick[0]][0]
        average_force_nick = average_arrays(grouped_forces_nick)
        force_dict_in[("704", 3, -1, bond)] = (base_idx_nick, average_force_nick)

    for bond in bond_names:
        for sequence in sequences:
            grouped_keys_nick_type_b = [key for key in force_dict_in.keys() if key[0] == sequence and
                                        key[1] == 1 and key[2] in [0, 1, 2] and key[3] == bond]
            grouped_keys_nick_type_c = [key for key in force_dict_in.keys() if key[0] == sequence and
                                        key[1] == 2 and key[2] in [0, 1, 2] and key[3] == bond]

            base_idx_no_nick = force_dict_in[grouped_keys_nick_type_b[0]][0]
            grouped_forces_nick_type_b = [force_dict_in[key][1] for key in grouped_keys_nick_type_b]
            grouped_forces_nick_type_c = [force_dict_in[key][1] for key in grouped_keys_nick_type_c]

            average_forces_nick_type_b = average_arrays(grouped_forces_nick_type_b)
            average_forces_nick_type_c = average_arrays(grouped_forces_nick_type_c)

            force_dict_in[(sequence, 1, -1, bond)] = (base_idx_no_nick, average_forces_nick_type_b)
            force_dict_in[(sequence, 2, -1, bond)] = (base_idx_no_nick, average_forces_nick_type_c)

    # Sequence + bond specific average over both nick types
    for bond in bond_names:
        for sequence in sequences_all:
            grouped_keys_nick_type_both = [key for key in force_dict_in.keys() if key[0] == sequence and
                                        key[1] in [1, 2] and key[2] in [0, 1, 2] and key[3] == bond]
            base_idx_no_nick = force_dict_in[grouped_keys_nick_type_both[0]][0]
            grouped_forces_nick_type_both = [force_dict_in[key][1] for key in grouped_keys_nick_type_both]
            average_forces_nick_type_both = average_arrays(grouped_forces_nick_type_both)
            force_dict_in[(sequence, 3, -1, bond)] = (base_idx_no_nick, average_forces_nick_type_both)

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


# Dictionary naming (SEQUENCE, NICK-TYPE, RUN-NUM),
# NICK-TYPE: 0 = no nick, 1 = nick type B, 2 = nick type C, 3 average over both nick types
# RUN-NUM 0,1,2 = run 0,1,2, -1 = Average


plotter = PlottingMD()
sequences_704 = ["AAT", "ACT", "AGT", "GAT", "GGT", "GCT", "GTT"]

# Read in Data
data_export_directory = "../ProcessedData/GCContent/"
base_distance_input_file_name = f"{data_export_directory}BaseDistance_all.csv"
end_to_end_distance_input_file_name = f"{data_export_directory}EndToEndDistance_fractional_all.csv"
backbone_bond_distance_input_file_name = f"{data_export_directory}BackboneBondDistance_all.csv"
graphics_directory = "../Graphics/GCContent/"
base_distance_data_dict = read_base_distances(base_distance_input_file_name)
end_to_end_distance_data_dict = read_end_to_end_distances(end_to_end_distance_input_file_name)
backbone_bond_distance_data_dict = read_backbone_bond_distances(backbone_bond_distance_input_file_name)

# End-to-End distance analysis

time_values = [0.02 * i for i in range(0, 5001)]
end_to_end_distance_data_dict = calculate_averages(end_to_end_distance_data_dict)

# Plot seperately for each sequence
plotter.plot_end_to_end_distance_all(end_to_end_distance_data_dict, time_values,
                                     f"{graphics_directory}End_to_end_distance_all.pdf")

# Plot for set 704* (AAT, ACT, AGT, GAT, GCT, GGT, GTT) nicked vs non-nicked
plotter.plot_end_to_end_distance_nick_no_nick(end_to_end_distance_data_dict, time_values,
                                              f"{graphics_directory}End_to_end_distance_Comparison_Nick_No_Nick.pdf")

# Export Figure data
exported_end_to_end_distances = {key: end_to_end_distance_data_dict[key] for key
                                 in end_to_end_distance_data_dict.keys() if key[0] in sequences_704 and
                                 key[1] in [0, 1, 2] and key[2] in [0, 1, 2]}
exported_end_to_end_distances[("704", 3, -1)] = end_to_end_distance_data_dict[("704", 3, -1)]
exported_end_to_end_distances[("704", 0, -1)] = end_to_end_distance_data_dict[("704", 0, -1)]

export_data(exported_end_to_end_distances, time_values, data_export_directory,
            "End_to_end_distance_compare_Nick_No_Nick.csv", "Simulation time [ns]")

# Base distance analysis

plot_base_pair_distances_all = calculate_averages(base_distance_data_dict)

# Plot seperately for each sequence
plotter.plot_base_pair_distances_all(base_distance_data_dict, f"{graphics_directory}Base_Pair_Distance_all.pdf")

# Plot for set 704* (AAT, ACT, AGT, GAT, GCT, GGT, GTT) nicked vs non-nicked
plotter.plot_base_pair_distance_nick_no_nick(base_distance_data_dict,
                                             f"{graphics_directory}Base_Pair_Distance_Comparison_Nick_No_Nick.pdf")

# Export Figure data
exported_base_pair_distances = {key: base_distance_data_dict[key] for key
                                in base_distance_data_dict.keys() if key[0] in sequences_704 and
                                key[1] in [0, 1, 2] and key[2] in [0, 1, 2]}
exported_base_pair_distances[("704", 3, -1)] = base_distance_data_dict[("704", 3, -1)]
exported_base_pair_distances[("704", 0, -1)] = base_distance_data_dict[("704", 0, -1)]
export_data(exported_base_pair_distances, [i for i in range(-49, 51)],
            data_export_directory, "Base_Pair_Distance_Comparison_Nick_No_Nick.csv",
            "Relative base index")

# Plot averages
plotter.plot_base_pair_distance_averages(base_distance_data_dict,
                                         f"{graphics_directory}Base_Pair_Distance_Comparison_GC_AT_704.pdf")

# Export Figure data
exported_base_pair_distances = {key: base_distance_data_dict[key] for key in base_distance_data_dict.keys()
                                if key[0] in ["GC", "AT", "704"] and key[1] == 3 and key[2] == -1}
export_data(exported_base_pair_distances, [i for i in range(-49, 51)],
            data_export_directory, "Base_Pair_Distance_Comparison_GC_AT_704.csv",
            "Relative base index")

# Calculate fraying start, end, width
left_704, right_704, width_704 = calculate_range(base_distance_data_dict,
                                                 ["AAT", "ACT", "AGT", "GAT", "GGT", "GCT", "GTT"])
left_at, right_at, width_at = calculate_range(base_distance_data_dict,
                                              ["AT"])
left_gc, right_gc, width_gc = calculate_range(base_distance_data_dict,
                                              ["GC"])

# Export to file
export_fraying_data(left_gc, left_704, left_at, data_export_directory, "Fraying_Start_Indices.csv")
export_fraying_data(right_gc, right_704, right_at, data_export_directory, "Fraying_End_Indices.csv")
export_fraying_data(width_gc, width_704, width_at, data_export_directory, "Fraying_Width.csv")

# Plot fraying data
plotter.plot_fraying(left_gc, left_704, left_at, "Left",
                     f"{graphics_directory}Fraying_Left_Comparison.pdf")
plotter.plot_fraying(right_gc, right_704, right_at, "Right",
                     f"{graphics_directory}Fraying_Right_Comparison.pdf")
plotter.plot_fraying(width_gc, width_704, width_at, "Width",
                     f"{graphics_directory}Fraying_Width_Comparison.pdf")

# Bond forces

force_dict = convert_to_forces(backbone_bond_distance_data_dict)
force_dict = calculate_force_average(force_dict)

# Plot for set 704* (AAT, ACT, AGT, GAT, GCT, GGT, GTT) nicked vs non-nicked
plotter.plot_forces_nick_no_nick(force_dict, f"{graphics_directory}Forces_Comparison_Nick_No_Nick.pdf")

# Export Figure data
exported_force_data = {key: force_dict[key][1] for key
                       in force_dict.keys() if key[0] in sequences_704 and
                       key[1] in [0, 1, 2] and key[2] in [0, 1, 2] and key[3] == "C3'-O3'"}
exported_force_data[("704", 3, -1, "C3'-O3'")] = force_dict[("704", 3, -1, "C3'-O3'")][1]
exported_force_data[("704", 0, -1, "C3'-O3'")] = force_dict[("704", 0, -1, "C3'-O3'")][1]

export_data(exported_force_data, [i for i in range(-49, 51)], data_export_directory,
            "Forces_Comparison_Nick_No_Nick.csv", "Relative base index")

# Compare nick types

# Plot all sequences besides GTT
plotter.plot_force_nick_type_all(force_dict, f"{graphics_directory}Forces_Comparison_Nick_Type.pdf")

exported_force_data = {key: force_dict[key][1] for key
                       in force_dict.keys() if key[0] in ["AAT", "ACT", "AGT", "GAT", "GCT", "GGT"] and
                       key[1] in [1, 2] and key[2] in [0, 1, 2, -1] and key[3] == "C3'-O3'"}
export_data_nick_type(exported_force_data, [i for i in range(-49, 51)], data_export_directory,
            "Forces_Comparison_Nick_Type.csv", "Relative base index",
                      ["AAT", "ACT", "AGT", "GAT", "GCT", "GGT"])

# Plot GTT seperate
plotter.plot_force_nick_type_single(force_dict, "GTT",
                                    f"{graphics_directory}Forces_Comparison_Nick_Type_GTT.pdf")

exported_force_data = {key: force_dict[key][1] for key
                       in force_dict.keys() if key[0] == "GTT" and
                       key[1] in [1, 2] and key[2] in [0, 1, 2, -1] and key[3] == "C3'-O3'"}
export_data_nick_type(exported_force_data, [i for i in range(-49, 51)], data_export_directory,
            "Forces_Comparison_Nick_Type_GTT.csv", "Relative base index", ["GTT"])

# Compare bond types

plotter.plot_force_bond_type_all(force_dict, f"{graphics_directory}Forces_Comparison_Bond_Type.pdf")
export_data_bond_type(force_dict, [i for i in range(-49, 51)], data_export_directory,
            "Forces_Comparison_Bond_Type.csv", "Relative base index",
                      ["AAT", "ACT", "AGT", "GAT", "GCT", "GGT", "GTT", "AT", "GC"])

plotter.plot_force_bond_type_single(force_dict, "GTT",
                                    f"{graphics_directory}Forces_Comparison_Bond_Type_GTT.pdf")

export_data_bond_type(force_dict, [i for i in range(-49, 51)], data_export_directory,
            "Forces_Comparison_Bond_Type_GTT.csv", "Relative base index",
                      ["GTT"])
