# Script to process data from distance extraction
# Author: Boris N. Schüpp

import sys


def obtain_data(file_path_in):
    file_in = open(file_path_in)

    data_x_out = []
    data_y_out = []

    for line in file_in:
        if line.startswith("#"):
            continue
        if line.startswith("@"):
            continue

        line = line.strip("\n")
        line = [float(a) for a in line.split(" ") if a != ""]
        data_x_out.append(line[0])
        data_y_out.append(line[1:])
    file_in.close()
    return data_x_out, data_y_out


def base_resolve_averages(base_list, time_averages):
    averages_base_resolved = {}
    for k, j in zip(time_averages, base_list):
        if j not in averages_base_resolved.keys():
            averages_base_resolved[j] = [k]
        else:
            averages_base_resolved[j].append(k)

    averages_base_list = [sum(averages_base_resolved[j]) / len(averages_base_resolved[j]) for j in
                          sorted(averages_base_resolved.keys())]

    return averages_base_list


def get_base_list(index_filename):
    file_idx = open(index_filename)
    bases_out = []
    for line in file_idx:
        if line.find("[") == -1:
            continue
        else:
            bases_out.append(int(line.split(".")[1]))
    file_idx.close()
    return bases_out


def get_time_average(cutoff, input_list):
    to_average = input_list[cutoff:]
    average = sum(to_average)/len(to_average)
    return average


def process_end_to_end_distances(out_put_dir_in, out_put_name_in, box_size_in):
    out_str_name = out_put_dir_in + out_put_name_in + "_EndToEndDistance.csv"
    out_str_name_frac = out_put_dir_in + out_put_name_in + "_EndToEndDistance_frac.csv"

    out_file = open(out_str_name, "w")
    out_file_frac = open(out_str_name_frac, "w")

    box_size = box_size_in
    file = "TempEndToEndDistances.xvg"
    data_x, data_y = obtain_data(file)
    data_y_corrected = []

    for i in data_y:
        i = i[0]
        while i >= (box_size + 10):
            i = i - box_size

        data_y_corrected.append(i)

    out_str = ",".join([str(i) for i in data_y_corrected])

    sequence, nick_type, run_num = out_put_name_in.split("_")

    out_str = sequence + "," + str(nick_type) + "," + str(run_num) + "," + out_str
    out_file.write(out_str)
    out_file.close()

    data_y_frac = [i / data_y_corrected[0] for i in data_y_corrected]
    out_str_frac = ",".join([str(i) for i in data_y_frac])
    out_str_frac = sequence + "," + str(nick_type) + "," + str(run_num) + "," + out_str_frac
    out_file_frac.write(out_str_frac)
    out_file_frac.close()


def process_base_distances(out_put_dir_in, out_put_name_in):
    out_str_name = out_put_dir_in + out_put_name_in + "_BaseDistance.csv"
    out_file = open(out_str_name, "w")

    file = f"TempBaseDistances.xvg"
    file_idx_name = f"TempBaseDist.ndx"

    data_x_all, data_y_all = obtain_data(file)
    bases = get_base_list(file_idx_name)
    averages = []
    for i in range(0, len(bases)):
        base_i = [j[i] for j in data_y_all]
        averages.append(get_time_average(1000, base_i))
    averages_bases = base_resolve_averages(bases, averages)

    sequence, nick_type, run_num = out_put_name_in.split("_")
    out_str = ",".join([str(i) for i in averages_bases])
    out_str = sequence + "," + str(nick_type) + "," + str(run_num) + "," + out_str
    out_file.write(out_str)
    out_file.close()


def process_backbone_bond_distances(out_put_dir_in, out_put_name_in):
    out_str_name = out_put_dir_in + out_put_name_in + "_BackboneBondDistances.csv"
    out_file = open(out_str_name, "w")
    file_name = f"TempBackboneDistances.txt"
    input_file = open(file_name)
    current_block = []
    all_blocks = []
    count = 0
    for line in input_file:
        line = line.strip("\n")
        current_block.append(line)
        count += 1
        if count == 4:
            count = 0
            all_blocks.append(current_block)
            current_block = []
    input_file.close()
    sequence, nick_type, run_num = out_put_name_in.split("_")

    for block in all_blocks:
        name = block[0].strip(":")
        bond_type = name.split(".")[0]
        if bond_type != "P-O5'" and bond_type != "O3'-P":
            base_index = int(name.split(".")[1])
        else:
            base_index = int(name.split(".")[1]) + 1
        distance_info = [a for a in block[2].split(" ") if a != ""][2]
        out_file.write(f"{sequence},{nick_type},{run_num},{base_index},{bond_type},{distance_info}\n")
    out_file.close()


if __name__ == "__main__":
    out_put_dir = sys.argv[1]
    out_put_name = sys.argv[2]
    box_size = int(sys.argv[3])
    process_end_to_end_distances(out_put_dir, out_put_name, box_size)
    process_base_distances(out_put_dir, out_put_name)
    process_backbone_bond_distances(out_put_dir, out_put_name)
