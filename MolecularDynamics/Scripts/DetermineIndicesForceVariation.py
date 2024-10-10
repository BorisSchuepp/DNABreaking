# Script to determine atom indices that are relevant for end-to-end distance, force and basepair distance analysis
# Author: Boris N. Schüpp

import os
import sys


def determine_indices_end_to_end(gro_file_name_in, structure_directory_in, file_name_out):
    gro_file = open(f"{structure_directory}{gro_file_name_in}", "r")
    strand_index = 0
    current_res_num = 1
    all_res = []
    all_res_cur_strand = []
    current_res = []

    for line in gro_file:
        if line.find("D") == -1:
            continue
        line = line.strip("\n")
        line = [a for a in line.split(" ") if a != ""][0:6]

        res_num = int(line[0][0:line[0].find("D")])
        line = list(line)
        line[0] = res_num
        if res_num == current_res_num:
            current_res.append(line)
        else:
            all_res_cur_strand.append(current_res)
            current_res = [line]
            current_res_num = res_num
            if res_num == 1:
                all_res.append(all_res_cur_strand)
                all_res_cur_strand = []

    all_res_cur_strand.append(current_res)
    all_res.append(all_res_cur_strand)
    gro_file.close()

    three_prime_ends = []
    for strand in all_res:
        for res in strand:
            for atom in res:
                if atom[1] == "H3T":
                    three_prime_ends.append(atom)

    if len(three_prime_ends) > 2:
        # print("WARNING: Found more than two 3'-Ends, will determine the two with the "
        #      "largest distance for end to end distance calculations!")
        max_distance = -1
        max_pair = [-1, -1]
        for i in range(0, len(three_prime_ends)):
            for j in range(i + 1, len(three_prime_ends)):
                coord_a = [float(i) for i in three_prime_ends[i][3:]]
                coord_b = [float(i) for i in three_prime_ends[j][3:]]
                distance = (sum([(k - l) ** 2 for k, l in zip(coord_a, coord_b)])) ** (1 / 2)
                if distance > max_distance:
                    max_distance = distance
                    max_pair = [i, j]

        atom_index_a = three_prime_ends[max_pair[0]][2]
        atom_index_b = three_prime_ends[max_pair[1]][2]

    elif len(three_prime_ends) != 2:
        print("ERROR: Found less than two 3'-Ends!")
        atom_index_a = -1
        atom_index_b = -1
        exit(2)
    else:
        atom_index_a = three_prime_ends[0][2]
        atom_index_b = three_prime_ends[1][2]
    output_file = open(file_name_out, "w")
    output_file.write(f"[ EndToEndDistance ]\n")
    output_file.write(f"{atom_index_a} {atom_index_b}\n")
    output_file.close()


def find_input_file(sequence_in, nick_type_in):
    files_list = []
    for root, dirs, files in os.walk(structure_directory):
        for file in files:
            relative_path = os.path.relpath(os.path.join(root, file), structure_directory)
            files_list.append(relative_path)
    search_string = "/" + sequence_in + "_" + str(nick_type_in)

    gro_file_name_out = ""
    for file_name in files_list:
        if file_name.find(search_string) != -1 and file_name.find(".gro") != -1:
            gro_file_name_out = file_name

    if gro_file_name_out == "":
        print("ERROR: Could not find matching .gro file!")
        exit(2)
    return gro_file_name_out


def determine_indices_base_distance(gro_file_name_in, structure_directory_in, file_name_out):
    gro_file = open(f"{structure_directory}{gro_file_name_in}", "r")
    output_file = open(file_name_out, "w")

    strand_index = 0
    current_res_num = 1
    all_res = []
    all_res_cur_strand = []
    current_res = []

    for line in gro_file:
        if line.find("D") == -1:
            continue
        line = line.strip("\n")
        line = [a for a in line.split(" ") if a != ""][0:6]

        res_num = int(line[0][0:line[0].find("D")])
        line.append(line[0][line[0].find("D") + 1])
        line = list(line)
        line[0] = res_num
        if res_num == current_res_num:
            current_res.append(line)
        else:
            all_res_cur_strand.append(current_res)
            current_res = [line]
            current_res_num = res_num
            if res_num == 1:
                all_res.append(all_res_cur_strand)
                all_res_cur_strand = []

    all_res_cur_strand.append(current_res)
    all_res.append(all_res_cur_strand)
    gro_file.close()

    principal_strand = (all_res[[len(i) for i in all_res].index(max([len(i) for i in all_res]))])
    principal_idx = [len(i) for i in all_res].index(max([len(i) for i in all_res]))
    info_strand = ""
    for i in range(0, 10):
        info_strand = info_strand + str(principal_strand[i][0][-1]) + " "

    g_h = ["O6", "H1", "H21"]
    c_h = ["H41", "N3", "O2"]
    a_h = ["N1", "H61"]
    t_h = ["O4", "H3"]
    base_h_dict = {"G": g_h, "C": c_h, "A": a_h, "T": t_h}

    principal_h = []

    for res in principal_strand:
        base = res[0][-1]
        for atom in res:
            if atom[1] in base_h_dict[base]:
                principal_h.append(atom)

    non_principal_h = []
    for i in range(0, len(all_res)):
        if i == principal_idx:
            continue
        cur_strand = all_res[i]
        for res in cur_strand:
            base = res[0][-1]
            for atom in res:
                if atom[1] in base_h_dict[base]:
                    non_principal_h.append(atom)

    partners = {"O6": ["H41"], "H41": ["O6"], "H1": ["N3"], "N3": ["H1"],
                "H21": ["O2"], "O2": ["H21"], "N1": ["H3"], "H3": ["N1"], "H61": ["O4"], "O4": ["H61"]}

    for atom in principal_h:
        possible_partners = partners[atom[1]]
        partner_list = [i for i in non_principal_h if i[1] in possible_partners]
        coords = [[float(i[3]), float(i[4]), float(i[5])] for i in partner_list]
        coord_a = [float(atom[3]), float(atom[4]), float(atom[5])]

        distances = [((i[0] - coord_a[0]) ** 2 + (i[1] - coord_a[1]) ** 2 + (i[2] - coord_a[2]) ** 2) ** (1 / 2) for i
                     in coords]
        partner = partner_list[distances.index(min(distances))]
        output_file.write(f"[ BaseDistance.{atom[0]}.{atom[1]} ]\n")
        output_file.write(f"{atom[2]} {partner[2]}\n")
    output_file.close()


def determine_indices_backbone_bonds(gro_file_name_in, structure_directory_in, file_name_out):
    gro_file = open(f"{structure_directory}{gro_file_name_in}", "r")
    output_file = open(file_name_out, "w")

    current_res_num = 1
    current_res = []
    all_res = []

    for line in gro_file:
        if line.find("D") == -1:
            continue
        line = line.strip("\n")
        line = [a for a in line.split(" ") if a != ""][0:3]

        res_num = int(line[0][0:line[0].find("D")])
        line = list(line)
        line[0] = res_num

        if res_num == current_res_num:
            current_res.append(line)
        else:
            all_res.append(current_res)
            current_res = [line]
            current_res_num = res_num
            if res_num == 1:
                break

    gro_file.close()

    interesting_atoms = ["O5'", "C5'", "C4'", "C3'", "O3'", "P"]
    atom_dict = {"P": [], "O5'": [], "C5'": [], "C4'": [], "C3'": [], "O3'": []}
    for res in all_res:
        for atom in res:
            if atom[1] in interesting_atoms:
                atom_dict[atom[1]].append((atom[0], int(atom[2])))

    p_atoms = atom_dict["P"]
    o5_atoms = atom_dict["O5'"]
    c5_atoms = atom_dict["C5'"]
    c4_atoms = atom_dict["C4'"]
    c3_atoms = atom_dict["C3'"]
    o3_atoms = atom_dict["O3'"]

    p_o5_bond_pairs = []
    o5_c5_bond_pairs = []
    c5_c4_bond_pairs = []
    c4_c3_bond_pairs = []
    c3_o3_bond_pairs = []
    o3_p_bond_pairs = []

    for atom_list_idx in range(0, 99):
        p_o5_bond_pairs.append((p_atoms[atom_list_idx][1], o5_atoms[atom_list_idx + 1][1]))
        o3_p_bond_pairs.append((o3_atoms[atom_list_idx][1], p_atoms[atom_list_idx][1]))

    for atom_list_idx in range(0, 100):
        o5_c5_bond_pairs.append((o5_atoms[atom_list_idx][1], c5_atoms[atom_list_idx][1]))
        c5_c4_bond_pairs.append((c5_atoms[atom_list_idx][1], c4_atoms[atom_list_idx][1]))
        c4_c3_bond_pairs.append((c4_atoms[atom_list_idx][1], c3_atoms[atom_list_idx][1]))
        c3_o3_bond_pairs.append((c3_atoms[atom_list_idx][1], o3_atoms[atom_list_idx][1]))

    bonds = [p_o5_bond_pairs, o5_c5_bond_pairs, c5_c4_bond_pairs, c4_c3_bond_pairs, c3_o3_bond_pairs, o3_p_bond_pairs]
    bond_names = ["P-O5'", "O5'-C5'", "C5'-C4'", "C4'-C3'", "C3'-O3'", "O3'-P"]

    total_index = 0
    for bond, bond_name in zip(bonds, bond_names):
        for pair_idx, pair in enumerate(bond):
            output_file.write(f"[ {bond_name}.{pair_idx + 1} ]\n")
            output_file.write(f"{pair[0]} {pair[1]}\n")
            total_index += 1

    output_file.close()


if __name__ == "__main__":
    gro_file_name = sys.argv[1]
    structure_directory = sys.argv[2]
    determine_indices_end_to_end(gro_file_name, structure_directory, "TempEtE.ndx")
    determine_indices_base_distance(gro_file_name, structure_directory, "TempBaseDist.ndx")
    determine_indices_backbone_bonds(gro_file_name, structure_directory, "TempBackboneDist.ndx")
