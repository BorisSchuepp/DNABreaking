# Nick creation script
# Author: Boris N. Schüpp, marked functions have been adapted from
# Benedikt Rennekamp (https://github.com/HITS-MBM/kimmdy) denoted as **

import numpy as np
from MolecularUtilities import find_d_coordinates


# Adapted from **
def get_data_from_file(filepath):
    file = open(filepath, 'r')
    data_all = []  # array with each entry corresponding to one (string) line
    data_array = []  # array of all lines with each subarray containing one value
    settings = []
    for line in file:
        data_all.append(line)
        line_array = np.asarray(line.split())
        data_array.append(line_array)
    file.close()
    return data_all, data_array


# Adapted from **
def store_linelist_to_file(data, filepath):
    file = open(filepath, "w")
    for line in data:
        file.write(line)
    file.close()


# Adapted from **
def modify_top(oldtop, newtop, breakpair):
    # function that cuts topology into the parts at the breakpair deleting all interactions
    data_all, data_array = get_data_from_file(oldtop)

    proper_set = False
    list_of_pairs = []  # pairs
    list_of_dihedrals = []
    number_of_atm = 0
    possible_lower_pairpartner = []
    possible_higher_pairpartner = []
    deleted = 0
    reached_pairs = False
    data_new = data_all[:]  # get an independent copy with slice
    # remove bond, angles and dihedrals where breakpair was involved
    # save pairs, dihedrals and their respective positions in new data set for next step
    pairs = []
    dihedrals = []
    for i in range(len(data_all)):
        if breakpair[0] in data_array[i][0:4] and breakpair[1] in data_array[i][0:4]:
            print("Deleted: " + data_new[i - deleted])  # shift -
            data_new.remove(data_new[i - deleted])
            deleted = deleted + 1

        if "[ bonds ]" in data_all[i]:
            bonds = i
            number_of_atms = data_array[bonds - 2][0]
        if "[ pairs ]" in data_all[i]:
            pairs = i
            reached_pairs = True
        if reached_pairs and i > pairs + 1:
            list_of_pairs.append(list(data_array[i][0:2]))
        if "[ angles ]" in data_all[i]:
            angles = i
            reached_pairs = False
            list_of_pairs = list_of_pairs[:-2]
        if "[ dihedrals ]" in data_all[i]:
            if not proper_set:
                dihedrals = i
                proper_set = True
            else:
                improper = i
                proper_set = False
        if proper_set and i > dihedrals + 1:
            list_of_dihedrals.append(data_array[i])

    list_of_dihedrals = list_of_dihedrals[:-1]
    deleted_pairs = 0

    # go through all dihedrals and thus find pairs to be deleted if breakpair is in middle of dihedral
    pairs_to_be_deleted = []
    for j in range(len(list_of_dihedrals)):
        if (breakpair[0] in list_of_dihedrals[j][0:4]) and (breakpair[1] in list_of_dihedrals[j][0:4]):
            if float(list_of_dihedrals[j][0]) < float(list_of_dihedrals[j][3]):  # pairs are sorted
                pair = [list_of_dihedrals[j][0], list_of_dihedrals[j][3]]
            else:
                pair = [list_of_dihedrals[j][3], list_of_dihedrals[j][0]]
            pairs_to_be_deleted.append(pair)
    for k in range(len(list_of_pairs)):
        if list_of_pairs[k] in pairs_to_be_deleted:
            print("Deleted: " + data_new[k + pairs - deleted_pairs + 1])
            data_new.remove(data_new[k + pairs - deleted_pairs + 1])
            deleted_pairs += 1

    store_linelist_to_file(data_new, newtop)


def read_in_residues(pdb_file_in):
    pdb_file = open(pdb_file_in)

    residue_list_a = []
    residue_list_b = []
    atom_list = []
    for line in pdb_file:
        if line.find("ATOM") == -1:
            continue

        values = [line[0:4], line[6:11], line[12:16], line[16], line[17:20], line[21],
                  line[22:26], line[27], line[30:38], line[38:46], line[46:54], line[54:60], line[60:66]]
        atom_list.append(values)
    pdb_file.close()
    number_of_residues = max([int(i[6]) for i in atom_list])
    for strand, strand_list in zip(["A", "B"], [residue_list_a, residue_list_b]):
        for i in range(1, number_of_residues + 1):
            residue_i = [atom for atom in atom_list if int(atom[6]) == i and atom[5] == strand]
            strand_list.append(residue_i)
    residue_list_b = list(reversed(residue_list_b))

    return residue_list_a, residue_list_b, number_of_residues


def print_middle_sequence(residue_list_a, residue_list_b, number_of_bases, display_range):
    residue_names_a = [residue[0][4] for residue in residue_list_a]
    residue_names_b = [residue[0][4] for residue in residue_list_b]
    bases_a = [name[1] for name in residue_names_a]
    bases_b = [name[1] for name in residue_names_b]

    displayed_bases_a = bases_a[int(number_of_bases / 2) - display_range:int(number_of_bases / 2) + display_range]
    displayed_bases_b = bases_b[int(number_of_bases / 2) - display_range:int(number_of_bases / 2) + display_range]
    display_indices = [str(i) for i in range(int(number_of_bases / 2) - display_range + 1,
                                             int(number_of_bases / 2) + display_range + 1)]
    display_index_string = "      " + " ".join(display_indices)
    display_string_a = "A 5'  " + "  ".join(displayed_bases_a) + "  3'"
    display_string_b = "B 3'  " + "  ".join(displayed_bases_b) + "  5'"
    print("The sequence in the middle is:")
    print(display_index_string)
    print(display_string_a)
    print(display_string_b)


def remove_phosphorous(information, strand_a_in, strand_b_in):
    strand = information[0]
    first_base = int(information[1])
    second_base = int(information[2])

    if strand == "A":
        remove_base = second_base - 1
        remove_strand = strand_a_in
    else:
        remove_base = first_base - 1
        remove_strand = strand_b_in
    new_res = [a for a in remove_strand[remove_base] if a[2].strip(" ") not in ["P", "O1P", "O2P"]]
    remove_strand[remove_base] = new_res


def add_hydrogen(information, strand_a_in, strand_b_in):
    strand = information[0]
    first_base = int(information[1])
    second_base = int(information[2])

    if strand == "A":
        three_prime_add = first_base - 1
        five_prime_add = second_base - 1
        add_strand = strand_a_in
    else:
        three_prime_add = second_base - 1
        five_prime_add = first_base - 1
        add_strand = strand_b_in

    c_four_prime = [a for a in add_strand[five_prime_add] if a[2].strip(" ") == "C4'"][0]
    c_five_prime = [a for a in add_strand[five_prime_add] if a[2].strip(" ") == "C5'"][0]
    o_five_prime = [a for a in add_strand[five_prime_add] if a[2].strip(" ") == "O5'"][0]
    coords_c_four_prime = [float(i.strip(" ")) for i in c_four_prime[8:11]]
    coords_c_five_prime = [float(i.strip(" ")) for i in c_five_prime[8:11]]
    coords_o_five_prime = [float(i.strip(" ")) for i in o_five_prime[8:11]]

    hydrogen_five_prime_coords = find_d_coordinates(coords_c_four_prime, coords_c_five_prime,
                                                    coords_o_five_prime, BOND_DISTANCE,
                                                    BOND_ANGLE, DIHEDRAL_ANGLE_FIVE_PRIME)
    hydrogen_five_prime_coords_padded = []
    for a in hydrogen_five_prime_coords:
        a = round(a, 3)
        a_padded = f"{a:>{8}}"
        hydrogen_five_prime_coords_padded.append(a_padded)

    h_five_prime = ["ATOM", str(10000), "HO5'", " ",
                    c_four_prime[4], c_four_prime[5], c_four_prime[6], c_four_prime[7],
                    hydrogen_five_prime_coords_padded[0], hydrogen_five_prime_coords_padded[1],
                    hydrogen_five_prime_coords_padded[2], c_four_prime[11], c_four_prime[12]]

    add_strand[five_prime_add].append(h_five_prime)

    c_two_prime = [a for a in add_strand[three_prime_add] if a[2].strip(" ") == "C2'"][0]
    c_three_prime = [a for a in add_strand[three_prime_add] if a[2].strip(" ") == "C3'"][0]
    o_three_prime = [a for a in add_strand[three_prime_add] if a[2].strip(" ") == "O3'"][0]
    coords_c_two_prime = [float(i.strip(" ")) for i in c_two_prime[8:11]]
    coords_c_three_prime = [float(i.strip(" ")) for i in c_three_prime[8:11]]
    coords_o_three_prime = [float(i.strip(" ")) for i in o_three_prime[8:11]]

    hydrogen_three_prime_coords = find_d_coordinates(coords_c_two_prime, coords_c_three_prime,
                                                     coords_o_three_prime, BOND_DISTANCE,
                                                     BOND_ANGLE, DIHEDRAL_ANGLE_THREE_PRIME)
    hydrogen_three_prime_coords_padded = []
    for a in hydrogen_three_prime_coords:
        a = round(a, 3)
        a_padded = f"{a:>{8}}"
        hydrogen_three_prime_coords_padded.append(a_padded)

    h_three_prime = ["ATOM", str(10000), "HO3'", " ",
                     c_two_prime[4], c_two_prime[5], c_two_prime[6], c_two_prime[7],
                     hydrogen_three_prime_coords_padded[0], hydrogen_three_prime_coords_padded[1],
                     hydrogen_three_prime_coords_padded[2], c_two_prime[11], c_two_prime[12]]

    add_strand[three_prime_add].append(h_three_prime)


def rename_residues(information, strand_a_in, strand_b_in):
    strand = information[0]
    first_base = int(information[1])
    second_base = int(information[2])

    if strand == "A":
        new_three_prime = first_base - 1
        new_five_prime = second_base - 1
        rename_strand = strand_a_in
    else:
        new_three_prime = second_base - 1
        new_five_prime = first_base - 1
        rename_strand = strand_b_in

    old_three_prime_res = rename_strand[new_three_prime]
    for a in old_three_prime_res:
        a[4] = a[4].strip(" ") + "3"

    old_three_prime_res = rename_strand[new_five_prime]
    for a in old_three_prime_res:
        a[4] = a[4].strip(" ") + "5"


def adjust_numbering(*strands):
    running_number = 1
    for strand in strands:
        for res in strand:
            for atom in res:
                atom[1] = f"{running_number:>{5}}"
                running_number += 1


def split_strand(information, strand_a_in, strand_b_in):
    strand = information[0]
    first_base = int(information[1])

    if strand == "A":
        last_base = first_base - 1
        to_be_split_strand = strand_a_in
    else:
        last_base = first_base - 1
        to_be_split_strand = strand_b_in

    new_strand_one = []
    new_strand_two = []

    for residue_num, residue in enumerate(to_be_split_strand):
        if residue_num <= last_base:
            new_strand_one.append(residue)
        else:
            new_strand_two.append(residue)

    if strand == "A":
        return new_strand_one, new_strand_two, list(reversed(strand_b_in))
    else:
        return strand_a_in, list(reversed(new_strand_two)), list(reversed(new_strand_one))


def rename_strands_residues(*strands):
    strand_ids = ["A", "B", "C", "D", "E", "F", "G"]
    for strand_idx, strand in enumerate(strands):
        strand_id = strand_ids[strand_idx]
        for residue_id, residue in enumerate(strand):
            for atom in residue:
                atom[5] = strand_id
                atom[6] = f"{residue_id + 1:>{4}}"


def write_to_file(*strands, out_file_name_in):
    out_file = open(out_file_name_in, "w")

    for strand in strands:
        for res in strand:
            for atom in res:
                out_string = atom[0] + "  " + atom[1] + " " + atom[2] + atom[3] + atom[4] + " " + atom[5] + atom[6] + \
                             atom[7] + "   " + atom[8] + atom[9] + atom[10] + atom[11] + atom[12] + "\n"
                out_file.write(out_string)
        out_file.write("TER\n")


DIHEDRAL_ANGLE_FIVE_PRIME = 65.0  # Determined from preamtive studies, in Degrees
DIHEDRAL_ANGLE_THREE_PRIME = 179.0  # Determined from preamtive studies, in Degrees
BOND_ANGLE = 108.5  # In Degrees
BOND_DISTANCE = 0.96  # In Angstrom
STRUCTURE_BASE_PATH = "../StructureGeneration/TestStructures/"

if __name__ == "__main__":

    nick_type = input("Which nick type is desired? 1 = Nick Type B, 2 = "
                      "Nick Type C (see Supplementary Information). \n")

    if int(nick_type) == 1:
        print("Selected nick type B.\n")
        file_name = input("Please enter the pdb file name \n")
        pdb_file_path = STRUCTURE_BASE_PATH + file_name
        residues_strand_a, residues_strand_b, num_of_bases = read_in_residues(pdb_file_path)
        print_middle_sequence(residues_strand_a, residues_strand_b, num_of_bases, 5)
        nick_information = input("Please select on which strand, and between which two "
                                 "indices a nick should be placed. Enter as A or B and "
                                 "two integers (example A 23 24).\n")
        nick_information = [a for a in nick_information.split(" ")]

        remove_phosphorous(nick_information, residues_strand_a, residues_strand_b)
        add_hydrogen(nick_information, residues_strand_a, residues_strand_b)
        rename_residues(nick_information, residues_strand_a, residues_strand_b)
        chain_a, chain_b, chain_c = split_strand(nick_information, residues_strand_a, residues_strand_b)
        rename_strands_residues(chain_a, chain_b, chain_c)
        adjust_numbering(chain_a, chain_b, chain_c)
        out_file_name = STRUCTURE_BASE_PATH + file_name.strip(".pdb") + "_modified.pdb"
        write_to_file(chain_a, chain_b, chain_c, out_file_name_in=out_file_name)
        print(f"Success! The modified structure has been written to {out_file_name}.")

    else:
        print("Selected nick type C.\n")
        input_file_name = STRUCTURE_BASE_PATH + input("Please enter the itp file name. \n")
        out_file_name = input_file_name[0:len(input_file_name)-4] + "_modified.itp"

        break_input = input("Please enter two atom numbers (integers) between the nick should be placed, "
                            "seperated by a space (example 1586 1587) \n")

        first_index = break_input.split(" ")[0]
        second_index = break_input.split(" ")[1]

        modify_top(input_file_name, out_file_name, [first_index, second_index])
        print("Topology changed. Make sure to change the file name in your .top file!\n")
