import os
import sys 


def combine_end_to_end_distance_files(directory_in):
    out_file_frac = open(f"{directory_in}/EndToEndDistance_fractional_all.csv", "w")
    out_file_absolute = open(f"{directory_in}/EndToEndDistance_all.csv", "w")
    for filename in os.listdir(directory_in):
        if filename.find("all") == -1 and filename.find("EndToEndDistance") != -1:
            full_path = os.path.join(directory_in, filename)
            cur_file = open(full_path)
            for line in cur_file:
                line = line.strip("\n") + "\n"
                if full_path.find("frac") != -1:
                    out_file_frac.write(line)
                else:
                    out_file_absolute.write(line)
            cur_file.close()
    out_file_frac.close()
    out_file_absolute.close()


def combine_base_distance_files(directory_in):
    out_file = open(f"{directory_in}/BaseDistance_all.csv", "w")
    for filename in os.listdir(directory_in):
        if filename.find("all") == -1 and filename.find("BaseDistance") != -1:
            full_path = os.path.join(directory_in, filename)
            cur_file = open(full_path)
            for line in cur_file:
                line = line.strip("\n") + "\n"
                out_file.write(line)
            cur_file.close()
    out_file.close()


def combine_backbone_bond_distance_files(directory_in):
    out_file = open(f"{directory_in}/BackboneBondDistance_all.csv", "w")
    for filename in os.listdir(directory_in):
        if filename.find("all") == -1 and filename.find("BackboneBondDistances") != -1:
            full_path = os.path.join(directory_in, filename)
            cur_file = open(full_path)
            for line in cur_file:
                line = line.strip("\n") + "\n"
                out_file.write(line)
            cur_file.close()
    out_file.close()


if __name__ == "__main__":
    out_put_dir = sys.argv[1]
    combine_end_to_end_distance_files(out_put_dir)
    combine_base_distance_files(out_put_dir)
    combine_backbone_bond_distance_files(out_put_dir)
