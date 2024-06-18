# Configuration handeler for input configurations of breaking analysis via NGS data
# Author: Boris N. Schuepp

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


# Check if a filename contains chars that prevent the saving of files with this name, returns a bool
def contains_invalid_chars(filename):
    invalid_chars = set('<>:"/\\|?*')
    return any((char in invalid_chars) for char in filename) or filename.endswith('.') or filename.endswith(' ')


# Replace invalid charactes in a given filename, returns a sanitized file name as a string
def sanitize_filename(filename):
    invalid_chars = '<>:"/\\|?*'
    sanitized_filename = ''.join(char for char in filename if char not in invalid_chars)
    sanitized_filename = sanitized_filename.rstrip('. ')  # Remove trailing dots and spaces
    return sanitized_filename


# Parse and check the configuration file provided. Return success parameter and a list of extracted fields in the shape
# success param (bool), [sample file names (list of tuples of str), sample names checked (list of str),
# sample sequences (list of str), sample_nicks (list of int), data_dir (str), graphics_dir (str), data_export_dir (str)]
def read_config(config_file_path_in):
    valid_keys = ["data_directory", "graphics_directory", "data_export_directory",
                  "sample_sequences", "nicking_pos", "sample_files", "sample_names",
                  "secondary_sequences", "secondary_nicking_pos"]
    config_file = open(config_file_path_in, "r")

    config_data_dict = {}

    # Read data
    for line in config_file:
        # Ignore comments
        if line.startswith("#"):
            continue
        current_key = line.strip("\n").split(" ")[0]
        if current_key not in valid_keys:
            print(f"ERROR: Invalid key {current_key} found. Please try again with different configuration file!")
            return False, []

        config_data_dict[current_key] = line.strip("\n").split(" ")[1:]

    for key in valid_keys:
        if key not in config_data_dict.keys():
            print(f"ERROR: Required key {key} not found. Please try again with different configuration file!")
            return False, []

    # Parse data
    sample_file_names = config_data_dict["sample_files"]
    sample_names = config_data_dict["sample_names"]
    sample_sequences = config_data_dict["sample_sequences"]
    sample_sequences_second = config_data_dict["secondary_sequences"]
    sample_sequences_second_final = []
    for sequence_sec, sequence in zip(sample_sequences_second, sample_sequences):
        if sequence_sec == "c":
            sample_sequences_second_final.append(get_complementary_strand(sequence))
        else:
            sample_sequences_second_final.append(sequence_sec)

    try:
        data_dir = config_data_dict["data_directory"][0]
        graphics_dir = config_data_dict["graphics_directory"][0]
        data_export_dir = config_data_dict["data_export_directory"][0]
    except IndexError:
        print("ERROR: Seems like one of the required directories is empty. "
              "Please try again with different configuration file!")
        return False, []

    # Create paired read files (only differ in R1/R2)
    sample_file_names_paired_reads = []
    for sample_file_name in sample_file_names:
        s1 = sample_file_name
        s2 = sample_file_name.replace("R1", "R2")
        sample_file_names_paired_reads.append((s1, s2))

    # Convert possible nicking configurations to numerical values
    nicking_pos_translation = {"top": 0, "bot": 1, "none": 2}
    sample_nicks = []
    for nicking_pos, nicking_pos_secondary in zip(config_data_dict["nicking_pos"], config_data_dict["secondary_nicking_pos"]):
        strand_a = nicking_pos.split("_")[0]
        strand_b = nicking_pos_secondary.split("_")[0]
        if (strand_a not in nicking_pos_translation.keys() or
                (strand_b not in nicking_pos_translation.keys() and strand_b != "c")):
            print(f"ERROR: Nicking position {nicking_pos} is not valid. "
                  f"Please try again with different configuration file!")
            return False, []
        else:
            try:
                location_a = round(float(nicking_pos.split("_")[1]),1)
                if strand_b != "c":
                    location_b = round(float(nicking_pos_secondary.split("_")[1]),1)
                else:
                    location_b = location_a
                    if strand_a == "top":
                        strand_b = "bot"
                    else:
                        strand_b = "top"
                sample_nicks.append((nicking_pos_translation[strand_a], location_a,
                                     nicking_pos_translation[strand_b], location_b))
            except ValueError:
                print(f"ERROR: Nicking position {nicking_pos} is not valid. "
                      f"Please try again with different configuration file!")

    # Check if all sample fields have the same number of values
    if len(set([len(i) for i in [sample_file_names, sample_names,
                                 sample_sequences, sample_sequences_second_final, sample_nicks]])) != 1:
        print(f"ERROR: Unequal number of samples across required parameters. "
              f"Please try again with different configuration file!")
        return False, []

    # Sanitize sample names
    samples_names_checked = []
    for sample_name in sample_names:
        if contains_invalid_chars(sample_name):
            new_sample_name = sanitize_filename(sample_name)
            print(f"WARNING: Seems like the sample name {sample_name} contains characters that will "
                  f"prevent the creation of output files. Name was changed to {new_sample_name}.")
            samples_names_checked.append(new_sample_name)
        else:
            samples_names_checked.append(sample_name)

    return True, [sample_file_names_paired_reads, samples_names_checked, sample_sequences, sample_sequences_second_final,
                  sample_nicks, data_dir, graphics_dir, data_export_dir]
