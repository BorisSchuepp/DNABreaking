How to create an input file for the matching algorithm:
Provide the following fields: data_directory, graphics_directory, data_export_directory, sample_names, sample_sequences, nicking_pos, sample_files, secondary_sequences, secondary_nicking_pos
-> data_directory: Directory path where the raw data is saved
-> graphics_directory: Directory where graphics are will be saved to
-> data_export_directory: Directory where processed data will be saved to
-> sample_names: Names of you samples
-> sample_sequences: DNA sequence of you DNA (5'->3' direction)
-> nicking_pos: Options (top, bot, none) on which strand is the nick (top = the strand you provided under sample_sequences)
-> sample_files: File names of the data samples within the data_directory
-> secondary_sequences, secondary_nicking_pos are only relevant for single read analysis, if your sequence is doublestranded (nicked/non-nicked) you don't need this can can put "c" for each sample
Write each of the fields in a separate line in the shape key value_1 value_2 value_3 ... value_n only seperated by spaces (don't use any extra spaces in file names etc.)
All directory names need to end with "/". For paired-end sequencing make sure that the paired files only differ in "R1" and "R2"!