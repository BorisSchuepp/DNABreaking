# Script for statistical analysis of breaking distributions
# Author: Boris N. Schüpp

import statistics
import scipy.stats as stats
import numpy as np
from PlottingNGS import *


# Bin the p-values into significacy categories
def determine_stars(pvalue_in):
    if pvalue_in < 0.0001:
        return 4
    elif pvalue_in < 0.001:
        return 3
    elif pvalue_in < 0.01:
        return 2
    elif pvalue_in < 0.05:
        return 1
    else:
        return 0


# Data read-in and preparation #
data_output = "../ProcessedData/704_AT_GC_400_1500_ATMi_Statistics/"
data_directory = "../ProcessedData/704_AT_GC_400_1500_ATMi/"
graphics_output = "../Graphics/704_AT_GC_400_1500_ATMi_Statistics/"
window = 30
file_names = ["704_1_Breaking_Dist_Full_Shifted.csv", "704_2_Breaking_Dist_Full_Shifted.csv",
              "704_3_Breaking_Dist_Full_Shifted.csv", "AT_1_Breaking_Dist_Full_Shifted.csv",
              "AT_2_Breaking_Dist_Full_Shifted.csv", "AT_3_Breaking_Dist_Full_Shifted.csv",
              "GC_1_Breaking_Dist_Full_Shifted.csv", "GC_2_Breaking_Dist_Full_Shifted.csv",
              "GC_3_Breaking_Dist_Full_Shifted.csv", "400_1_Breaking_Dist_Full_Shifted.csv",
              "400_2_Breaking_Dist_Full_Shifted.csv", "400_3_Breaking_Dist_Full_Shifted.csv",
              "1500_1_Breaking_Dist_Full_Shifted.csv", "1500_2_Breaking_Dist_Full_Shifted.csv",
              "1500_3_Breaking_Dist_Full_Shifted.csv",
              "ATMi0_1_Breaking_Dist_Full_Shifted.csv", "ATMi0_2_Breaking_Dist_Full_Shifted.csv",
              "ATMi0_3_Breaking_Dist_Full_Shifted.csv", "ATMi1_1_Breaking_Dist_Full_Shifted.csv",
              "ATMi1_2_Breaking_Dist_Full_Shifted.csv", "ATMi1_3_Breaking_Dist_Full_Shifted.csv",
              "ATMi2_1_Breaking_Dist_Full_Shifted.csv", "ATMi2_2_Breaking_Dist_Full_Shifted.csv",
              "ATMi2_3_Breaking_Dist_Full_Shifted.csv", "ATMi3_1_Breaking_Dist_Full_Shifted.csv",
              "ATMi3_2_Breaking_Dist_Full_Shifted.csv", "ATMi3_3_Breaking_Dist_Full_Shifted.csv",
              "ATMi4_1_Breaking_Dist_Full_Shifted.csv", "ATMi4_2_Breaking_Dist_Full_Shifted.csv",
              "ATMi4_3_Breaking_Dist_Full_Shifted.csv"]

# For each sequence combine the three replicas to one big dataset (aggregated data)
# Run this first with this set of prefixes, then swap for the other and run again
# prefixes = ["400", "1500", "704", "GC", "AT"]
# addtion = "704_AT_GC_1500_400"

prefixes = ["ATMi0", "ATMi1", "ATMi2", "ATMi3", "ATMi4"]
addtion = "ATMi"


breaking_data = DataUtilities.read_breaking_results([a for a in file_names if
                                                     a.startswith(tuple([f"{b}_" for b in prefixes]))], data_directory)
all_histogram_data_full = {i.split("#")[0]: breaking_data[i] for i in breaking_data.keys() if i.find("TOP") != -1}

# Restrict histogram to desired window (here = 30 bases in both direcitons)
all_histogram_data = {}
for i in all_histogram_data_full.keys():
    temp_dict = {}
    for index in all_histogram_data_full[i]:
        if -window <= index <= window:
            temp_dict[index] = all_histogram_data_full[i][index]
    all_histogram_data[i] = temp_dict

# Recreate full data set from histogram
all_list_data = {}
for i in all_histogram_data.keys():
    temp_list = []
    for index in all_histogram_data[i]:
        for k in range(0, int(all_histogram_data[i][index])):
            temp_list.append(index)
    all_list_data[i] = temp_list

name_translate = {'704_1': '704* - S1', '704_2': '704* - S2', '704_3': '704* - S3',
                  'AT_1': "AT* - S1", 'AT_2': "AT* - S2", 'AT_3': "AT* - S3",
                  'GC_1': "GC* - S1", 'GC_2': "GC* - S2", 'GC_3': "GC* - S3",
                  '400_1': "400* - S1", '400_2': "400* - S2", '400_3': "400* - S3",
                  '1500_1': "1500* - S1", '1500_2': "1500* - S2", '1500_3': "1500* - S3",
                  '704': "704*", "AT": "AT*", "GC": "GC*", "400": "400*", "1500": "1500*",
                  "HP0_1": "HP0 - S1", "HP0_2": "HP0 - S2", "HP1_1": "HP1 - S1", "HP1_2": "HP1 - S2",
                  "HP2_1": "HP2 - S1", "HP2_2": "HP2 - S2", "HP2_3": "HP2 -S3", "HP0": "HP0", "HP1": "HP1",
                  "HP2": "HP2", "ATMi0_1": "ATMi0 - S1", "ATMi0_2": "ATMi0 - S2", "ATMi0_3": "ATMi0 - S3",
                  "ATMi1_1": "ATMi1 - S1", "ATMi1_2": "ATMi1 - S2", "ATMi1_3": "ATMi1 - S3",
                  "ATMi2_1": "ATMi2 - S1", "ATMi2_2": "ATMi2 - S2", "ATMi2_3": "ATMi2 - S3",
                  "ATMi3_1": "ATMi3 - S1", "ATMi3_2": "ATMi3 - S2", "ATMi3_3": "ATMi3 - S3",
                  "ATMi4_1": "ATMi4 - S1", "ATMi4_2": "ATMi4 - S2", "ATMi4_3": "ATMi4 - S3",
                  "ATMi0": "ATMi0", "ATMi1": "ATMi1", "ATMi2": "ATMi2", "ATMi3": "ATMi3", "ATMi4": "ATMi4"}


aggregated_lists = {}
for prefix in prefixes:
    names_sample = [f"{prefix}_{i}" for i in range(1, 4)]
    aggregated_data_list = []
    for name_sample in names_sample:
        try:
            current_sample_histogram = all_histogram_data[name_sample]
            current_sample_data = all_list_data[name_sample]
            aggregated_data_list += current_sample_data
        except FileNotFoundError:
            print(f"WARNING: No sample {name_sample} found!")
    aggregated_lists[prefix] = aggregated_data_list

# Calculate mean and standard deviation of the breaking distributions
all_mu = {}
all_sigma = {}
order = [f"{sample}_{i}" for sample in prefixes for i in [1, 2, 3]]

for sample in order:
    mu = statistics.mean(all_list_data[sample])
    sigma = statistics.stdev(all_list_data[sample])
    all_mu[sample] = mu
    all_sigma[sample] = sigma

groups_dict = {sample: sample.split("_")[0] for sample in all_list_data.keys()}

DataUtilities.print_fit_parameters(all_mu, groups_dict,
                                   f"Data_Mu_All_{addtion}.csv",
                                   data_output)
DataUtilities.print_fit_parameters(all_sigma, groups_dict,
                                   f"Data_Sigma_All_{addtion}.csv",
                                   data_output)

plotter = PlottingNGS()

plotter.plot_fit_paramameter_scatter_single([f"Data_Mu_All_{addtion}.csv"],
                                            data_output, "Mu",
                                            f"{graphics_output}Mu_Shifted_Scatter_Plot_{addtion}.pdf")
plotter.plot_fit_paramameter_scatter_single([f"Data_Sigma_All_{addtion}.csv"],
                                            data_output, "Sigma",
                                            f"{graphics_output}Sigma_Scatter_Plot_{addtion}.pdf")

# Perform the Welch's t-test for all pairs of single samples

out_matrix_single_samples = np.zeros((len(all_list_data.keys()), len(all_list_data.keys())))
out_matrix_single_samples_values = np.zeros((len(all_list_data.keys()), len(all_list_data.keys())))

for i, name in enumerate(all_list_data.keys()):
    for j, name_2 in enumerate(all_list_data.keys()):
        pvalue = stats.ttest_ind(all_list_data[name], all_list_data[name_2], equal_var=False).pvalue
        out_matrix_single_samples_values[i][j] = pvalue
        out_matrix_single_samples[i][j] = determine_stars(pvalue)

DataUtilities.print_p_matrix(out_matrix_single_samples_values, [name_translate[i] for i in all_list_data.keys()],
                             f"Welch_Single_Exact_Values_{addtion}.csv",
                             data_output)
DataUtilities.print_p_matrix(out_matrix_single_samples, [name_translate[i] for i in all_list_data.keys()],
                             f"Welch_Single_Binned_Values_{addtion}.csv",
                             data_output)

# Create a plot of the p-value bins (red = n.s., orange = *, green = **)
plotter.plot_p_value_matrix(f"Welch_Single_Binned_Values_{addtion}.csv", data_output,
                            f"{graphics_output}Welch_Single_Binned_{addtion}.pdf")

# Perform the Welch's t-test for all pairs of aggregated samples

out_matrix_aggregated = np.zeros((len(aggregated_lists.keys()), len(aggregated_lists.keys())))
out_matrix_aggregated_values = np.zeros((len(aggregated_lists.keys()), len(aggregated_lists.keys())))

for i, name in enumerate(aggregated_lists.keys()):
    for j, name_2 in enumerate(aggregated_lists.keys()):
        pvalue = stats.ttest_ind(aggregated_lists[name], aggregated_lists[name_2], equal_var=False).pvalue
        out_matrix_aggregated_values[i][j] = pvalue
        out_matrix_aggregated[i][j] = determine_stars(pvalue)

DataUtilities.print_p_matrix(out_matrix_aggregated_values, [name_translate[i] for i in aggregated_lists.keys()],
                             f"Welch_Aggregated_Exact_Values_{addtion}.csv",
                             data_output)
DataUtilities.print_p_matrix(out_matrix_aggregated, [name_translate[i] for i in aggregated_lists.keys()],
                             f"Welch_Aggregated_Binned_Values_{addtion}.csv",
                             data_output)

# Create a plot of the p-value bins (red = n.s., orange = *, green = **)
plotter.plot_p_value_matrix(f"Welch_Aggregated_Binned_Values_{addtion}.csv", data_output,
                            f"{graphics_output}Welch_Aggregated_Binned_{addtion}.pdf")

# Perform the Levene's test for all pairs of single samples

out_matrix_single_samples_sigma = np.zeros((len(all_list_data.keys()), len(all_list_data.keys())))
out_matrix_single_samples_values_sigma = np.zeros((len(all_list_data.keys()), len(all_list_data.keys())))

for i, name in enumerate(all_list_data.keys()):
    for j, name_2 in enumerate(all_list_data.keys()):
        pvalue = stats.levene(all_list_data[name], all_list_data[name_2], center='mean').pvalue
        out_matrix_single_samples_values_sigma[i][j] = pvalue
        out_matrix_single_samples_sigma[i][j] = determine_stars(pvalue)

DataUtilities.print_p_matrix(out_matrix_single_samples_values_sigma, [name_translate[i] for i in all_list_data.keys()],
                             f"Levene_Single_Exact_Values_{addtion}.csv",
                             data_output)
DataUtilities.print_p_matrix(out_matrix_single_samples_sigma, [name_translate[i] for i in all_list_data.keys()],
                             f"Levene_Single_Binned_Values_{addtion}.csv",
                             data_output)

# Create a plot of the p-value bins (red = n.s., orange = *, green = **)
plotter.plot_p_value_matrix(f"Levene_Single_Binned_Values_{addtion}.csv", data_output,
                            f"{graphics_output}/Levene_Single_Binned_{addtion}.pdf")

# Perform the Levene's test for all pairs of aggregated samples

out_matrix_aggregated_sigma = np.zeros((len(aggregated_lists.keys()), len(aggregated_lists.keys())))
out_matrix_aggregated_values_sigma = np.zeros((len(aggregated_lists.keys()), len(aggregated_lists.keys())))

for i, name in enumerate(aggregated_lists.keys()):
    for j, name_2 in enumerate(aggregated_lists.keys()):
        pvalue = stats.levene(aggregated_lists[name], aggregated_lists[name_2], center='mean').pvalue
        out_matrix_aggregated_values_sigma[i][j] = pvalue
        out_matrix_aggregated_sigma[i][j] = determine_stars(pvalue)

DataUtilities.print_p_matrix(out_matrix_aggregated_values_sigma, [name_translate[i] for i in aggregated_lists.keys()],
                             f"Levene_Aggregated_Exact_Values_{addtion}.csv",
                             data_output)
DataUtilities.print_p_matrix(out_matrix_aggregated_sigma, [name_translate[i] for i in aggregated_lists.keys()],
                             f"Levene_Aggregated_Binned_Values_{addtion}.csv",
                             data_output)

# Create a plot of the p-value bins (red = n.s., orange = *, green = **)
plotter.plot_p_value_matrix(f"Levene_Aggregated_Binned_Values_{addtion}.csv", data_output,
                            f"{graphics_output}Levene_Aggregated_Binned_{addtion}.pdf")
