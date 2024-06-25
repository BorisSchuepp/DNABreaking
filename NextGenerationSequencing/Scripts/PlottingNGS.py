# Utilities for plotting results from the NGS analysis section
# Author: Boris N. Schüpp

import matplotlib.pyplot as plt
import DataUtilities
import Fitting
import matplotlib.colors as mcolors

from matplotlib.colors import LinearSegmentedColormap

class PlottingNGS:

    def __init__(self):
        self.font_size = 20
        self.font = "Roboto"
        self.font_weight = "bold"
        self.marker_size = 400
        self.colors = {"AT": (0.9450980392156862, 0.7019607843137254, 0.25882352941176473, 1.0),
                       "GC": (0.19215686274509805, 0.49411764705882355, 0.7607843137254902, 1.0),
                       "704": (0.27058823529411763, 0.27058823529411763, 0.27058823529411763, 1.0),
                       "1500": (0.35294117647058826, 0.6705882352941176, 0.27058823529411763, 1.0),
                       "400": (0.6745098039215687, 0.8313725490196079, 0.6352941176470588, 1.0),
                       "BlackOpaque": (0, 0, 0, 0.8),
                       "BlackFaded": (0, 0, 0, 0.5),
                       "DarkGrey": (0.43921569, 0.43921569, 0.44313725, 1),
                       "LightGrey": (0.85490196, 0.85490196, 0.85882352, 1),
                       "HP0": "green",
                       "HP1": "orange",
                       "HP2": "red",
                       "ATMi0": (1, 0, 0, 1),
                       "ATMi1": (1, 0.25, 0, 1),
                       "ATMi2": (1, 0.5, 0, 1),
                       "ATMi3": (1, 0.75, 0, 1),
                       "ATMi4": (1, 1, 0, 1)
                       }
        plt.rcParams['mathtext.fontset'] = 'custom'
        plt.rcParams['mathtext.it'] = f'{self.font}:italic'
        plt.rcParams['mathtext.bf'] = f'{self.font}:italic:bold'
        plt.rcParams['axes.linewidth'] = 2

    @staticmethod
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

    # Plot the filtering results, optionally provide a directory where the graphics are saved to
    def plot_filtering_results(self, file_names_in, directory_in, save_path=None):
        filtering_results_dict = DataUtilities.read_filtering_results(file_names_in, directory_in)
        plt.figure(figsize=(9, 9))

        bar_width = 0.2
        indices = []
        mult_index = 0

        # Obtain correct position for bar plots
        for sample, filtering_result in filtering_results_dict.items():
            current_indices = [bar_width * i + mult_index * (0.8 + bar_width * len(filtering_result.keys()))
                               for i in range(0, len(filtering_result.keys()))]
            mult_index += 1
            indices.append(current_indices)

        keys = list(filtering_results_dict[list(filtering_results_dict.keys())[0]].keys())

        # Plot all bars at the respective spots
        for cur_key_id, cur_key in enumerate(keys):
            cur_values = [filtering_results_dict[sample_i][cur_key] for sample_i in filtering_results_dict.keys()]
            cur_indices = [indices[i][cur_key_id] for i in range(len(indices))]
            plt.bar(cur_indices, cur_values, bar_width, label=cur_key, edgecolor='black')

        legend_properties = {'weight': self.font_weight, 'family': self.font, 'size': f"{self.font_size - 5}"}
        plt.legend(prop=legend_properties)

        x_tick_list = [0.5 + 2 * i for i in range(0, len(filtering_results_dict.keys()))]
        plt.xticks(x_tick_list, [i.split("_")[0] + "*" + " Sample " + i.split("_")[1] for i in list(filtering_results_dict.keys())],
                   font=self.font, size=self.font_size, weight=self.font_weight)

        ax = plt.gca()
        ax.tick_params(axis="x", which='both', length=0)

        y_ticks = ax.get_yticks()
        plt.yticks(y_ticks, [str(round(i, 1)) for i in y_ticks],
                   font=self.font, size=self.font_size, weight=self.font_weight)

        plt.ylabel("Fraction of single reads",
                   font=self.font, size=self.font_size + 10, weight=self.font_weight)
        plt.tick_params(width=2)

        if save_path:
            plt.savefig(save_path)
        else:
            plt.show()
            plt.clf()
        plt.close()

    # Plot the length distribution from a provided csv file path, optionally save the pdf image to a file path
    def plot_lengths(self, file_name_in, directory_in, save_path=None):
        fragment_length_dict = DataUtilities.read_lengths(file_name_in, directory_in)
        plt.figure(figsize=(9, 9))
        lengths = sorted(list(fragment_length_dict.keys()))
        fraction_of_reads = [fragment_length_dict[i] for i in lengths]

        # Determine the largest consecutive group of fragment lengths with each of them above 0.1% of total reads
        cutoff = 0.001
        cutoff_range = []
        for idx, fraction_of_read in enumerate(fraction_of_reads):
            if fraction_of_read > cutoff:
                cutoff_range.append(idx)

        groups = []
        cur_group = []
        for idx in range(0, len(cutoff_range) - 1):
            if len(cur_group) == 0:
                cur_group.append(cutoff_range[idx])
            if cutoff_range[idx + 1] != cutoff_range[idx] + 1:
                groups.append(cur_group)
                cur_group = []
            else:
                cur_group.append(cutoff_range[idx + 1])

        if len(cur_group) > 0:
            groups.append(cur_group)

        max_group = groups[[len(i) for i in groups].index(max([len(i) for i in groups]))]

        plt.xlim([lengths[max_group[0]], lengths[max_group[-1]]])

        plt.bar(lengths, fraction_of_reads)
        plt.ylabel("Fraction of reads",
                   font=self.font, size=self.font_size + 10, weight=self.font_weight)
        plt.xlabel("Fragment length [bases]",
                   font=self.font, size=self.font_size + 10, weight=self.font_weight)
        plt.tick_params(width=2)
        ax = plt.gca()

        y_ticks = ax.get_yticks()
        x_ticks = ax.get_xticks()

        plt.yticks(y_ticks, [str(round(i, 6)) for i in y_ticks],
                   font=self.font, size=self.font_size, weight="bold")
        plt.xticks(x_ticks, [str(int(i)) for i in x_ticks],
                   font=self.font, size=self.font_size, weight="bold")

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
        else:
            plt.show()
            plt.clf()
        plt.close()

    # Plot the basecounts from a list of provided csv file paths, optionally save the pdf image to a file path
    def plot_basecounts(self, file_names_in, directory_in, save_path=None):
        basecount_results = DataUtilities.read_basecount_results(file_names_in, directory_in)
        plt.figure(figsize=(9, 9))
        x_max = -1
        for sample_name in basecount_results.keys():
            cur_sample = basecount_results[sample_name]
            current_x = sorted(list(cur_sample.keys()))
            current_y = [cur_sample[i] for i in current_x]
            if sample_name.find("_") != -1:
                label = sample_name.split("_")[0] + "*" + " Sample " + sample_name.split("_")[1]
            else:
                label = sample_name
            plt.plot(current_x, current_y, label=label, linewidth=2.5)
            plt.xlim(current_x[0], current_x[-1])
            x_max = current_x[-1]

        plt.ylabel("Fractional basecounts",
                   font=self.font, size=self.font_size + 10, weight=self.font_weight)
        plt.xlabel("Base index",
                   font=self.font, size=self.font_size + 10, weight=self.font_weight)
        plt.tick_params(width=2)
        ax = plt.gca()

        y_ticks = [i for i in ax.get_yticks() if i >= 0]
        x_ticks = [i for i in ax.get_xticks() if i <= x_max]

        plt.yticks(y_ticks, [str(round(i, 1)) for i in y_ticks],
                   font=self.font, size=self.font_size, weight="bold")
        plt.xticks(x_ticks, [str(int(i)) for i in x_ticks],
                   font=self.font, size=self.font_size, weight="bold")

        legend_properties = {'weight': self.font_weight, 'family': self.font, 'size': f"{self.font_size - 5}"}
        plt.legend(prop=legend_properties)

        if save_path:
            plt.savefig(save_path)
        else:
            plt.show()
            plt.clf()
        plt.close()

    def plot_single_breaking_bases(self, breaking_data, lims, save_path=None):
        plt.figure(figsize=(9, 9))

        for color, base_values, base in zip(["blue", "red", "orange", "green"],
                                            [breaking_data[i] for i in breaking_data.keys()],
                                            list(breaking_data.keys())):
            plt.bar([x[0] for x in base_values], [y[1] for y in base_values], color=color, label=base)

        plt.xlim(lims[0], lims[1])
        plt.ylim(0, 2)
        plt.xlabel("Fragment break index $\\mathbf{j}$", size=self.font_size + 10, font=self.font,
                   weight=self.font_weight)
        plt.ylabel("Relative breaking rate", size=self.font_size + 10, font=self.font, weight=self.font_weight)
        legend_properties = {'weight': self.font_weight, 'family': self.font, 'size': f"{self.font_size - 5}"}
        plt.legend(prop=legend_properties)
        plt.axhline(1, color='black', linewidth=2.5)
        plt.xticks(weight=self.font_weight, fontname=self.font, size=self.font_size)
        plt.yticks(weight=self.font_weight, fontname=self.font, size=self.font_size)
        plt.tick_params(width=2)

        if save_path:
            plt.savefig(save_path)
        else:
            plt.show()
            plt.clf()
        plt.close()

    def plot_single_breaking_data(self, breaking_data, top, gauss, shift, average, centroid=352,
                                  color=(0, 0, 0, 1),
                                  save_path=None):
        # Select desired subset
        data = []
        data_er = []
        if top and average:
            data = [breaking_data[i] for i in breaking_data.keys() if i.find("TOP") != -1 and i.find("AVERAGE") != -1][
                0]
            data_er = [breaking_data[i] for i in breaking_data.keys() if i.find("TOP") != -1 and i.find("STD") != -1][0]
        if not top and average:
            data = [breaking_data[i] for i in breaking_data.keys() if i.find("BOT") != -1 and i.find("AVERAGE") != -1][
                0]
            data_er = [breaking_data[i] for i in breaking_data.keys() if i.find("BOT") != -1 and i.find("STD") != -1][0]

        if top and not average:
            data = [breaking_data[i] for i in breaking_data.keys() if i.find("TOP") != -1][0]
        if not top and not average:
            data = [breaking_data[i] for i in breaking_data.keys() if i.find("BOT") != -1][0]
        if len(data) == 0:
            print("Data is not correct, please try again")
        data = {i: data[i] for i in data.keys() if data[i] != 0}
        plt.figure(figsize=(9, 9))
        x_data = list(data.keys())
        y_data = [data[i] for i in x_data]

        if average:
            plt.ylim(0, max(y_data) + max([data_er[i] for i in x_data]))
        else:
            plt.ylim(0, max(y_data) + 0.0001)

        if average:
            plt.bar(x_data, y_data, yerr=[data_er[i] for i in x_data], capsize=5, color=color,
                    edgecolor='black', error_kw={'ecolor': self.colors["BlackOpaque"]})
        else:
            plt.bar(x_data, y_data, color="blue")  # , edgecolor="black")

        if gauss == 0:
            pass
        elif gauss == 1:
            mu, sigma, x_gauss, y_gauss = Fitting.fit_gauss(data, True)
            plt.plot(x_gauss, y_gauss, color="red", linewidth=3.5)
        else:
            mu_1, sigma_1, a_1, mu_2, sigma_2, a_2, x_gauss, y_gauss = Fitting.fit_bimodal(data, True)
            print("Mu left:", mu_1, "Sigma left:", sigma_1, "Mu right:", mu_2, "Sigma right:", sigma_2,
                  "Left to right ratio:", a_1 / a_2)
            plt.plot(x_gauss, y_gauss, color="red", linewidth=3.5)

        plt.ylabel("Fraction of reads", size=self.font_size + 10, font=self.font, weight=self.font_weight)
        if top:
            if shift:
                plt.xlabel("Relative fragment break index $\\mathbf{j}$",
                           size=self.font_size + 10, font=self.font, weight=self.font_weight)
            else:
                plt.xlabel("Fragment break index $\\mathbf{j}$",
                           size=self.font_size + 10, font=self.font, weight=self.font_weight)
        else:
            if shift:
                plt.xlabel("Relative fragment break index $\\mathbf{i}$",
                           size=self.font_size + 10, font=self.font, weight=self.font_weight)
            else:
                plt.xlabel("Fragment break index $\\mathbf{i}$",
                           size=self.font_size + 10, font=self.font, weight=self.font_weight)

        plt.xlim(min(x_data), max(x_data))
        if not shift:
            x_ticks = [i for i in range(int(min(x_data)), int(max(x_data))) if i % 5 == centroid % 5]
            x_tick_labels = [str(round(int(i))) for i in x_ticks]
            # x_tick_labels[x_tick_labels.index((str(float(centroid))))] = str(centroid) + "*"
            plt.xticks(x_ticks, x_tick_labels, weight=self.font_weight, fontname=self.font, size=self.font_size)
        else:
            x_ticks = [i for i in range(int(min(x_data)), int(max(x_data))) if i % 10 == 0]
            x_tick_labels = []
            for x_tick in x_ticks:
                if x_tick > 0:
                    x_tick_labels.append("+" + str(int(x_tick)))
                else:
                    x_tick_labels.append(str(int(x_tick)))
            plt.xticks(x_ticks, x_tick_labels, weight=self.font_weight, fontname=self.font, size=self.font_size)

        plt.yticks(weight="bold", fontname=self.font, size=self.font_size)
        plt.tick_params(width=2)
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
        else:
            plt.show()
            plt.clf()
        plt.close()

    # Plot a single breaking diagram, provide if it is on the upper strand (the one provided in the configuration) via
    # parameter top, provide if a Gaussian fit should be performed (gauss), if you want to display a shilfed coordinate
    # system (shift) and if an averages is plotted (average). Optionally provide a color for the plpt and a path to
    # save it to
    def plot_single_breaking(self, file_names_in, directory_in, top, gauss, shift, average, centroid=352,
                             color=(0, 0, 0, 1),
                             save_path=None):

        breaking_data = DataUtilities.read_breaking_results(file_names_in, directory_in)
        self.plot_single_breaking_data(breaking_data, top, gauss, shift, average, centroid, color, save_path)

    # Plot three breaking distribtions side-by-side, provide the files and the directory, indicate if this is the
    # upper strand and if a gaussian fit should be plotted, for different nick positions provide a centriod and
    # optionally provide a path to save the figure to
    def plot_three_breaking(self, file_names_in, directory_in, top, gauss, centroid=352, save_path=None):

        breaking_data = DataUtilities.read_breaking_results(file_names_in, directory_in)

        if top:
            data = [breaking_data[i] for i in breaking_data.keys() if i.find("TOP") != -1]
        else:
            data = [breaking_data[i] for i in breaking_data.keys() if i.find("BOT") != -1]

        if len(data) == 0:
            data = [i for i in breaking_data.values()]

        total_max = max([max(list(data[i].values())) for i in range(0, 3)])
        fig, axs = plt.subplots(1, 3, figsize=(20, 5))

        for i in range(0, 3):
            axis = axs[i]
            plt.sca(axis)
            plt.yticks(weight="bold", fontname=self.font, size=self.font_size-5)
            plt.tick_params(width=2)

            fraction_of_reads = list(data[i].values())
            base_indices = list(data[i].keys())

            axis.bar(base_indices, fraction_of_reads)

            axis.set_ylabel("Fraction of reads", size=self.font_size + 10, font=self.font, weight=self.font_weight)
            if top:
                axis.set_xlabel("Fragment break index $\\mathbf{j}$",
                                size=self.font_size + 10, font=self.font, weight=self.font_weight)
            else:
                axis.set_xlabel("Fragment break index $\\mathbf{i}$",
                                size=self.font_size + 10, font=self.font, weight=self.font_weight)

            axis.set_xlim(min(base_indices), max(base_indices))

            x_ticks = [i for i in range(int(min(base_indices)), int(max(base_indices))) if i % 10 == centroid % 10]
            x_tick_labels = [str(int(i)) for i in x_ticks]
            x_tick_labels[x_tick_labels.index((str(int(centroid))))] = str(int(centroid)) + "*"
            plt.xticks(x_ticks, x_tick_labels, weight='bold', fontname=self.font, size=self.font_size-5)

            axis.set_ylim(0, total_max + 0.005)

            if gauss:
                mu, sigma, x_gauss, y_gauss = Fitting.fit_gauss(data[i], True)
                plt.plot(x_gauss, y_gauss, color="red", linewidth=2.5)
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
        else:
            plt.show()
            plt.clf()
        plt.close()

    # Create a scatter plot of Gaussian fit parameter, provide the input files and directory and the name of the
    # parameter as a string ("Mu" or "Sigma"). Optionally provide a path to save the figure to.
    def plot_fit_paramameter_scatter_single(self, file_names_in, directory_in, parameter_name, save_path=None):

        fit_parameters = DataUtilities.read_fit_parameters(file_names_in, directory_in)
        num_of_samples = len(list(fit_parameters.keys()))
        positions = range(1, num_of_samples + 1)

        plt.figure(figsize=(7, 7))

        for sample, position in zip(fit_parameters.keys(), positions):
            sample_id = sample.strip("*")
            if sample_id in self.colors.keys():
                current_color = self.colors[sample_id]
            else:
                current_color = 'black'
            for sample_num in fit_parameters[sample].keys():
                if fit_parameters[sample][sample_num] == -100 or fit_parameters[sample][sample_num] == 0.0:
                    continue
                if sample_num != "Average":
                    plt.scatter(position, fit_parameters[sample][sample_num], alpha=0.5, color=current_color,
                                marker="o", s=self.marker_size-250, edgecolor='black', linewidth=1.5)
                else:
                    plt.scatter(position, fit_parameters[sample][sample_num], alpha=1, color=current_color,
                                marker="*", s=self.marker_size, edgecolor='black', linewidth=1.5)
        if len(list(fit_parameters.keys())) <= 5:
            plt.xticks(positions, [i.strip("*") + "*" for i in fit_parameters.keys()], weight=self.font_weight,
                       fontname=self.font, size=self.font_size)
        else:
            plt.xticks(positions, [i.strip("*") + "*" for i in fit_parameters.keys()], weight=self.font_weight,
                       fontname=self.font, size=self.font_size-10)

        if parameter_name == "Mu":
            plt.yticks([-1.0, -0.5, 0, 0.5, 1.0], ["-1.0", "-0.5", "0", "0.5", "1.0"], weight=self.font_weight,
                       fontname=self.font, size=self.font_size)
            plt.ylabel("$\\mathbf{\\mu}$ difference from center position", weight=self.font_weight, fontname=self.font,
                       size=self.font_size + 10)
        elif parameter_name == "Sigma":
            plt.ylabel("Breaking distribution width $\\mathbf{\\sigma}$", weight=self.font_weight, fontname=self.font,
                       size=self.font_size + 10)

            total_min = min([min(fit_parameters[i].values()) for i in fit_parameters.keys()])
            total_max = max([max(fit_parameters[i].values()) for i in fit_parameters.keys()])
            y_ticks = [i for i in range(int(round(total_min, 0)), int(round(total_max, 0)) + 1)]
            # y_ticks = [4, 4.5, 5, 5.5, 6]
            plt.yticks(y_ticks, [str(i) for i in y_ticks],
                       font=self.font, size=self.font_size, weight=self.font_weight)
        plt.tick_params(width=2)

        plt.xlabel("Sample", weight=self.font_weight, fontname=self.font,
                   size=self.font_size + 10)
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
        else:
            plt.show()
            plt.clf()
        plt.close()

        # Create a scatter plot of Gaussian fit parameter, provide the input files and directory and the name of the
        # parameter as a string ("Mu" or "Sigma"). Optionally provide a path to save the figure to.

    def plot_fit_paramameter_scatter_single_fit(self, file_names_in, directory_in, parameter_name, save_path=None):

        fit_parameters = DataUtilities.read_fit_parameters(file_names_in, directory_in)
        num_of_samples = len(list(fit_parameters.keys()))
        positions = range(0, num_of_samples)

        plt.figure(figsize=(7, 7))

        for sample, position in zip(fit_parameters.keys(), positions):
            sample_id = sample.strip("*")
            if sample_id in self.colors.keys():
                current_color = self.colors[sample_id]
            else:
                current_color = 'black'
            for sample_num in fit_parameters[sample].keys():
                if fit_parameters[sample][sample_num] == -100 or fit_parameters[sample][sample_num] == 0.0:
                    continue
                if sample_num != "Average":
                    plt.scatter(position, fit_parameters[sample][sample_num], alpha=0.5, color=current_color,
                                marker="o", s=self.marker_size-250, edgecolor='black', linewidth=1.5)
                else:
                    plt.scatter(position, fit_parameters[sample][sample_num], alpha=1, color=current_color,
                                marker="*", s=self.marker_size, edgecolor='black', linewidth=1.5)
        if len(list(fit_parameters.keys())) <= 5:
            plt.xticks(positions, [i.strip("*") + "*" for i in fit_parameters.keys()], weight=self.font_weight,
                       fontname=self.font, size=self.font_size)
        else:
            plt.xticks(positions, [i.strip("*") + "*" for i in fit_parameters.keys()], weight=self.font_weight,
                       fontname=self.font, size=self.font_size - 10)

        significance = {('0', '4'): '*', ('0', '3'): '**', ('0','2'): '*', ('1','3'): '*', ('2','3'): '*'} 
        ax = plt.gca()
        y_max = ax.get_ylim()[1] 
        y = y_max + 0.1
        for (pos1, pos2), stars in significance.items(): 
            
            ax.annotate('', xy=(int(pos1), y), xycoords='data', xytext=(int(pos2), y), textcoords='data', arrowprops=dict(arrowstyle='|-|', mutation_scale = 4, lw=2.5, color='black')) 
            ax.text((int(pos1)+int(pos2))/2, y+0.01, stars, horizontalalignment='center', verticalalignment='center', fontsize=self.font_size, fontweight=self.font_weight, fontname=self.font)
            y += 0.15

        x = [0, 1, 2, 3, 4]
        y = [[fit_parameters[sample][a] for a in fit_parameters[sample] if a != "Average"] for sample in fit_parameters.keys()]
        if parameter_name == "Sigma":
            result, prediction, prediction_with_ci, x_values_large = Fitting.fit_linear(x, y)
            m = result.params[1]
            b = result.params[0]
            r2 = result.rsquared
            label=str(round(m, 2)) + "$\\mathbf{\\cdot x +}$" + str(round(b, 2)) + ", with " + "$\\mathbf{R^2} =$" + " " + str(round(r2, 3))
            plt.plot(x_values_large, prediction.predicted_mean, linestyle='--', color="black", linewidth=2.5,
                     label=label)
            plt.plot(x_values_large, prediction_with_ci.obs_ci_lower, linestyle='--', color="grey", linewidth=1.5)
            plt.plot(x_values_large, prediction_with_ci.obs_ci_upper, linestyle='--', color="grey", linewidth=1.5,)
            plt.fill_between(x_values_large, prediction_with_ci.obs_ci_lower, prediction_with_ci.obs_ci_upper, color = "grey", alpha = 0.2, label="Confidence interval of fit")
            plt.xlim(-0.15, 4.15)
        if parameter_name == "Mu":
            plt.yticks([-1.0, -0.5, 0, 0.5, 1.0], ["-1.0", "-0.5", "0", "0.5", "1.0"], weight=self.font_weight,
                       fontname=self.font, size=self.font_size)
            plt.ylabel("$\\mathbf{\\mu}$ difference from center position", weight=self.font_weight, fontname=self.font,
                       size=self.font_size + 10)
        elif parameter_name == "Sigma":
            plt.ylabel("Breaking distribution width $\\mathbf{\\sigma}$", weight=self.font_weight, fontname=self.font,
                       size=self.font_size + 10)
            # total_min = min([min(fit_parameters[i].values()) for i in fit_parameters.keys()])
            # total_max = max([max(fit_parameters[i].values()) for i in fit_parameters.keys()])
            # y_ticks = [i for i in range(int(round(total_min, 0)), int(round(total_max, 0)) + 1)]
            y_ticks = [4, 4.5, 5, 5.5, 6, 6.5]
            plt.yticks(y_ticks, [str(i) for i in y_ticks],
                       font=self.font, size=self.font_size, weight=self.font_weight)
        plt.tick_params(width=2)
        plt.xlabel("Sample", weight=self.font_weight, fontname=self.font,
                   size=self.font_size + 10)
        legend_properties = {'weight': self.font_weight, 'family': self.font, 'size': f"{self.font_size - 5}"}
        plt.legend(prop=legend_properties, loc='lower right')
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
        else:
            plt.show()
            plt.clf()
        plt.close()

    # Create a side-by-side scatter plot of Gaussian fit parameter, provide the input files and directory and the
    # name of the parameter as a string ("Mu" or "Sigma"). Optionally provide a path to save the figure to.
    def plot_fit_paramameter_scatter_double(self, file_names_in, directory_in, parameter_name, save_path=None):
        all_values_dict = DataUtilities.read_fit_parameters(file_names_in, directory_in)
        gaussian_values = {}
        raw_values = {}
        for sample in all_values_dict.keys():
            if sample.find("fit") != -1:
                gaussian_values[sample.split(" ")[0]] = all_values_dict[sample]
            if sample.find("data") != -1:
                raw_values[sample.split(" ")[0]] = all_values_dict[sample]
        print(raw_values)

        fit_parameters_all = [gaussian_values, raw_values]
        fig, axs = plt.subplots(1, 2, figsize=(14, 7))  # Adjust figsize as needed
        for i in range(0, 2):
            fit_parameters = fit_parameters_all[i]
            current_ax = axs[i]
            num_of_samples = len(list(fit_parameters.keys()))
            positions = range(1, num_of_samples + 1)
            namespace = fit_parameters.keys() #["400", "1500", "704", "GC", "AT"]
            for sample, position in zip(namespace, positions):
                sample_id = sample.strip("*")
                if sample_id in self.colors.keys():
                    current_color = self.colors[sample_id]
                else:
                    current_color = 'black'
                for sample_num in fit_parameters[sample].keys():
                    if sample_num != "Average":
                        current_ax.scatter(position, fit_parameters[sample][sample_num], alpha=0.5, color=current_color,
                                           marker="o", s=self.marker_size, edgecolor='black', linewidth=1.5)
                    else:
                        current_ax.scatter(position, fit_parameters[sample][sample_num], alpha=1, color=current_color,
                                           marker="*", s=self.marker_size, edgecolor='black', linewidth=1.5)
            current_ax.set_xticks(positions, [i.strip("*") + "*" for i in namespace],
                                  weight=self.font_weight,
                                  fontname=self.font, size=self.font_size)

            if parameter_name == "Mu":
                current_ax.set_yticks([-1, 0, 1], ["-1", "0", "+1"], weight=self.font_weight, fontname=self.font,
                                      size=self.font_size)
                current_ax.set_ylabel("$\\mathbf{\\mu}$ difference from nick position", weight=self.font_weight,
                                      fontname=self.font,
                                      size=self.font_size + 10)
            elif parameter_name == "Sigma":
                current_ax.set_ylabel("Breaking distribution width $\\mathbf{\\sigma}$", weight=self.font_weight,
                                      fontname=self.font,
                                      size=self.font_size + 10)

                total_min = min([min(fit_parameters[i].values()) for i in fit_parameters.keys()])
                total_max = max([max(fit_parameters[i].values()) for i in fit_parameters.keys()])
                y_ticks = [i for i in range(int(round(total_min, 0)), int(round(total_max, 0)) + 1)]
                current_ax.set_yticks(y_ticks, [str(int(i)) for i in y_ticks],
                                      font=self.font, size=self.font_size, weight=self.font_weight)
                plt.tick_params(width=2)
            current_ax.set_xlabel("Sample", weight=self.font_weight, fontname=self.font, size=self.font_size)
        axs[0].set_title(weight=self.font_weight, fontname=self.font, size=self.font_size, label="A: From fit")
        axs[1].set_title(weight=self.font_weight, fontname=self.font, size=self.font_size, label="B: From data")
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
        else:
            plt.show()
            plt.clf()
        plt.close()

    # Create overlay histogramm of three breaking distributions, make sure only three are provided that are all
    # normalized and shifted. Provide in order from top to bot and provide a set of colors for the plot. Optionally
    # provide a path to save the figure to.
    def plot_overlay_histogramm(self, file_names_in, directory_in, colors, save_path=None):
        data_breaking_all = DataUtilities.read_breaking_results(file_names_in, directory_in)
        plotting_dict = {}

        # Select only top strand and seperate values and standard deviation
        for key in data_breaking_all.keys():
            sample = key.split("#")[0]
            position = key.split("#")[1]
            data_type = key.split('#')[2]
            if position != "TOP":
                continue
            if sample not in plotting_dict.keys():
                plotting_dict[sample] = {}
            plotting_dict[sample][data_type] = data_breaking_all[key]

        fig, ax = plt.subplots(figsize=(7, 9))

        # Plot the horizontal lines first (behind the distributions)

        ax.axhline(y=0.03, color='black', linewidth=1.0, zorder=1)
        ax.axhline(y=0.06, color='black', linewidth=1.0, zorder=1)
        shift = 0.06
        x_min = 1000
        x_max = -1
        for sample, color in zip(plotting_dict.keys(), colors):
            x = list(plotting_dict[sample]["AVERAGE"].keys())
            y = list(plotting_dict[sample]["AVERAGE"].values())
            y_er = list(plotting_dict[sample]["STD"].values())
            ax.bar(x, y, bottom=shift, label=sample, zorder=2, yerr=y_er, capsize=5, ecolor=self.colors["BlackFaded"],
                   color=color)
            if min(x) < x_min:
                x_min = min(x)
            if max(x) > x_max:
                x_max = max(x)
            mu, sigma, x_gauss, y_gauss = Fitting.fit_gauss(plotting_dict[sample]["AVERAGE"], True)
            y_gauss = [i + shift for i in y_gauss]
            ax.plot(x_gauss, y_gauss, color="black", linewidth=2.5)
            shift -= 0.03

            ax.set_ylim(0)
        ax.set_xlim(left=-40, right=40)

        ax.axvline(x=0, color='red', linestyle='--', linewidth=2.5)

        ax.yaxis.set_ticks([])
        ax.spines['top'].set_visible(True)
        ax.spines['right'].set_visible(True)
        ax.spines['bottom'].set_visible(True)
        ax.spines['left'].set_visible(True)
        ax.set_ylabel("Fraction of reads", fontname=self.font, weight=self.font_weight, size=self.font_size + 10)

        x_ticks = [i for i in range(int(x_min), int(x_max)) if i % 20 == 0]
        x_tick_labels = []
        for x_tick in x_ticks:
            if x_tick > 0:
                x_tick_labels.append("+" + str(int(x_tick)))
            else:
                x_tick_labels.append(str(int(x_tick)))
        plt.xticks(x_ticks, x_tick_labels, weight=self.font_weight, fontname=self.font, size=self.font_size)

        # legend_properties = {'weight': self.font_weight, 'family': self.font, 'size': f"{self.font_size - 5}"}
        # plt.legend(prop=legend_properties)

        ax.set_xlabel("Relative break index", fontname=self.font, weight=self.font_weight, size=self.font_size + 10)
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
        else:
            plt.show()
            plt.clf()
        plt.close()

    # Plot the quality values. Provide the file names and directory and optionally a path to save the figure to.
    def plot_quality_values(self, file_names_in, directory_in, save_path=None):
        quality_results = DataUtilities.read_quality_values(file_names_in, directory_in)
        plt.figure(figsize=(9, 9))
        for sample_name in quality_results.keys():
            cur_sample = quality_results[sample_name]
            current_x = sorted(list(cur_sample.keys()))
            current_y = [cur_sample[i] for i in current_x]
            plt.plot(current_x, current_y, label=sample_name, linewidth=2.5)

        plt.xlabel("Distance from 5'-end [bases]",
                   font=self.font, weight=self.font_weight, size=self.font_size + 10)
        plt.ylabel("Phred quality value",
                   font=self.font, weight=self.font_weight, size=self.font_size + 10)
        plt.xticks(weight=self.font_weight, fontname=self.font, size=self.font_size)
        plt.yticks(weight=self.font_weight, fontname=self.font, size=self.font_size)
        plt.tick_params(width=2)
        legend_properties = {'weight': self.font_weight, 'family': self.font, 'size': f"{self.font_size - 5}"}
        plt.legend(prop=legend_properties)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
        else:
            plt.show()
            plt.clf()
        plt.close()

    # Plot a p matrix from a file and directory as a nxn heatmap, only works with binned p matrices.
    def plot_p_value_matrix(self, file_name_in, directory_in, save_path=None):
        matrix, names_fixed = DataUtilities.read_p_matrix(file_name_in, directory_in)
        names_fixed = [name + "*" for name in names_fixed]
        fig, ax = plt.subplots(figsize=(9, 9))
        
        colors = ['white', 'white']
        cmap = LinearSegmentedColormap.from_list("Custom_cmap", colors, N=256)
        cax = ax.matshow(matrix, interpolation='nearest', cmap=cmap, vmin=0, vmax=1)
        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                addition = ""
                if matrix[i][j]<0.05:
                    addition = "*"
                if matrix[i][j]<0.01:
                    addition = "**"
                if matrix[i][j]<0.001:
                    addition = "***"
                if matrix[i][j]<0.0001:
                    addition = "****"
                if addition != "":
                    if matrix[i][j] < 0.001:
                        ax.text(j,i, f"{matrix[i][j]:.1e}{addition}", va='center', ha='center', color='black', fontsize=self.font_size-5, fontweight=self.font_weight)
                    else:
                        ax.text(j,i, f"{matrix[i][j]:.3f}{addition}", va='center', ha='center', color='black', fontsize=self.font_size-5, fontweight=self.font_weight)

                else: 
                    ax.text(j,i, f"{matrix[i][j]:.3f}", va='center', ha='center', color='black', fontsize=self.font_size-5)

        ax.grid(which='minor', color='black',linewidth=2)
        ax.set_xticks([a-0.5 for a in range(len(names_fixed))], minor=True)
        ax.set_yticks([a-0.5 for a in range(len(names_fixed))], minor=True)
        ax.set_xticks(range(len(names_fixed)))
        ax.set_yticks(range(len(names_fixed)))
        ax.tick_params(width=0, which='minor')
        ax.tick_params(width=2, which='major')

        ax.set_xticklabels(names_fixed, rotation=90,
                           fontsize=self.font_size, fontname=self.font, fontweight=self.font_weight)
        ax.set_yticklabels(names_fixed,
                           fontsize=self.font_size, fontname=self.font, fontweight=self.font_weight)
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
        else:
            plt.show()
            plt.clf()
        plt.close()

    def plot_relative_breakpairs(self, data_dict_in, save_path=None):
        a_bases = ["AA", "CA", 'GA', "TA"]
        c_bases = ["AC", "CC", 'GC', "TC"]
        g_bases = ["AG", "CG", 'GG', "TG"]
        t_bases = ["AT", "CT", 'GT', "TT"]
        plt.figure(figsize=(9, 9))
        plt.errorbar([1, 2, 3, 4], [data_dict_in[o][0] for o in a_bases], marker="o", ls="--",
                     yerr=[data_dict_in[o][1] for o in a_bases], color='blue', capsize=3, capthick=1)
        plt.errorbar([6, 7, 8, 9], [data_dict_in[o][0] for o in c_bases], marker="o", ls="--",
                     yerr=[data_dict_in[o][1] for o in c_bases], color='red', capsize=3, capthick=1)
        plt.errorbar([11, 12, 13, 14], [data_dict_in[o][0] for o in g_bases], marker="o", ls="--",
                     yerr=[data_dict_in[o][1] for o in g_bases], color='orange', capsize=3, capthick=1)
        plt.errorbar([16, 17, 18, 19], [data_dict_in[o][0] for o in t_bases], marker="o", ls="--",
                     yerr=[data_dict_in[o][1] for o in t_bases], color='green', capsize=3, capthick=1)
        plt.xticks([1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14, 16, 17, 18, 19], a_bases + c_bases + g_bases + t_bases,
                   weight=self.font_weight, fontname=self.font, size=self.font_size - 5)
        plt.yticks(weight=self.font_weight, fontname=self.font, size=self.font_size - 5)
        plt.xlabel("Break pair in 5'-XY-3' direction", weight=self.font_weight, fontname=self.font, size=self.font_size)
        plt.ylabel("Relative breaking intensity", weight=self.font_weight, fontname=self.font, size=self.font_size)
        plt.tick_params(width=2)

        if save_path:
            plt.savefig(save_path)
        else:
            plt.show()
            plt.clf()
        plt.close()

    def plot_relative_breakpairs_per_sample_scatterplot(self, data_dict_in, name_dict_in, save_path=None):
        markers = ["x", "o", "*", "^", "h", "D", "+", "s"]
        a_bases = ["AA", "CA", 'GA', "TA"]
        c_bases = ["AC", "CC", 'GC', "TC"]
        g_bases = ["AG", "CG", 'GG', "TG"]
        t_bases = ["AT", "CT", 'GT', "TT"]
        # bases_ordered = a_bases + c_bases + g_bases + t_bases

        plt.figure(figsize=(9, 9))
        for sample_key, marker in zip(list(data_dict_in.keys()), markers):
            sample_dict = data_dict_in[sample_key]
            plt.scatter([1, 2, 3, 4], [sample_dict[o] for o in a_bases], color='blue', marker=marker)
            plt.scatter([6, 7, 8, 9], [sample_dict[o] for o in c_bases], color='red', marker=marker)
            plt.scatter([11, 12, 13, 14], [sample_dict[o] for o in g_bases], color='orange', marker=marker)
            plt.scatter([16, 17, 18, 19], [sample_dict[o] for o in t_bases], color='green', marker=marker)
            label_clean = name_dict_in[sample_key]
            plt.scatter([], [], marker=marker, label=label_clean, color="black")

        plt.xticks([1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14, 16, 17, 18, 19], a_bases + c_bases + g_bases + t_bases,
                   weight=self.font_weight, fontname=self.font, size=self.font_size - 5)
        plt.yticks(weight=self.font_weight, fontname=self.font, size=self.font_size)
        plt.xlabel("Break pair in 5'-XY-3' direction", weight=self.font_weight, fontname=self.font,
                   size=self.font_size + 10)
        plt.ylabel("Relative breaking intensity", weight=self.font_weight, fontname=self.font, size=self.font_size + 10)
        plt.tick_params(width=2)
        legend_properties = {'weight': self.font_weight, 'family': self.font, 'size': f"{self.font_size - 5}"}
        plt.legend(prop=legend_properties)
        if save_path:
            plt.savefig(save_path)
        else:
            plt.show()
            plt.clf()
        plt.close()

    def plot_relative_breakpairs_per_sample_boxplot(self, data_dict_in, save_path=None):
        a_bases = ["AA", "CA", 'GA', "TA"]
        c_bases = ["AC", "CC", 'GC', "TC"]
        g_bases = ["AG", "CG", 'GG', "TG"]
        t_bases = ["AT", "CT", 'GT', "TT"]
        bases_ordered = a_bases + c_bases + g_bases + t_bases

        plt.figure(figsize=(9, 9))
        all_vals = []
        for base in bases_ordered:
            values_current_base = []
            for sample_key in data_dict_in.keys():
                sample_dict = data_dict_in[sample_key]
                values_current_base.append(sample_dict[base])
            all_vals.append(values_current_base)

        bplot = plt.boxplot(all_vals, positions=[1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14, 16, 17, 18, 19],
                            patch_artist=True)
        for patch, color in zip(bplot['boxes'], ["blue", "blue", "blue", "blue", "red", "red", "red", "red",
                                                 "orange", "orange", "orange", "orange", "green", "green", "green",
                                                 "green"]):
            patch.set_facecolor(color)
        plt.xticks([1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14, 16, 17, 18, 19], a_bases + c_bases + g_bases + t_bases,
                   weight=self.font_weight, fontname=self.font, size=self.font_size - 5)
        plt.yticks(weight=self.font_weight, fontname=self.font, size=self.font_size)
        plt.xlabel("Break pair in 5'-XY-3' direction", weight=self.font_weight, fontname=self.font,
                   size=self.font_size + 10)
        plt.ylabel("Relative breaking intensity", weight=self.font_weight, fontname=self.font, size=self.font_size + 10)
        plt.tick_params(width=2)
        if save_path:
            plt.savefig(save_path)
        else:
            plt.show()
            plt.clf()
        plt.close()


if __name__ == "__main__":
    # Example on how to use this script to create figures from data created in FilteringAlgorithm.py
    plotter = PlottingNGS()
    plotter.plot_overlay_histogramm(
        ["GC_Breaking_Dist_Norm_Shifted_Average.csv", "704_Breaking_Dist_Norm_Shifted_Average.csv",
         "AT_Breaking_Dist_Norm_Shifted_Average.csv"], "../ProcessedData/704_AT_GC_400_1500_ATMi/",
        [plotter.colors["GC"], plotter.colors["704"], plotter.colors["AT"]])