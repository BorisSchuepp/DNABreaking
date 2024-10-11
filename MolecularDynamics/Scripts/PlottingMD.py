# Script to plot the processed data from MD simulation analysis
# Author: Boris N. Schüpp

import matplotlib.pyplot as plt
from matplotlib import gridspec


class PlottingMD:

    def __init__(self):
        self.font_size = 20
        self.font = "Arial"
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
                       "Teal": (0, 0.5, 0.5, 1),
                       "P-O5'": "red",
                       "O5'-C5'": "orange",
                       "C5'-C4'": "blue",
                       "C4'-C3'": "green",
                       "C3'-O3'": "yellow",
                       "O3'-P": (0, 0.5, 0.5, 1)
                       }
        self.marker_size = 400
        plt.rcParams['mathtext.fontset'] = 'custom'
        plt.rcParams['mathtext.it'] = f'{self.font}:italic'
        plt.rcParams['mathtext.bf'] = f'{self.font}:italic:bold'
        plt.rcParams['axes.linewidth'] = 2

    def plot_base_pair_distances_all_force_variation(self, data_dict_in, save_path=None):
        sequences = set([key[0] for key in data_dict_in.keys() if key[2] != -1])
        sequences = sorted(sequences)
        sequences_final = [i for i in sequences if len(i) == 3]
        sequences_final += [i for i in sequences if len(i) == 2]
        fig, axs = plt.subplots(2, 3, figsize=(18, 12))
        x_values = [i for i in range(-49, 51)]
        for sequence_id, sequence in enumerate(sequences_final):
            row = int(sequence_id / 3)
            col = int(sequence_id % 3)
            for key in data_dict_in.keys():
                if key[0] != sequence:
                    continue

                if key[2] == -1:
                    alpha = 1
                    width = 2.5
                else:
                    alpha = 0.2
                    width = 1

                if key[1] == 0:
                    color = self.colors["Teal"]
                else:
                    color = self.colors["704"]

                axs[row, col].plot(x_values, data_dict_in[key], color=color, alpha=alpha)
                axs[row, col].set_title(sequence[0] + "." + sequence[1] + " nN",
                                        fontname=self.font, weight=self.font_weight, size=self.font_size + 10)
        x_tick = 0
        x_ticks = [-40, -20, 0, 20, 40]
        x_tick_labels = ["-40", "-20", "0", "+20", "+40"]

        for ax_ix, ax in enumerate(axs.flatten()):
            ax.set_ylabel(f"Average basepair distance [nm]",
                          fontname=self.font, weight=self.font_weight, size=self.font_size + 5)
            ax.set_xlabel(f"Relative base index",
                          fontname=self.font, weight=self.font_weight, size=self.font_size + 5)
            ax.set_xticks(x_ticks)
            y_tick_labels = ax.get_yticklabels()
            ax.set_xticklabels(x_tick_labels,
                               fontname=self.font, weight=self.font_weight, size=self.font_size)
            for label in ax.get_yticklabels():
                label.set_fontproperties(self.font)
                label.set_fontsize(self.font_size)
                label.set_fontweight(self.font_weight)
            ax.tick_params(width=2)
            ax.axhline(y=1, color="black", linestyle='--', linewidth=2.5)
            ax.set_xlim(-41, 41)
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
            plt.clf()
        else:
            plt.show()
        plt.close()

    def plot_force_bond_type_all_force_variation(self, forces_dict_in, save_path=None):
        sequences_final = ["05", "10", "15", "20", "25", "30"]
        fig, axs = plt.subplots(2, 3, figsize=(18, 12))
        for sequence_id, sequence in enumerate(sequences_final):
            row = int(sequence_id / 3)
            col = int(sequence_id % 3)
            for key in forces_dict_in.keys():
                if key[0] != sequence or key[2] != -1 or key[1] != 1:
                    continue
                axs[row, col].plot([i - 50 for i in forces_dict_in[key][0]], forces_dict_in[key][1],
                                   label=key[3], linewidth=2.5, color=self.colors[key[3]])
                title = sequence[0] + "." + sequence[1] + " nN"
                axs[row, col].set_title(title,
                                        fontname=self.font, weight=self.font_weight, size=self.font_size + 10)
        x_tick = 0
        x_ticks = [-40, -20, 0, 20, 40]
        x_tick_labels = ["-40", "-20", "0", "+20", "+40"]

        for ax_ix, ax in enumerate(axs.flatten()):
            ax.set_ylabel(f"Average force [nN]",
                          fontname=self.font, weight=self.font_weight, size=self.font_size + 5)
            ax.set_xlabel(f"Relative base index",
                          fontname=self.font, weight=self.font_weight, size=self.font_size + 5)
            ax.set_xticks(x_ticks)
            y_tick_labels = ax.get_yticklabels()
            ax.set_xticklabels(x_tick_labels,
                               fontname=self.font, weight=self.font_weight, size=self.font_size)
            for label in ax.get_yticklabels():
                label.set_fontproperties(self.font)
                label.set_fontsize(self.font_size)
                label.set_fontweight(self.font_weight)
            ax.tick_params(width=2)
            ax.axvline(x=0, color='black', linestyle='--', linewidth=2.5)
            ax.set_xlim(-41, 41)
        plt.tight_layout()
        plt.legend()
        if save_path:
            plt.savefig(save_path)
            plt.clf()
        else:
            plt.show()
        plt.close()

    def plot_end_to_end_distance_all_force_variation(self, data_dict_in, x_data_in, save_path=None):
        sequences = set([key[0] for key in data_dict_in.keys() if key[2] != -1])
        sequences_final = sorted(sequences)

        fig, axs = plt.subplots(2, 3, figsize=(18, 12))

        for sequence_id, sequence in enumerate(sequences_final):
            row = int(sequence_id / 3)
            col = int(sequence_id % 3)
            for key in data_dict_in.keys():
                if key[0] != sequence:
                    continue

                if key[2] == -1:
                    alpha = 1
                    width = 2.5
                else:
                    alpha = 0.2
                    width = 1
                if key[1] == 0:
                    color = self.colors["Teal"]
                else:
                    color = self.colors["704"]
                axs[row, col].plot(x_data_in[0:len(data_dict_in[key])], data_dict_in[key],
                                   color=color, alpha=alpha, linewidth=width)
                axs[row, col].set_title(sequence[0] + "." + sequence[1] + " nN",
                                        fontname=self.font, weight=self.font_weight, size=self.font_size + 10)
        x_tick = 0
        x_ticks = []
        max_x = x_data_in[-1]

        while x_tick * 20 <= max_x:
            x_ticks.append(x_tick * 20)
            x_tick += 1
        y_ticks = [round(1 + 0.2 * i, 2) for i in range(0, 8)]
        for ax in axs.flatten():
            ax.set_ylabel(f"Fractional extension",
                          fontname=self.font, weight=self.font_weight, size=self.font_size + 10)
            ax.set_xlabel(f"Simulation time [ns]",
                          fontname=self.font, weight=self.font_weight, size=self.font_size + 10)
            ax.set_xticks(x_ticks)
            ax.set_yticks(y_ticks)
            ax.set_xticklabels([str(i) for i in x_ticks],
                               fontname=self.font, weight=self.font_weight, size=self.font_size)
            ax.set_yticklabels([str(i) for i in y_ticks],
                               fontname=self.font, weight=self.font_weight, size=self.font_size)
            ax.tick_params(width=2)
            ax.set_ylim(1, 2.5)
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
            plt.clf()
        else:
            plt.show()
        plt.close()

    def plot_end_to_end_distance_all(self, data_dict_in, x_data_in, save_path=None):
        sequences = set([key[0] for key in data_dict_in.keys() if key[2] != -1])
        sequences = sorted(sequences)
        sequences_final = [i for i in sequences if len(i) == 3]
        sequences_final += [i for i in sequences if len(i) == 2]
        fig, axs = plt.subplots(3, 3, figsize=(18, 18))

        for sequence_id, sequence in enumerate(sequences_final):
            row = int(sequence_id / 3)
            col = int(sequence_id % 3)
            for key in data_dict_in.keys():
                if key[0] != sequence:
                    continue

                if key[2] == -1:
                    alpha = 1
                    width = 2.5
                else:
                    alpha = 0.2
                    width = 1
                if key[1] == 0:
                    color = self.colors["Teal"]
                else:
                    color = self.colors["704"]

                axs[row, col].plot(x_data_in[0:len(data_dict_in[key])], data_dict_in[key],
                                   color=color, alpha=alpha, linewidth=width)
                axs[row, col].set_title(sequence + "*",
                                        fontname=self.font, weight=self.font_weight, size=self.font_size + 10)
        x_tick = 0
        x_ticks = []
        max_x = x_data_in[-1]

        while x_tick * 20 < max_x:
            x_ticks.append(x_tick * 20)
            x_tick += 1
        y_ticks = [1 + 0.2 * i for i in range(0, 7)]
        for ax in axs.flatten():
            ax.set_ylabel(f"Fractional extension",
                          fontname=self.font, weight=self.font_weight, size=self.font_size + 10)
            ax.set_xlabel(f"Simulation time [ns]",
                          fontname=self.font, weight=self.font_weight, size=self.font_size + 10)
            ax.set_xticks(x_ticks)
            ax.set_yticks(y_ticks)
            ax.set_xticklabels([str(i) for i in x_ticks],
                               fontname=self.font, weight=self.font_weight, size=self.font_size)
            ax.set_yticklabels([str(i) for i in y_ticks],
                               fontname=self.font, weight=self.font_weight, size=self.font_size)
            ax.tick_params(width=2)
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
            plt.clf()
        else:
            plt.show()
        plt.close()

    def plot_base_pair_distances_all(self, data_dict_in, save_path=None):
        sequences = set([key[0] for key in data_dict_in.keys() if key[2] != -1])
        sequences = sorted(sequences)
        sequences_final = [i for i in sequences if len(i) == 3]
        sequences_final += [i for i in sequences if len(i) == 2]
        fig, axs = plt.subplots(3, 3, figsize=(18, 18))
        x_values = [i for i in range(-49, 51)]
        for sequence_id, sequence in enumerate(sequences_final):
            row = int(sequence_id / 3)
            col = int(sequence_id % 3)
            for key in data_dict_in.keys():
                if key[0] != sequence:
                    continue

                if key[2] == -1:
                    alpha = 1
                    width = 2.5
                else:
                    alpha = 0.2
                    width = 1

                if key[1] == 0:
                    color = self.colors["Teal"]
                else:
                    color = self.colors["704"]

                axs[row, col].plot(x_values, data_dict_in[key], color=color, alpha=alpha)
                axs[row, col].set_title(sequence + "*",
                                        fontname=self.font, weight=self.font_weight, size=self.font_size + 10)
        x_tick = 0
        x_ticks = [-40, -20, 0, 20, 40]
        x_tick_labels = ["-40", "-20", "0", "+20", "+40"]

        for ax_ix, ax in enumerate(axs.flatten()):
            ax.set_ylabel(f"Average basepair distance [nm]",
                          fontname=self.font, weight=self.font_weight, size=self.font_size + 5)
            ax.set_xlabel(f"Relative base index",
                          fontname=self.font, weight=self.font_weight, size=self.font_size + 5)
            ax.set_xticks(x_ticks)
            y_tick_labels = ax.get_yticklabels()
            ax.set_xticklabels(x_tick_labels,
                               fontname=self.font, weight=self.font_weight, size=self.font_size)
            for label in ax.get_yticklabels():
                label.set_fontproperties(self.font)
                label.set_fontsize(self.font_size)
                label.set_fontweight(self.font_weight)
            ax.tick_params(width=2)
            ax.axvline(x=0, color='black', linestyle='--', linewidth=2.5)
            ax.axhline(y=1, color="black", linestyle='--', linewidth=2.5)
            ax.set_xlim(-41, 41)
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
            plt.clf()
        else:
            plt.show()
        plt.close()

    def plot_end_to_end_distance_nick_no_nick(self, data_dict_in, x_data_in, save_path=None):
        fig = plt.figure(figsize=(15, 12))

        gs = gridspec.GridSpec(2, 1, height_ratios=[3.7, 1])
        ax = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])

        sequences = ["AAT", "ACT", "AGT", "GAT", "GGT", "GCT", "GTT"]
        for key in data_dict_in.keys():
            if key[0] not in sequences or key[1] not in [0, 1, 2] or key[2] not in [0, 1, 2]:
                continue

            if key[1] == 0:
                color = self.colors["Teal"]
            else:
                color = self.colors["704"]
            ax.plot(x_data_in, data_dict_in[key], color=color, alpha=0.1)
            ax2.plot(x_data_in, data_dict_in[key], color=color, alpha=0.1)

        ax2.set_ylim(0.99, 1.15)
        ax.set_ylim(1.75, 2.3)
        ax.spines['bottom'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax.plot(x_data_in, data_dict_in[("704", 3, -1)], color=self.colors["704"])
        ax2.plot(x_data_in, data_dict_in[("704", 3, -1)], color=self.colors["704"])
        ax.plot(x_data_in, data_dict_in[("704", 0, -1)], color=self.colors["Teal"])
        ax2.plot(x_data_in, data_dict_in[("704", 0, -1)], color=self.colors["Teal"])
        ax2.set_yticks([1.0, 1.1], ["1.0", "1.1"],
                       weight=self.font_weight, fontname=self.font, size=self.font_size)
        ax.set_yticks([1.8, 1.9, 2.0, 2.1, 2.2, 2.3], ["1.8", "1.9", "2.0", "2.1", "2.2", "2.3"],
                      weight=self.font_weight, fontname=self.font, size=self.font_size)

        ax.tick_params(labeltop=False)
        ax2.xaxis.tick_bottom()
        ax.xaxis.tick_top()
        ax.set_xticks([])
        ax2.set_xticks([0, 20, 40, 60, 80, 100], ["0", "20", "40", "60", "80", "100"],
                       weight=self.font_weight, fontname=self.font, size=self.font_size)

        d = .015
        kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
        ax.plot((-d, +d), (-d, +d), **kwargs, linewidth=2)  # top-left diagonal
        ax.plot((1 - d, 1 + d), (-d, +d), **kwargs, linewidth=2)  # top-right diagonal

        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1 - d * 4, 1 + d * 4), **kwargs, linewidth=2)
        ax2.plot((1 - d, 1 + d), (1 - d * 4, 1 + d * 4), **kwargs, linewidth=2)
        fig.subplots_adjust(hspace=0.05)
        label_x_position = 0.06
        label_y_position = 0.5

        fig.text(label_x_position, label_y_position, f"Fractional extension", va='center', rotation='vertical',
                 fontname=self.font, fontweight=self.font_weight, fontsize=self.font_size + 10)
        ax2.set_xlabel("Simulation time [ns]", fontname=self.font, weight=self.font_weight, size=self.font_size + 10)
        ax.tick_params(width=2)
        ax2.tick_params(width=2)
        if not save_path:
            plt.show()
        else:
            plt.savefig(save_path)
            plt.clf()
        plt.close()

    def plot_base_pair_distance_nick_no_nick(self, data_dict_in, save_path=None):
        fig = plt.figure(figsize=(9, 9))
        x_values = [i for i in range(-49, 51)]
        sequences = ["AAT", "ACT", "AGT", "GGT", "GCT", "GTT", "GAT"]
        for key in data_dict_in.keys():
            if key[0] not in sequences or key[1] not in [0, 1, 2] or key[2] not in [0, 1, 2]:
                continue
            if key[1] == 0:
                color = self.colors["Teal"]
            else:
                color = self.colors["704"]
            plt.plot(x_values, data_dict_in[key], color=color, alpha=0.1)

        plt.plot(x_values, data_dict_in[("704", 3, -1)], color=self.colors["704"], linewidth=2.5)
        plt.plot(x_values, data_dict_in[("704", 0, -1)], color=self.colors["Teal"], linewidth=2.5)

        plt.axhline(y=1, color="black", linestyle='--', linewidth=3)
        plt.axvline(x=0, color='black', linestyle='--', linewidth=3)

        plt.tick_params(width=3)
        plt.xlim(-48, 48)

        plt.ylabel(f"Average basepair distance [nm] ",
                   fontname=self.font, weight=self.font_weight, size=self.font_size + 10)
        plt.xlabel("Relative base index",
                   fontname=self.font, weight=self.font_weight, size=self.font_size + 10)

        plt.yticks([0.5, 1, 1.5, 2, 2.5, 3, 3.5], fontname=self.font, weight=self.font_weight, size=self.font_size)
        plt.xticks([-40, -20, 0, 20, 40], labels=["-40", "-20", "0", "+20", "+40"],
                   fontname=self.font, weight=self.font_weight, size=self.font_size)
        plt.ylim(0, 3.7)
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
        else:
            plt.show()
            plt.clf()
        plt.close()

    def plot_base_pair_distance_averages(self, data_dict_in, save_path=None):
        fig = plt.figure(figsize=(7, 7))
        x_values = [i for i in range(-49, 51)]
        sequences = ["AT", "GC", "704"]
        for key in data_dict_in.keys():
            if key[0] not in sequences or key[1] != 3 or key[2] != -1:
                continue
            color = self.colors[key[0]]
            plt.plot(x_values, data_dict_in[key], color=color, linewidth=2.5)

        plt.axhline(y=1, color="black", linestyle='--', linewidth=3)
        plt.axvline(x=0, color='black', linestyle='--', linewidth=3)

        plt.tick_params(width=3)
        plt.xlim(-48, 48)

        plt.ylabel(f"Average basepair distance [nm] ",
                   fontname=self.font, weight=self.font_weight, size=self.font_size + 10)
        plt.xlabel("Relative base index",
                   fontname=self.font, weight=self.font_weight, size=self.font_size + 10)

        plt.yticks(fontname=self.font, weight=self.font_weight, size=self.font_size)
        plt.xticks([-40, -20, 0, 20, 40], labels=["-40", "-20", "0", "+20", "+40"],
                   fontname=self.font, weight=self.font_weight, size=self.font_size)
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
        else:
            plt.show()
            plt.clf()
        plt.close()

    def plot_fraying_force_variation(self, in_1, in_2, in_3, in_4, in_5, in_6, description, save_path=None):

        plt.figure(figsize=(7, 7))
        plt.scatter([0.5 for _ in range(0, len(in_1))], in_1, color=self.colors["P-O5'"], alpha=0.3, marker="o",
                    s=self.marker_size, edgecolor='black', linewidth=1.5)
        plt.scatter(0.5, sum(in_1) / len(in_1), color=self.colors["P-O5'"], alpha=1, marker="*", s=self.marker_size,
                    edgecolor='black', linewidth=1.5)

        plt.scatter([1 for _ in range(0, len(in_2))], in_2, color="orange", alpha=0.3, marker="o",
                    s=self.marker_size, edgecolor='black', linewidth=1.5)
        plt.scatter(1, sum(in_2) / len(in_2), color="orange", alpha=1, marker="*", s=self.marker_size,
                    edgecolor='black', linewidth=1.5)

        plt.scatter([1.5 for _ in range(0, len(in_3))], in_3, color="blue", alpha=0.3, marker="o",
                    s=self.marker_size, edgecolor='black', linewidth=1.5)
        plt.scatter(1.5, sum(in_3) / len(in_3), color="blue", alpha=1, marker="*", s=self.marker_size,
                    edgecolor='black', linewidth=1.5)

        plt.scatter([2 for _ in range(0, len(in_4))], in_4, color="yellow", alpha=0.3, marker="o",
                    s=self.marker_size, edgecolor='black', linewidth=1.5)
        plt.scatter(2, sum(in_4) / len(in_4), color="yellow", alpha=1, marker="*", s=self.marker_size,
                    edgecolor='black', linewidth=1.5)

        plt.scatter([2.5 for _ in range(0, len(in_5))], in_5, color=self.colors["O3'-P"], alpha=0.3, marker="o",
                    s=self.marker_size, edgecolor='black', linewidth=1.5)
        plt.scatter(2.5, sum(in_5) / len(in_5), color=self.colors["O3'-P"], alpha=1, marker="*", s=self.marker_size,
                    edgecolor='black', linewidth=1.5)

        plt.scatter([3 for _ in range(0, len(in_6))], in_6, color=self.colors["GC"], alpha=0.3, marker="o",
                    s=self.marker_size, edgecolor='black', linewidth=1.5)
        plt.scatter(3, sum(in_6) / len(in_6), color=self.colors["GC"], alpha=1, marker="*", s=self.marker_size,
                    edgecolor='black', linewidth=1.5)

        plt.xticks([i * 1 for i in [0.5, 1, 1.5, 2, 2.5, 3]], ["0.5", "1.0", "1.5", "2.0", "2.5", "3.0"],
                   weight=self.font_weight, fontname=self.font, size=self.font_size)
        plt.yticks(weight=self.font_weight, fontname=self.font, size=self.font_size)
        plt.xlabel("Force [nN]", fontname=self.font, weight=self.font_weight, size=self.font_size + 10)
        if description == "Left":
            plt.ylabel("Fraying start base index (left)",
                       fontname=self.font, weight=self.font_weight, size=self.font_size + 10)
            plt.ylim(18, 53)
        if description == "Right":
            plt.ylabel("Fraying end base index (right)",
                       fontname=self.font, weight=self.font_weight, size=self.font_size + 10)
            plt.ylim(50, 80)
        if description == "Width":
            plt.ylabel("Fraying width",
                       fontname=self.font, weight=self.font_weight, size=self.font_size + 10)
            plt.yticks([0, 10, 20, 30, 40, 50, 60],
                       weight=self.font_weight, fontname=self.font, size=self.font_size)
            plt.ylim(-1, 62)

        plt.tight_layout()
        plt.tick_params(width=3)
        if not save_path:
            plt.show()
            plt.clf()
        else:
            plt.savefig(save_path)
        plt.close()

    def plot_fraying(self, in_gc, in_704, in_at, description, save_path=None):
        plt.figure(figsize=(7, 7))
        plt.scatter([2 for i in range(0, len(in_704))], in_704, color=self.colors["704"], alpha=0.3, marker="o",
                    s=self.marker_size, edgecolor='black', linewidth=1.5)
        plt.scatter(2, sum(in_704) / len(in_704), color=self.colors["704"], alpha=1, marker="*", s=self.marker_size,
                    edgecolor='black', linewidth=1.5)
        plt.scatter([3 for i in range(0, len(in_at))], in_at, color=self.colors["AT"], alpha=0.3, marker="o",
                    s=self.marker_size, edgecolor='black', linewidth=1.5)
        plt.scatter(3, sum(in_at) / len(in_at), color=self.colors["AT"], alpha=1, marker="*", s=self.marker_size,
                    edgecolor='black', linewidth=1.5)
        plt.scatter([1 for i in range(0, len(in_gc))], in_gc, color=self.colors["GC"], alpha=0.3, marker="o",
                    s=self.marker_size, edgecolor='black', linewidth=1.5)
        plt.scatter(1, sum(in_gc) / len(in_gc), color=self.colors["GC"], alpha=1, marker="*", s=self.marker_size,
                    edgecolor='black', linewidth=1.5)

        plt.xticks([i * 1 for i in [1, 2, 3]], ["GC*", "704*", "AT*"],
                   weight=self.font_weight, fontname=self.font, size=self.font_size)
        plt.yticks(weight=self.font_weight, fontname=self.font, size=self.font_size)
        plt.xlabel("Sample", fontname=self.font, weight=self.font_weight, size=self.font_size + 10)
        if description == "Left":
            plt.ylabel("Fraying start base index (left)",
                       fontname=self.font, weight=self.font_weight, size=self.font_size + 10)
            plt.ylim(39, 53)
        if description == "Right":
            plt.ylabel("Fraying end base index (right)",
                       fontname=self.font, weight=self.font_weight, size=self.font_size + 10)
            plt.ylim(48, 72)
        if description == "Width":
            plt.ylabel("Fraying width",
                       fontname=self.font, weight=self.font_weight, size=self.font_size + 10)
            plt.yticks([2, 6, 10, 14, 18],
                       weight=self.font_weight, fontname=self.font, size=self.font_size)
            plt.ylim(-1, 21)

        plt.tight_layout()
        plt.tick_params(width=3)
        if not save_path:
            plt.show()
            plt.clf()
        else:
            plt.savefig(save_path)
        plt.close()

    def plot_forces_nick_no_nick(self, forces_dict_in, save_path=None):
        sequences = ["AAT", "ACT", "AGT", "GGT", "GCT", "GTT", "GAT"]
        plt.figure(figsize=(9, 9))
        for key in forces_dict_in.keys():
            if key[0] not in sequences or key[1] not in [0, 1, 2] or key[2] not in [0, 1, 2] or key[3] != "C3'-O3'":
                continue
            if key[1] == 0:
                color = self.colors["Teal"]
            else:
                color = self.colors["704"]
            plt.plot([i - 50 for i in forces_dict_in[key][0]], forces_dict_in[key][1], color=color, alpha=0.1)
        plt.plot([i - 50 for i in forces_dict_in[("704", 0, -1, "C3'-O3'")][0]],
                 forces_dict_in[("704", 0, -1, "C3'-O3'")][1], color=self.colors["Teal"], linewidth=3)
        plt.plot([i - 50 for i in forces_dict_in[("704", 3, -1, "C3'-O3'")][0]],
                 forces_dict_in[("704", 3, -1, "C3'-O3'")][1], color=self.colors["704"], linewidth=3)

        plt.xlabel("Relative base index",
                   fontsize=self.font_size + 10, weight=self.font_weight, font=self.font)
        plt.ylabel("Average force [nN]",
                   fontsize=self.font_size + 10, weight=self.font_weight, font=self.font)

        plt.xlim(-48, 49)
        plt.ylim(0.59, 1.41)

        plt.yticks([0.6, 0.8, 1.0, 1.2, 1.4], fontname=self.font, weight=self.font_weight, size=self.font_size)
        plt.xticks([-40, -20, 0, 20, 40], labels=["-40", "-20", "0", "+20", "+40"],
                   fontname=self.font, weight=self.font_weight, size=self.font_size)
        plt.axvline(x=0, color='black', linestyle='--', linewidth=2.5)
        plt.tick_params(width=2)
        plt.tight_layout()
        if not save_path:
            plt.show()
            plt.clf()
        else:
            plt.savefig(save_path)
        plt.close()

    def plot_force_nick_type_all(self, forces_dict_in, save_path=None):
        sequences = set([key[0] for key in forces_dict_in.keys() if key[2] != -1 and key[0]
                         not in ["GC", "AT", "GTT"]])
        sequences = sorted(sequences)
        sequences_final = [i for i in sequences if len(i) == 3]
        sequences_final += [i for i in sequences if len(i) == 2]
        fig, axs = plt.subplots(2, 3, figsize=(18, 12))

        for sequence_id, sequence in enumerate(sequences_final):
            row = int(sequence_id / 3)
            col = int(sequence_id % 3)
            for key in forces_dict_in.keys():
                if key[0] != sequence or key[1] in [0, 3] or key[3] != "C3'-O3'":
                    continue

                if key[2] == -1:
                    alpha = 1
                    width = 2.5
                else:
                    alpha = 0.2
                    width = 1

                if key[1] == 1:
                    color = "Purple"
                else:
                    color = "Orange"
                axs[row, col].plot([i - 50 for i in forces_dict_in[key][0]], forces_dict_in[key][1],
                                   color=color, alpha=alpha, linewidth=width)
                axs[row, col].set_title(sequence + "*",
                                        fontname=self.font, weight=self.font_weight, size=self.font_size + 10)
        x_tick = 0
        x_ticks = [-40, -20, 0, 20, 40]
        x_tick_labels = ["-40", "-20", "0", "+20", "+40"]

        for ax_ix, ax in enumerate(axs.flatten()):
            ax.set_ylabel(f"Average force [nN]",
                          fontname=self.font, weight=self.font_weight, size=self.font_size + 5)
            ax.set_xlabel(f"Relative base index",
                          fontname=self.font, weight=self.font_weight, size=self.font_size + 5)
            ax.set_xticks(x_ticks)
            y_tick_labels = ax.get_yticklabels()
            ax.set_xticklabels(x_tick_labels,
                               fontname=self.font, weight=self.font_weight, size=self.font_size)
            for label in ax.get_yticklabels():
                label.set_fontproperties(self.font)
                label.set_fontsize(self.font_size)
                label.set_fontweight(self.font_weight)
            ax.tick_params(width=2)
            ax.axvline(x=0, color='black', linestyle='--', linewidth=2.5)
            ax.set_xlim(-41, 41)
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
            plt.clf()
        else:
            plt.show()
        plt.close()

    def plot_force_nick_type_single(self, forces_dict_in, sequence_in, save_path=None):

        plt.figure(figsize=(9, 9))
        x_values = [i for i in range(-49, 51)]

        for key in forces_dict_in.keys():
            if key[0] != sequence_in or key[1] in [0, 3] or key[3] != "C3'-O3'":
                continue

            if key[2] == -1:
                alpha = 1
                width = 2.5
            else:
                alpha = 0.2
                width = 1

            if key[1] == 1:
                color = "Purple"
            else:
                color = "Orange"
            plt.plot([i - 50 for i in forces_dict_in[key][0]], forces_dict_in[key][1],
                     color=color, alpha=alpha, linewidth=width)

        x_tick = 0
        x_ticks = [-40, -20, 0, 20, 40]
        x_tick_labels = ["-40", "-20", "0", "+20", "+40"]

        plt.ylabel(f"Average force [nN]",
                   fontname=self.font, weight=self.font_weight, size=self.font_size + 5)
        plt.xlabel(f"Relative base index",
                   fontname=self.font, weight=self.font_weight, size=self.font_size + 5)
        plt.xticks(x_ticks)
        ax = plt.gca()
        y_tick_labels = ax.get_yticklabels()
        ax.set_xticklabels(x_tick_labels,
                           fontname=self.font, weight=self.font_weight, size=self.font_size)
        for label in ax.get_yticklabels():
            label.set_fontproperties(self.font)
            label.set_fontsize(self.font_size)
            label.set_fontweight(self.font_weight)
        ax.tick_params(width=2)
        ax.axvline(x=0, color='black', linestyle='--', linewidth=2.5)
        ax.set_xlim(-41, 41)
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
            plt.clf()
        else:
            plt.show()
        plt.close()

    def plot_force_bond_type_all(self, forces_dict_in, save_path=None):
        sequences = set([key[0] for key in forces_dict_in.keys() if key[2] != -1 and key[0]
                         not in ["704"]])
        sequences = sorted(sequences)
        sequences_final = [i for i in sequences if len(i) == 3]
        sequences_final += [i for i in sequences if len(i) == 2]
        fig, axs = plt.subplots(3, 3, figsize=(18, 18))
        for sequence_id, sequence in enumerate(sequences_final):
            row = int(sequence_id / 3)
            col = int(sequence_id % 3)
            for key in forces_dict_in.keys():
                if key[0] != sequence or key[2] != -1 or key[1] != 3:
                    continue
                axs[row, col].plot([i - 50 for i in forces_dict_in[key][0]], forces_dict_in[key][1],
                                   label=key[3], linewidth=2.5, color=self.colors[key[3]])
                axs[row, col].set_title(sequence + "*",
                                        fontname=self.font, weight=self.font_weight, size=self.font_size + 10)
        x_tick = 0
        x_ticks = [-40, -20, 0, 20, 40]
        x_tick_labels = ["-40", "-20", "0", "+20", "+40"]

        for ax_ix, ax in enumerate(axs.flatten()):
            ax.set_ylabel(f"Average force [nN]",
                          fontname=self.font, weight=self.font_weight, size=self.font_size + 5)
            ax.set_xlabel(f"Relative base index",
                          fontname=self.font, weight=self.font_weight, size=self.font_size + 5)
            ax.set_xticks(x_ticks)
            y_tick_labels = ax.get_yticklabels()
            ax.set_xticklabels(x_tick_labels,
                               fontname=self.font, weight=self.font_weight, size=self.font_size)
            for label in ax.get_yticklabels():
                label.set_fontproperties(self.font)
                label.set_fontsize(self.font_size)
                label.set_fontweight(self.font_weight)
            ax.tick_params(width=2)
            ax.axvline(x=0, color='black', linestyle='--', linewidth=2.5)
            ax.set_xlim(-41, 41)
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
            plt.clf()
        else:
            plt.show()
        plt.close()

    def plot_force_bond_type_single(self, forces_dict_in, sequence_in, save_path=None):

        plt.figure(figsize=(9, 9))
        x_values = [i for i in range(-49, 51)]

        for key in forces_dict_in.keys():
            if key[0] != sequence_in or key[2] != -1 or key[1] != 3:
                continue

            plt.plot([i - 50 for i in forces_dict_in[key][0]], forces_dict_in[key][1],
                     color=self.colors[key[3]], linewidth=2.5, label=key[3])

        legend_properties = {'weight': self.font_weight, 'family': self.font, 'size': f"{self.font_size - 5}"}
        plt.legend(prop=legend_properties)
        x_tick = 0
        x_ticks = [-40, -20, 0, 20, 40]
        x_tick_labels = ["-40", "-20", "0", "+20", "+40"]

        plt.ylabel(f"Average force [nN]",
                   fontname=self.font, weight=self.font_weight, size=self.font_size + 5)
        plt.xlabel(f"Relative base index",
                   fontname=self.font, weight=self.font_weight, size=self.font_size + 5)
        plt.xticks(x_ticks)
        ax = plt.gca()
        y_tick_labels = ax.get_yticklabels()
        ax.set_xticklabels(x_tick_labels,
                           fontname=self.font, weight=self.font_weight, size=self.font_size)
        for label in ax.get_yticklabels():
            label.set_fontproperties(self.font)
            label.set_fontsize(self.font_size)
            label.set_fontweight(self.font_weight)
        ax.tick_params(width=2)
        ax.axvline(x=0, color='black', linestyle='--', linewidth=2.5)
        ax.set_xlim(-41, 41)
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
            plt.clf()
        else:
            plt.show()
        plt.close()
