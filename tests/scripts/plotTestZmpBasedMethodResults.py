#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt


class PlotTestZmpBasedMethodResults(object):
    def __init__(self, plot_capture_point=False):
        self.fig = plt.figure()
        plt.rcParams["font.size"] = 10

        self.result_data_list = {}

        # Loop for methods
        self.method_list = ["PreviewControlZmp", "DdpZmp", "DcmTracking", "FootGuidedControl", "LinearMpcZmp"]
        for method_idx, method_str in enumerate(self.method_list):
            # Load file
            result_file_path = "/tmp/Test{}.txt".format(method_str)
            self.result_data_list[method_str] = np.genfromtxt(result_file_path, dtype=None, delimiter=None, names=True)
            result_data = self.result_data_list[method_str]
            print("[PlotTestZmpBasedMethodResults] Load {}".format(result_file_path))

            # Loop for axes
            self.axis_list = ["x", "y"]
            for axis_idx, axis_str in enumerate(self.axis_list):
                # Setup subplot
                ax = self.fig.add_subplot(len(self.method_list) + 1, len(self.axis_list), 2 * method_idx + axis_idx + 1)

                # Plot
                planned_zmp_key = "planned_zmp_{}".format(axis_str)
                ax.plot(result_data["time"], result_data[planned_zmp_key],
                        color="red", label="planned ZMP")
                ref_zmp_key = "ref_zmp_{}".format(axis_str)
                if ref_zmp_key in result_data.dtype.names:
                    ax.plot(result_data["time"], result_data[ref_zmp_key],
                            color="blue", linestyle="dashed", label="ref ZMP")
                ref_zmp_min_key = "ref_zmp_min_{}".format(axis_str)
                ref_zmp_max_key = "ref_zmp_max_{}".format(axis_str)
                if ref_zmp_min_key in result_data.dtype.names:
                    ax.plot(result_data["time"], result_data[ref_zmp_min_key],
                            color="blue", linestyle="dashed", label="min ZMP")
                    ax.plot(result_data["time"], result_data[ref_zmp_max_key],
                            color="blue", linestyle="dashed", label="max ZMP")
                com_pos_key = "com_pos_{}".format(axis_str)
                ax.plot(result_data["time"], result_data[com_pos_key],
                        color="green", label="CoM")
                capture_point_key = "capture_point_{}".format(axis_str)
                if plot_capture_point and capture_point_key in result_data.dtype.names:
                    ax.plot(result_data["time"], result_data[capture_point_key],
                            color="orange", label="capture point")

                # Set labels, etc.
                ax.set_title("{}-{}".format(method_str, axis_str.upper()))
                ax.set_xlabel("time [s]")
                ax.set_ylabel("pos [m]".format(axis_str))
                ax.grid()
                if method_idx == 0 and axis_idx == 0:
                    ax.legend(loc="lower right")

        # Plot computation time
        ax = self.fig.add_subplot(len(self.method_list) + 1, len(self.axis_list), len(self.method_list) * len(self.axis_list) + 1)
        ax.bar(np.arange(len(self.method_list)),
               [np.mean(self.result_data_list[method_str]["computation_time"]) for method_str in self.method_list],
               yerr=[np.std(self.result_data_list[method_str]["computation_time"]) for method_str in self.method_list],
               tick_label=self.method_list, ecolor="black", capsize=5, align="center", log=True)
        ax.set_title("Computation time")
        plt.xticks(rotation=20, ha="center")
        ax.set_ylabel("time [ms]".format(axis_str))
        ax.grid(axis="y")

        # Show
        # plt.tight_layout()
        plt.subplots_adjust(left=0.06, bottom=0.06, right=0.98, top=0.96, wspace=0.2, hspace=0.6)
        plt.show()


if __name__ == "__main__":
    plot = PlotTestZmpBasedMethodResults()
