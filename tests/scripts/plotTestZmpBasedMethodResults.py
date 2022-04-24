#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt


class PlotTestZmpBasedMethodResults(object):
    def __init__(self,
                 method,
                 plot_capture_point=False,
                 plot_comp_time=False):
        # Setup plot
        self.fig = plt.figure()
        plt.rcParams["font.size"] = 10
        if isinstance(method, list):
            method_list = method
        else:
            method_list = [method]
        axis_list = ["x", "y"]
        subplot_args = [len(method_list), len(axis_list), 0]
        if plot_comp_time:
            subplot_args[0] += 1

        # Loop for methods
        self.result_data_list = {}
        for method_idx, method_str in enumerate(method_list):
            # Load file
            result_file_path = "/tmp/Test{}.txt".format(method_str)
            print("[PlotTestZmpBasedMethodResults] Load {}".format(result_file_path))
            self.result_data_list[method_str] = np.genfromtxt(result_file_path, dtype=None, delimiter=None, names=True)
            result_data = self.result_data_list[method_str]

            # Loop for axes
            for axis_idx, axis_str in enumerate(axis_list):
                # Setup subplot
                subplot_args[2] = 2 * method_idx + axis_idx + 1
                ax = self.fig.add_subplot(*subplot_args)

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
                ax.set_xlabel("time [s]", labelpad=-2)
                ax.set_ylabel("pos [m]".format(axis_str))
                ax.grid()
                if method_idx == 0 and axis_idx == 0:
                    ax.legend(loc="lower right")

        # Plot computation time
        if plot_comp_time:
            subplot_args[2] = len(method_list) * len(axis_list) + 1
            ax = self.fig.add_subplot(*subplot_args)
            ax.bar(np.arange(len(method_list)),
                   [np.mean(self.result_data_list[method_str]["computation_time"]) for method_str in method_list],
                   yerr=[np.std(self.result_data_list[method_str]["computation_time"]) for method_str in method_list],
                   tick_label=method_list, ecolor="black", capsize=5, align="center", log=True)
            ax.set_title("Computation time")
            plt.xticks(rotation=20, ha="center")
            ax.set_ylabel("time [ms]".format(axis_str))
            ax.grid(axis="y")

        # Show
        # plt.tight_layout()
        if len(method_list) == 1:
            self.fig.set_size_inches((12, 4))
            plt.subplots_adjust(left=0.06, bottom=0.14, right=0.98, top=0.90, wspace=0.2, hspace=0.6)
        else:
            self.fig.set_size_inches((12, 13))
            plt.subplots_adjust(left=0.06, bottom=0.06, right=0.98, top=0.96, wspace=0.2, hspace=0.7)
        plt.show()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Plot test results of ZMP-based methods.")
    all_methods = ["PreviewControlZmp", "DdpZmp", "DcmTracking", "FootGuidedControl", "LinearMpcZmp", "IntrinsicallyStableMpc"]
    parser.add_argument("--method", "-m", type=str, default="PreviewControlZmp",
                        choices=all_methods+["All"])
    parser.add_argument("--plot-capture-point", "-pcp", action="store_true")
    parser.add_argument("--plot-comp-time", "-pct", action="store_true")

    args = parser.parse_args()
    if args.method == "All":
        args.method = all_methods
    plot = PlotTestZmpBasedMethodResults(method=args.method,
                                         plot_capture_point=args.plot_capture_point,
                                         plot_comp_time=args.plot_comp_time)
