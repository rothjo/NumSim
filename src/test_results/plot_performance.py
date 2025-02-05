import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
from matplotlib.ticker import ScalarFormatter

def extract_performance_data(file_path):
    """Extracts GridSize, TotalTime, PressureSolver, CycleType, and AvgCycleIterations from a given file."""
    with open(file_path, "r") as f:
        lines = [line.strip() for line in f.readlines()]
    
    grid_size = int(lines[lines.index("GridSize") + 1])  # GridSize value
    total_time = float(lines[lines.index("TotalTime") + 1])  # TotalTime value
    pressure_solver = lines[lines.index("PressureSolver") + 1]  # Pressure Solver type
    cycle_type = lines[lines.index("CycleType") + 1] if "CycleType" in lines else ""  # Cycle type (V/W)
    
    cycle_start = lines.index("CycleIterations") + 1
    cycle_iterations = [int(lines[i]) for i in range(cycle_start, len(lines)) if lines[i].isdigit()]
    avg_cycle_iterations = np.mean(cycle_iterations) if cycle_iterations else 0
    
    return grid_size, total_time, pressure_solver, cycle_type, avg_cycle_iterations

def plot_performance(data_folder):
    """Reads all performance files and generates the required plots."""
    performance_files = glob.glob(os.path.join(data_folder, "parameters_*_performance.csv"))
    if not performance_files:
        print("No performance files found in the specified directory.")
        return
    
    data = [extract_performance_data(f) for f in performance_files]
    df = pd.DataFrame(data, columns=["GridSize", "TotalTime", "PressureSolver", "CycleType", "AvgCycleIterations"])
    df = df.sort_values(by="GridSize")
    
    x_ticks = [2**i for i in range(4, 10)]

    # Plot TotalTime vs GridSize with different lines for each PressureSolver
    plt.figure(figsize=(10, 5))
    ax = plt.gca()

    for solver, sub_df in df.groupby("PressureSolver"):
        if solver == "Multigrid":
            for cycle, cycle_df in sub_df.groupby("CycleType"):
                slope, _, _, _, _ = linregress(np.log2(cycle_df["GridSize"]), np.log2(cycle_df["TotalTime"]))
                plt.plot(cycle_df["GridSize"], cycle_df["TotalTime"], marker='o', linestyle='-', label=f"{solver} ({cycle}-cycle), slope={slope:.2f}")
        else:
            slope, _, _, _, _ = linregress(np.log2(sub_df["GridSize"]), np.log2(sub_df["TotalTime"]))
            plt.plot(sub_df["GridSize"], sub_df["TotalTime"], marker='o', linestyle='-', label=f"{solver}, slope={slope:.2f}")
    
    plt.xlabel("Grid Size", fontsize=13)
    plt.ylabel("Total Time (s)", fontsize=13)
    plt.xscale("log", base=2)
    plt.yscale("log")
    plt.xticks(ticks=x_ticks)
    ax.set_xticklabels([str(x) for x in x_ticks])

    handles, labels = ax.get_legend_handles_labels()
    order = [1, 4, 0, 2, 3]
    plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], fontsize=13)

    plt.grid()
    plt.title("Runtime Comparison (Log-Log)", fontsize=16)
    plt.savefig("TotalTime_vs_GridSize.png", dpi=200)
    
    # Plot AvgCycleIterations vs GridSize for Multigrid only (V and W cycle)
    plt.figure(figsize=(10, 5))
    ax = plt.gca()

    multigrid_df = df[df["PressureSolver"] == "Multigrid"]
    for cycle, cycle_df in multigrid_df.groupby("CycleType"):
        plt.plot(cycle_df["GridSize"], cycle_df["AvgCycleIterations"], marker='o', linestyle='-', label=f"Multigrid ({cycle}-cycle)")

    plt.xlabel("Grid Size", fontsize=13)
    plt.ylabel("Average Cycle Iterations", fontsize=13)
    plt.xscale("log", base=2)
    plt.xticks(ticks=x_ticks)
    ax.set_xticklabels([str(x) for x in x_ticks])

    plt.legend(fontsize=13)
    plt.title("Cycle Iterations for Multigrid", fontsize=16)
    plt.grid()
    plt.savefig("Iterations_vs_GridSize.png", dpi=200)

if __name__ == "__main__":
    data_folder = "./runtimeIterationsTest"  # Change this to the folder containing the performance files
    plot_performance(data_folder)
