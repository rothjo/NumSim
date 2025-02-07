import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
from matplotlib.ticker import ScalarFormatter

def extract_performance_data_lowest_level(file_path):
    """Extracts GridSize, TotalTime, PressureSolver, CycleType, and AvgCycleIterations from a given file."""
    with open(file_path, "r") as f:
        lines = [line.strip() for line in f.readlines()]
    
    grid_size = int(lines[lines.index("GridSize") + 1])  # GridSize value
    total_time = float(lines[lines.index("TotalTime") + 1])  # TotalTime value
    pressure_solver = lines[lines.index("PressureSolver") + 1]  # Pressure Solver type
    cycle_type = lines[lines.index("CycleType") + 1] if "CycleType" in lines else ""  # Cycle type (V/W)
    lowest_level = int(lines[lines.index("LowestLevel") + 1])  # Lowest level in the multigrid hierarchy

    cycle_start = lines.index("CycleIterations") + 1
    cycle_iterations = [int(lines[i]) for i in range(cycle_start, len(lines)) if lines[i].isdigit()]
    avg_cycle_iterations = np.mean(cycle_iterations) if cycle_iterations else 0
    
    return grid_size, total_time, pressure_solver, cycle_type, lowest_level, avg_cycle_iterations

def extract_performance_iterations(file_path):
    """Extracts GridSize, TotalTime, PressureSolver, CycleType, and AvgCycleIterations from a given file."""
    with open(file_path, "r") as f:
        lines = [line.strip() for line in f.readlines()]
    
    grid_size = int(lines[lines.index("GridSize") + 1])  # GridSize value
    pressure_solver = lines[lines.index("PressureSolver") + 1]  # Pressure Solver type
    cycle_type = lines[lines.index("CycleType") + 1] if "CycleType" in lines else ""  # Cycle type (V/W)

    if pressure_solver == "GaussSeidel":
        iterations_start = lines.index("CycleIterations") + 1
        iterations = [int(lines[i]) for i in range(iterations_start, len(lines)) if lines[i].isdigit()]
        total_iterations = sum(iterations)
        timesteps = len(iterations)

    else:
        total_iterations = int(lines[lines.index("totalGSIter") + 1])
        timesteps = int(lines[lines.index("Timesteps") + 1])

    return grid_size, pressure_solver, cycle_type, total_iterations, timesteps

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

def extract_performance_data_residuals(file_path):
    """Extracts GridSize, TotalTime, PressureSolver, CycleType, and AvgCycleIterations from a given file."""
    with open(file_path, "r") as f:
        lines = [line.strip() for line in f.readlines()]
    
    pressure_solver = lines[lines.index("PressureSolver") + 1]  # Pressure Solver type
    cycle_type = lines[lines.index("CycleType") + 1] if "CycleType" in lines else ""  # Cycle type (V/W)

    cycle_start = lines.index("CycleIterations") + 1
    residuals_start = lines.index("Residuals") + 1

    cycle_iterations = [int(lines[i]) for i in range(cycle_start, residuals_start-1) if lines[i].isdigit()]
    residuals = [float(lines[i]) for i in range(residuals_start, len(lines))]
        
    return pressure_solver, cycle_type, cycle_iterations, residuals

def plot_performance_residuals(data_folder):
    """Reads all performance files and generates the required plots."""
    performance_files = glob.glob(os.path.join(data_folder, "parameters_*_performance.csv"))
    if not performance_files:
        print("No performance files found in the specified directory.")
        return
    
    data = [extract_performance_data_residuals(f) for f in performance_files]
    df = pd.DataFrame(data, columns=["PressureSolver", "CycleType", "CycleIterations", "Residuals"])

    plt.figure(figsize=(10, 5))
    ax = plt.gca()

    for solver, sub_df in df.groupby("PressureSolver"):
        if solver == "Multigrid":
            for cycle, cycle_df in sub_df.groupby("CycleType"):
                for _,row in cycle_df.iterrows():
                    plt.plot(row["Residuals"][-row["CycleIterations"][-1]:], label=f"{solver} ({cycle}-cycle)")
        # else:
        #     for _,row in sub_df.iterrows():
        #         plt.plot(row["Residuals"][-row["CycleIterations"][-1]:], label=f"{solver}")
        
    plt.xlabel("Iterations", fontsize=13)
    plt.ylabel("Residuals", fontsize=13)
    # plt.xscale("log")
    plt.yscale("log")

    # handles, labels = ax.get_legend_handles_labels()
    # order = [1, 4, 0, 2, 3]
    # plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], fontsize=13)
    plt.legend(fontsize=10)
    plt.grid()
    plt.title("Residuals, Grid: 512x512", fontsize=16)
    plt.savefig("residuals_multigrid.png", dpi=200)
    
def plot_performance_lowest_level(data_folder):
    """Reads all performance files and generates the required plots."""
    performance_files = glob.glob(os.path.join(data_folder, "parameters_*_performance.csv"))
    if not performance_files:
        print("No performance files found in the specified directory.")
        return
    
    data = [extract_performance_data_lowest_level(f) for f in performance_files]
    df = pd.DataFrame(data, columns=["GridSize", "TotalTime", "PressureSolver", "CycleType", "LowestLevel", "AvgCycleIterations"])
    df = df.sort_values(by=["GridSize", "LowestLevel"])
    
    x_ticks = sorted(df["GridSize"].unique())

    # Plot TotalTime vs GridSize with different lines for each PressureSolver
    plt.figure(figsize=(7, 5))
    ax = plt.gca()

    multigrid_df = df[df["PressureSolver"] == "Multigrid"]
    v_df = multigrid_df[multigrid_df["CycleType"] == "V"]
    w_df = multigrid_df[multigrid_df["CycleType"] == "W"]

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    plt.gca().set_prop_cycle(color=colors)

    for (cycle, lowest_level), cycle_df in v_df.groupby(["CycleType", "LowestLevel"]):
        plt.plot(cycle_df["GridSize"], cycle_df["TotalTime"], marker='o', linestyle='-', label=f"L={lowest_level}")
        
    log_log = True

    plt.xlabel("Grid Size", fontsize=13)
    plt.ylabel("Total Time (s)", fontsize=13)
    
    if log_log:
        plt.xscale("log", base=2)
        plt.yscale("log")
    plt.xticks(ticks=x_ticks)
    ax.set_xticklabels([str(x) for x in x_ticks])

    # handles, labels = ax.get_legend_handles_labels()
    # order = [1, 4, 0, 2, 3]
    # plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], fontsize=13)
    plt.legend(fontsize=10)
    plt.grid()
    
    if log_log:
        plt.title("Cycle Depth Comparison V-cycle (Log-Log)", fontsize=16)
        plt.savefig("depth_v_log_log.png", dpi=200)
    else:
        plt.title("Cycle Depth Comparison V-cycle", fontsize=16)
        plt.savefig("depth_v.png", dpi=200)

    plt.figure(figsize=(7, 5))
    ax = plt.gca()

    multigrid_df = df[df["PressureSolver"] == "Multigrid"]
    v_df = multigrid_df[multigrid_df["CycleType"] == "V"]
    w_df = multigrid_df[multigrid_df["CycleType"] == "W"]

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    plt.gca().set_prop_cycle(color=colors)

    for (cycle, lowest_level), cycle_df in w_df.groupby(["CycleType", "LowestLevel"]):
        slope, _, _, _, _ = linregress(np.log2(cycle_df["GridSize"]), np.log2(cycle_df["TotalTime"]))
        plt.plot(cycle_df["GridSize"], cycle_df["TotalTime"], marker='o', linestyle='-', label=f"L={lowest_level}")
    
    
    plt.xlabel("Grid Size", fontsize=13)
    plt.ylabel("Total Time (s)", fontsize=13)
    if log_log:
        plt.xscale("log", base=2)
        plt.yscale("log")
    plt.xticks(ticks=x_ticks)
    ax.set_xticklabels([str(x) for x in x_ticks])

    # handles, labels = ax.get_legend_handles_labels()
    # order = [1, 4, 0, 2, 3]
    # plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], fontsize=13)
    plt.legend(fontsize=10)
    plt.grid()
    if log_log:
        plt.title("Cycle Depth Comparison W-cycle (Log-Log)", fontsize=16)
        plt.savefig("depth_w_log_log.png", dpi=200)
    else:
        plt.title("Cycle Depth Comparison W-cycle", fontsize=16)
        plt.savefig("depth_w.png", dpi=200)

    
# Backup plot_performance function
def plot_performance(data_folder):
    """Reads all performance files and generates the required plots."""
    performance_files = glob.glob(os.path.join(data_folder, "parameters_*_performance.csv"))
    if not performance_files:
        print("No performance files found in the specified directory.")
        return
    
    data = [extract_performance_data(f) for f in performance_files]
    df = pd.DataFrame(data, columns=["GridSize", "TotalTime", "PressureSolver", "CycleType", "AvgCycleIterations"])
    df = df.sort_values(by="GridSize")
    
    x_ticks = sorted(df["GridSize"].unique())
    
    # Plot TotalTime vs GridSize with different lines for each PressureSolver
    solver_cycle_colors = {
        "CG": "tab:blue",
        ("Multigrid", "V"): "tab:orange",
        ("Multigrid", "W"): "tab:green",
        "SOR": "tab:purple",
        "GaussSeidel": "tab:red"
    }

    plt.figure(figsize=(7, 5))
    ax = plt.gca()

    df_small = df[df["GridSize"] <= 128]
    
    for solver, sub_df in df_small.groupby("PressureSolver"):
        if solver == "Multigrid":
            for cycle, cycle_df in sub_df.groupby("CycleType"):
                slope, _, _, _, _ = linregress(np.log(cycle_df["GridSize"]), np.log(cycle_df["TotalTime"]))
                color = solver_cycle_colors.get((solver, cycle), "black")
                plt.plot(cycle_df["GridSize"], cycle_df["TotalTime"], marker='o', color=color, linestyle='-', label=f"{solver} ({cycle}-cycle)")
        else:
            slope, yaxis, _, _, _ = linregress(np.log(sub_df["GridSize"]), np.log(sub_df["TotalTime"]))
            color = solver_cycle_colors.get(solver, "black")
            plt.plot(sub_df["GridSize"], sub_df["TotalTime"], marker='o', color=color, linestyle='-', label=f"{solver}")
    
    plt.xlabel("Grid Size", fontsize=13)
    plt.ylabel("Total Time (s)", fontsize=13)
    # plt.xscale("log", base=2)
    # plt.yscale("log", base=10)

    x_ticks_small = [16, 32, 64, 128]
    plt.xticks(ticks=x_ticks_small)
    ax.set_xticklabels([str(x) for x in x_ticks_small])

    # handles, labels = ax.get_legend_handles_labels()
    # order = [1, 4, 0, 2, 3]
    # plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], fontsize=13)
    plt.legend()
    plt.grid()
    plt.title("Runtime Comparison", fontsize=16)
    plt.savefig("TotalTime_vs_GridSize_small.png", dpi=200)
    

    plt.figure(figsize=(7, 5))
    ax = plt.gca()

    df_big = df[df["PressureSolver"].isin(["CG", "Multigrid"])]
    for solver, sub_df in df_big.groupby("PressureSolver"):
        if solver == "Multigrid":
            for cycle, cycle_df in sub_df.groupby("CycleType"):
                slope, _, _, _, _ = linregress(np.log(cycle_df["GridSize"]), np.log(cycle_df["TotalTime"]))
                color = solver_cycle_colors.get((solver, cycle), "black")
                plt.plot(cycle_df["GridSize"], cycle_df["TotalTime"], marker='o', color=color, linestyle='-', label=f"{solver} ({cycle}-cycle)")
        else:
            last_index = sub_df.index[-1]  # Get the last row index
            slope, yaxis, _, _, _ = linregress(np.log(sub_df["GridSize"]), np.log(sub_df["TotalTime"]))
            color = solver_cycle_colors.get(solver, "black")
            plt.plot(sub_df["GridSize"], sub_df["TotalTime"], marker='o', color=color, linestyle='-', label=f"{solver}")
            
            # Highlight the last unsure point
            plt.scatter(sub_df.loc[last_index, "GridSize"], sub_df.loc[last_index, "TotalTime"], 
                        edgecolor='red', facecolor='none', s=100)
            
            # Add text annotation
            plt.text(sub_df.loc[last_index, "GridSize"], sub_df.loc[last_index, "TotalTime"] * 1.1,
                    "Did not finish", color='red', fontsize=10)    

    plt.xlabel("Grid Size", fontsize=13)
    plt.ylabel("Total Time (s)", fontsize=13)
    # plt.xscale("log", base=2)
    # plt.yscale("log", base=10)
    plt.xticks(ticks=x_ticks)
    ax.set_xticklabels([str(x) for x in x_ticks])

    # handles, labels = ax.get_legend_handles_labels()
    # order = [1, 4, 0, 2, 3]
    # plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], fontsize=13)
    plt.legend()
    plt.grid()
    plt.title("Runtime Comparison", fontsize=16)
    plt.savefig("TotalTime_vs_GridSize_big.png", dpi=200)


def plot_performance_iterations(data_folder):
    """Reads all performance files and generates the required plots."""
    performance_files = glob.glob(os.path.join(data_folder, "parameters_*_performance.csv"))
    if not performance_files:
        print("No performance files found in the specified directory.")
        return

    data = [extract_performance_iterations(f) for f in performance_files]
    df = pd.DataFrame(data, columns=["GridSize", "PressureSolver", "CycleType", "TotalIterations", "Timesteps"])
    df = df.sort_values(by="GridSize")
    
    x_ticks = sorted(df["GridSize"].unique())


    # Plot AvgCycleIterations vs GridSize for Multigrid only (V and W cycle)
    plt.figure(figsize=(7, 5))
    ax = plt.gca()


    for solver, sub_df in df.groupby("PressureSolver"):
        if solver == "Multigrid":
            for cycle, cycle_df in sub_df.groupby("CycleType"):
                # slope,_,_,_,_ = linregress(np.log2(cycle_df["GridSize"]), np.log2(cycle_df["AvgCycleIterations"]))
                plt.plot(cycle_df["GridSize"], cycle_df["TotalIterations"]/cycle_df["Timesteps"], marker='o', linestyle='-', label=f"{solver} ({cycle}-cycle)")
        # else:
        #     # slope,_,_,_,_ = linregress(np.log2(sub_df["GridSize"]), np.log2(sub_df["AvgCycleIterations"]))
        #     plt.plot(sub_df["GridSize"], sub_df["TotalIterations"]/sub_df["Timesteps"], marker='o', linestyle='-', label=f"{solver}")
    
    plt.xlabel("Grid Size", fontsize=13)
    plt.ylabel("Iterations per Timestep", fontsize=13)
    plt.xscale("log", base=2)
    # plt.yscale("log")
    plt.xticks(ticks=x_ticks)
    ax.set_xticklabels([str(x) for x in x_ticks])

    plt.legend(fontsize=13)
    plt.title("Gauss Seidel Iterations per Timestep", fontsize=16)
    plt.grid()
    plt.savefig("GS_Iterations_multigrid.png", dpi=200)


if __name__ == "__main__":
    # data_folder = "./lowest_level"  # Change this to the folder containing the performance files
    # plot_performance_lowest_level(data_folder)

    # data_folder = "./runtimeIterationsTest"
    # plot_performance(data_folder)

    # data_folder = "./final_iterations"
    # plot_performance_iterations(data_folder)

    data_folder = "./residualOverIterations"
    plot_performance_residuals(data_folder)


