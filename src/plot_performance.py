import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def extract_performance_data(file_path):
    """Extracts GridSize, TotalTime, and AvgCycleIterations from a given file."""
    with open(file_path, "r") as f:
        lines = [line.strip() for line in f.readlines()]
    
    grid_size = int(lines[1])  # GridSize value
    total_time = float(lines[13])  # TotalTime value
    
    cycle_start = lines.index("CycleIterations") + 1
    cycle_iterations = [int(lines[i]) for i in range(cycle_start, len(lines)) if lines[i].isdigit()]
    avg_cycle_iterations = np.mean(cycle_iterations) if cycle_iterations else 0
    
    return grid_size, total_time, avg_cycle_iterations

def plot_performance(data_folder):
    """Reads all performance files and generates the required plots."""
    performance_files = glob.glob(os.path.join(data_folder, "parameters_*_performance.csv"))
    
    if not performance_files:
        print("No performance files found in the specified directory.")
        return
    
    data = [extract_performance_data(f) for f in performance_files]
    df = pd.DataFrame(data, columns=["GridSize", "TotalTime", "AvgCycleIterations"])
    df = df.sort_values(by="GridSize")
    
    # Plot TotalTime vs GridSize
    plt.figure(figsize=(10, 5))
    plt.plot(df["GridSize"], df["TotalTime"], marker='o', linestyle='-', label="Total Time")
    plt.xlabel("Grid Size")
    plt.ylabel("Total Time (s)")
    plt.title("Total Time vs Grid Size")
    plt.legend()
    plt.grid()
    plt.show()
    
    # Plot AvgCycleIterations vs GridSize
    plt.figure(figsize=(10, 5))
    plt.plot(df["GridSize"], df["AvgCycleIterations"], marker='s', linestyle='-', color='r', label="Avg Cycle Iterations")
    plt.xlabel("Grid Size")
    plt.ylabel("Average Cycle Iterations")
    plt.title("Avg Cycle Iterations vs Grid Size")
    plt.legend()
    plt.grid()
    plt.show()

if __name__ == "__main__":
    data_folder = "."  # Change this to the folder containing the performance files
    plot_performance(data_folder)
