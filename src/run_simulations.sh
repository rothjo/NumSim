#!/bin/bash

# Loop through all parameter files in the current directory
for param_file in parameters_*.txt; do
    echo "Running simulation with $param_file..."
    mpirun -np 1 ./numsim "$param_file"
    echo "Finished simulation with $param_file."
done

echo "All simulations completed."