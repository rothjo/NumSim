# NumSim
implementation for NumSim

To execute all parameter files in build in a single run: do ./run_simulations.sh in build folder
All parameter files must be named parameters_*.txt

Otherwise run as usual: mpirun -np 1 ./numsim "$param_file" or ./numsim "$param_file"

Multigrid parameters can be adjusted in parameters file, CG refers to PCG with Jacobi