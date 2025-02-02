#pragma once

#include "computation.h"
#include "partitioning/partitioning.h"
#include "pressureSolver/parallelPressureSolver.h"
#include "pressureSolver/parallelGaussSeidel.h"
#include "pressureSolver/parallel_sor.h"
#include "pressureSolver/parallel_cg.h"
#include "output_writer/output_writer_paraview_parallel.h"
#include "output_writer/output_writer_text_parallel.h"
#include "output_writer/output_writer.h"
#include "mpi.h"

class ParallelComputation : public Computation {
public:

    /**
     * Initialize the parallel computation object, parse the settings from the file that is given as the only command line argument
     * @param argc number of command line arguments
     * @param argv command line arguments, as given in main
     */
    virtual void initialize(int argc, char* argv[]);

    /**
     * run the whole simulation until t_end
     */
    virtual void runSimulation();
    
protected:

    /**
     * Compute the time step width dt from the CFL-condition including communication
     */
    virtual void computeTimeStepWidth();

    /**
     * Set boundary values of u and v or communicate to the neighbouring ranks
     */
    virtual void applyBoundaryValues();

    /**
     * set initial boundary values, which do not change during the simulation
     */
    virtual void applyInitalBoundaryValues();

    /**
     * set boundary values of T or communicate to the neighbouring ranks
     */
    virtual void communicateTemperature();


    std::unique_ptr<OutputWriterParaviewParallel> outputWriterParaview_;
    std::unique_ptr<OutputWriterTextParallel> outputWriterText_;
};