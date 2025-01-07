#pragma once

#include "discretization/discretization.h"
#include "discretization/centralDifferences.h"
#include "discretization/donorCell.h"
#include "pressureSolver/pressureSolver.h"
#include "pressureSolver/gaussSeidel.h"
#include "pressureSolver/sor.h"
#include "pressureSolver/cg.h"
#include "output_writer/output_writer.h"
#include "output_writer/output_writer_paraview.h"
#include "output_writer/output_writer_text.h"
#include "output_writer/output_writer_paraview_parallel.h"
#include "output_writer/output_writer_text_parallel.h"
#include "partitioning/partitioning.h"
#include "settings.h"
#include <memory>
#include <iostream>



/**
 * This class handles the main simulation. 
 * It implements the time stepping scheme, computes all the terms and calls the pressure solver.
 */
class Computation {
public:

    /**
     * Initialize the computation object, parse the settings from the file that is given as the only command line argument
     * @param argc number of command line arguments
     * @param argv command line arguments, as given in main
     */
    void initialize(int argc, char* argv[]);
    
    /**
     * run the whole simulation until t_end
     */
    void runSimulation();
    
protected:

    /**
     * Compute the time step width dt from the CFL-condition
     */
    void computeTimeStepWidth();

    /**
     * Set boundary values of u and v
     */
    void applyBoundaryValues();

    /**
     * set initial boundary values, which do not change during the simulation
     */
    void applyInitalBoundaryValues();

    /**
     * set boundary values of preliminary velocities f and g
     */
    void applyPreliminaryBoundaryValues();

    /**
     * Compute the preliminary velocities F and G and set boundary values
     */
    void computePreliminaryVelocities();

    /**
     * Compute the right hand side of the Poisson equation for the pressure
     */
    void computeRightHandSide();


    /**
     * Solve the Poisson equation for the pressure
     */
    void computePressure();


    /**
     * Compute the new velocities u,v from the preliminary velocities F,G and the pressure p
     */
    void computeVelocities();

    /**
     * Compute the new temperature T
     */
    void computeTemperature();

    // set attributes
    Settings settings_;
    std::shared_ptr<Discretization> discretization_;
    std::unique_ptr<PressureSolver> pressureSolver_;
    // std::unique_ptr<OutputWriterParaview> outputWriterParaview_;
    // std::unique_ptr<OutputWriterText> outputWriterText_;
    std::array<double, 2> meshWidth_;
    double dt_;
    std::shared_ptr<Partitioning> partitioning_;
};