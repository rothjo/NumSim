#include "computation.h"
#include <cmath>

/**
 * Initialize the computation object, parse the settings from the file that is given as the only command line argument
 * 
 */
void Computation::initialize(int argc, char* argv[]) {
    // load settigns from file
    // settings_ = Settings();
    std::string filename = argv[1];
    settings_.loadFromFile(filename);

    // init meshWidth
    meshWidth_ = {settings_.physicalSize[0] / settings_.nCells[0], settings_.physicalSize[1] / settings_.nCells[1]};

    // init discretization
    if (settings_.useDonorCell) {
        discretization_ = std::make_shared<DonorCell>(settings_.nCells, meshWidth_, settings_.alpha);
    } else {
        discretization_ = std::make_shared<CentralDifferences>(settings_.nCells, meshWidth_);
    }

    // init output writers
    outputWriterParaview_ = std::make_unique<OutputWriterParaview>(discretization_);
    outputWriterText_ = std::make_unique<OutputWriterText>(discretization_);

    // init pressure solvers
    if (settings_.pressureSolver == "SOR") {
        pressureSolver_ = std::make_unique<SOR>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega);
    } else if (settings_.pressureSolver == "GaussSeidel") {
        pressureSolver_ = std::make_unique<GaussSeidel>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations);
    } else {
        std::cerr << "Unknown pressure solver: " << settings_.pressureSolver << std::endl;
        std::exit(1);
    }
}


/**
 * run the whole simulation until t_end
 */
void Computation::runSimulation() {
    double time = 0.0;
    int t_iter = 0;
    double time_epsilon = 1e-8;

    // Loop over all time steps until t_end is reached
    while (time < (settings_.endTime - time_epsilon)) {

        applyBoundaryValues();
        
        computeTimeStepWidth();
        // Check to choose the last time step differently to hit t_end
        if (time + dt_ > settings_.endTime - time_epsilon) {
            dt_ = settings_.endTime - time;
        }
        time += dt_;

        computePreliminaryVelocities();

        computeRightHandSide();

        computePressure();

        computeVelocities();

        outputWriterParaview_->writeFile(time); // Output
        outputWriterText_->writeFile(time); // Output
    }


}

void Computation::computeTimeStepWidth() {
    // Compute the time step width dt from maximum velocities
    // Idea: CFL-condition satisfied with factor 0.5
    const double dx2 = discretization_->dx() * discretization_->dx();
    const double dy2 = discretization_->dy() * discretization_->dy();

    // Compute CFL condition from diffusion operator
    const double dt_diffusion = (settings_.re / 2.0) * (dx2 * dy2)/(dx2 + dy2);

    // Compute CFL condition from convection operator
    // Compute u_max and v_max

    // Compute dt_convection
    const double dt_convection_x = discretization_->dx() / discretization_->u().computeMaxAbs();
    const double dt_convection_y = discretization_->dy() / discretization_->v().computeMaxAbs();
    const double dt_convection = std::min(dt_convection_x, dt_convection_y);

    const double dt = settings_.tau * std::min(dt_diffusion, dt_convection);
    if ( dt > settings_.maximumDt) {
        std::cout << "Warning: Time step width is larger than maximum time step width. Using maximum time step width instead." << std::endl;
        dt_ = settings_.maximumDt;
    } else {
        dt_ = dt;
    }
}

void Computation::applyBoundaryValues() {
    // Set boundary values of u, v and for F and G to correct values
    // F only needs to be set on left and right boundary
    // G only needs to be set on top and bottom boundary

    
    // Set top boundary 
    // velocity u
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
        discretization_->u(i, discretization_->uJEnd()) = 2 * settings_.dirichletBcTop[0] - discretization_->u(i, discretization_->uJEnd() - 1);
        // discretization_->f(i, discretization_->uJEnd()) = discretization_->u(i, discretization_->uJEnd());
    }
    // velocity v
    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
        discretization_->v(i, discretization_->vJEnd()) = settings_.dirichletBcTop[1];
        // discretization_->g(i, discretization_->vJEnd()) = discretization_->v(i, discretization_->vJEnd());
    }

    // Set bottom boundary
    // velocity u
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
        discretization_->u(i, discretization_->uJBegin() - 1) = 2 * settings_.dirichletBcBottom[0] - discretization_->u(i, discretization_->uJBegin());
        // discretization_->f(i, discretization_->uJBegin() - 1) = discretization_->u(i, discretization_->uJBegin() - 1);
    }
    // velocity v and preliminary velocity G
    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
        discretization_->v(i, discretization_->vJBegin() - 1) = settings_.dirichletBcBottom[1];
        // discretization_->g(i, discretization_->vJBegin() - 1) = discretization_->v(i, discretization_->vJBegin() - 1);
    }

    // Set left boundary
    // velocity u and preliminary velocity F
    for (int j = (discretization_->uJBegin() - 1); j < (discretization_->uJEnd() + 1); j++) {
        discretization_->u(discretization_->uIBegin() - 1, j) = settings_.dirichletBcLeft[0];
        // discretization_->f(discretization_->uIBegin() - 1, j) = discretization_->u(discretization_->uIBegin() - 1, j);
    }
    // velocity v
    for (int j = (discretization_->vJBegin() - 1); j < (discretization_->vJEnd() + 1); j++) {
        discretization_->v(discretization_->vIBegin() - 1, j) = 2 * settings_.dirichletBcLeft[1] -  discretization_->v(discretization_->vIBegin(), j);
        // discretization_->g(discretization_->vIBegin() - 1, j) = discretization_->v(discretization_->vIBegin() - 1, j);
    }

    // Set right boundary
    // velocity u and preliminary velocity F
    for (int j = (discretization_->uJBegin() - 1); j < (discretization_->uJEnd() + 1); j++) {
        discretization_->u(discretization_->uIEnd(), j) = settings_.dirichletBcRight[0];
        // discretization_->f(discretization_->uIEnd(), j) = discretization_->u(discretization_->uIEnd(), j);
    }
    // velocity v
    for (int j = (discretization_->vJBegin() - 1); j < (discretization_->vJEnd() + 1); j++) {
        discretization_->v(discretization_->vIEnd(), j) = 2 * settings_.dirichletBcRight[1] -  discretization_->v(discretization_->vIEnd() - 1, j);
        // discretization_->g(discretization_->vIEnd(), j) = discretization_->v(discretization_->vIEnd(), j);
    }
}

void Computation::computePreliminaryVelocities() {
    // Compute the preliminary velocities, F and G
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
        for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
            // Compute F
            double f_diffusion_term = (discretization_->computeD2uDx2(i,j) + discretization_->computeD2uDy2(i,j))/settings_.re;
            double f_convection_term = (discretization_->computeDu2Dx(i,j) + discretization_->computeDuvDy(i,j));
            discretization_->f(i,j) = discretization_->u(i,j) + dt_*(f_diffusion_term - f_convection_term + settings_.g[0]);
        }
    }
    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
        for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++) {
            // Compute G
            double g_diffusion_term = (discretization_->computeD2vDx2(i,j) + discretization_->computeD2vDy2(i,j))/settings_.re;
            double g_convection_term = (discretization_->computeDuvDx(i,j) + discretization_->computeDv2Dy(i,j));
            discretization_->g(i,j) = discretization_->v(i,j) + dt_*(g_diffusion_term - g_convection_term + settings_.g[1]);
        }
    }    



    // test: apply prelim BC
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
        discretization_->f(i, discretization_->uJEnd()) = discretization_->u(i, discretization_->uJEnd());
    }
    // velocity v
    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
        discretization_->g(i, discretization_->vJEnd()) = discretization_->v(i, discretization_->vJEnd());
    }

    // Set bottom boundary
    // velocity u
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
        discretization_->f(i, discretization_->uJBegin() - 1) = discretization_->u(i, discretization_->uJBegin() - 1);
    }
    // velocity v and preliminary velocity G
    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
        discretization_->g(i, discretization_->vJBegin() - 1) = discretization_->v(i, discretization_->vJBegin() - 1);
    }

    // Set left boundary
    // velocity u and preliminary velocity F
    for (int j = (discretization_->uJBegin() - 1); j < (discretization_->uJEnd() + 1); j++) {
        discretization_->f(discretization_->uIBegin() - 1, j) = discretization_->u(discretization_->uIBegin() - 1, j);
    }
    // velocity v
    for (int j = (discretization_->vJBegin() - 1); j < (discretization_->vJEnd() + 1); j++) {
        discretization_->g(discretization_->vIBegin() - 1, j) = discretization_->v(discretization_->vIBegin() - 1, j);
    }

    // Set right boundary
    // velocity u and preliminary velocity F
    for (int j = (discretization_->uJBegin() - 1); j < (discretization_->uJEnd() + 1); j++) {
        discretization_->f(discretization_->uIEnd(), j) = discretization_->u(discretization_->uIEnd(), j);
    }
    // velocity v
    for (int j = (discretization_->vJBegin() - 1); j < (discretization_->vJEnd() + 1); j++) {
        discretization_->g(discretization_->vIEnd(), j) = discretization_->v(discretization_->vIEnd(), j);
    } 
}

void Computation::computeRightHandSide() {
    // Compute the right hand side of the Poisson equation for the pressure
    for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {

            const double F_x = ((discretization_->f(i, j) - discretization_->f(i - 1, j)) / discretization_->dx());
            const double G_y = ((discretization_->g(i, j) - discretization_->g(i, j - 1)) / discretization_->dy());

            discretization_->rhs(i, j) = ((F_x + G_y) / dt_);
        }
    }
}

void Computation::computePressure() {
    // Solve the Poisson equation for the pressure
    pressureSolver_->solve();
}

void Computation::computeVelocities() {
    // Compute the new velocities, u,v from the preliminary velocities, F,G and the pressure, p
    
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
        for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
            discretization_->u(i, j) = discretization_->f(i, j) - dt_ * discretization_->computeDpDx(i, j);
        }
    }

    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
        for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++) {
            discretization_->v(i, j) = discretization_->g(i, j) - dt_ * discretization_->computeDpDy(i, j);
        }
    }
}