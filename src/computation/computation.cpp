#include "computation.h"
#include <cmath>


void Computation::initialize(int argc, char* argv[]) {
    // load settigns from file
    std::string filename = argv[1];
    settings_.loadFromFile(filename);
    partitioning_ = std::make_shared<Partitioning>();
    partitioning_->initialize(settings_.nCells);
    meshWidth_ = {settings_.physicalSize[0] / settings_.nCells[0], settings_.physicalSize[1] / settings_.nCells[1]};

    // init discretization
    if (settings_.useDonorCell) {
        discretization_ = std::make_shared<DonorCell>(settings_.nCells, meshWidth_, partitioning_, settings_.alpha);
    } else {
        discretization_ = std::make_shared<CentralDifferences>(settings_.nCells, meshWidth_, partitioning_);
    }

    // init temperature field
    if (settings_.computeHeat) {
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
            for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
                discretization_->t(i, j) = settings_.initialTemp;
            }
        }
    }

    // init output writers
    outputWriterParaview_ = std::make_unique<OutputWriterParaview>(discretization_, *partitioning_);
    outputWriterText_ = std::make_unique<OutputWriterText>(discretization_, *partitioning_);

    cg_ = std::make_unique<CG>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations);

    // init pressure solvers
    if (settings_.pressureSolver == "SOR") {
        pressureSolver_ = std::make_unique<SOR>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega);
    } else if (settings_.pressureSolver == "GaussSeidel") {
        pressureSolver_ = std::make_unique<GaussSeidel>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations);
    } else if (settings_.pressureSolver == "CG") {
        pressureSolver_ = std::make_unique<CG>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations);
    } else if(settings_.pressureSolver == "Multigrid") {
        pressureSolver_ = std::make_unique<Multigrid>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, partitioning_);
    }
    else {
        std::cerr << "Unknown pressure solver: " << settings_.pressureSolver << std::endl;
        std::exit(1);
    }
}

void Computation::runSimulation() {
    applyInitalBoundaryValues();

    double time = 0.0;
    int t_iter = 0;
    double time_epsilon = 1e-8;
    double output = 0.0;
    
    // Loop over all time steps until t_end is reached
    while (time < (settings_.endTime - time_epsilon)) {
        
        applyBoundaryValues();
        
        computeTimeStepWidth();

        // Check to choose the last time step differently to hit t_end exactly
        if (time + dt_ > settings_.endTime - time_epsilon) {
            dt_ = settings_.endTime - time;
        }
        time += dt_;
        // std::cout << std::endl;
        // std::cout << time << std::endl;

        if (settings_.computeHeat) {
            computeTemperature();
        }

        computePreliminaryVelocities();

        // applyPreliminaryBoundaryValues();

        computeRightHandSide();

        computePressure();

        // if (settings_.pressureSolver == "Multigrid") {
        //     computeRightHandSide();
        // }
        
        // cg_->solve();

        computeVelocities();

        // outputWriterParaview_->writeFile(time); // Output
        // outputWriterText_->writeFile(time); // Output

        // Output
        if (time >= output) {
            // (*outputWriterParaview_).writeFile(time); // Output
            // outputWriterText_->writeFile(time); // Output
            output = output + 0.1;
        }
    }
}

void Computation::computeTimeStepWidth() {
    const double dx = discretization_->dx();
    const double dy = discretization_->dy();

    const double dx2 = dx * dx;
    const double dy2 = dy * dy;

    const double dt_diffusion = (settings_.re / 2.0) * (dx2 * dy2)/(dx2 + dy2);
    const double dt_diffusion_temp = settings_.pr * dt_diffusion;

    const double dt_convection_x = dx / discretization_->u().computeMaxAbs();
    const double dt_convection_y = dy / discretization_->v().computeMaxAbs();
    const double dt_convection = std::min(dt_convection_x, dt_convection_y);

    double dt = settings_.tau * std::min(dt_diffusion, dt_convection);
    if (settings_.computeHeat) {
        dt = std::min(dt, dt_diffusion_temp);
    }

    if (dt > settings_.maximumDt) {
        std::cout << "Warning: Time step width is larger than maximum time step width. Using maximum time step width instead." << std::endl;
        dt_ = settings_.maximumDt;
    } else {
        dt_ = dt;
    }
}

void Computation::applyInitalBoundaryValues() {
    // u and f boundary values
    for (int j = (discretization_->uJBegin() - 1); j < (discretization_->uJEnd() + 1); j++) {
        // left boundary
        discretization_->u(discretization_->uIBegin() - 1, j) = settings_.dirichletBcLeft[0];
        discretization_->f(discretization_->uIBegin() - 1, j) = settings_.dirichletBcLeft[0];

        // right boundary
        discretization_->u(discretization_->uIEnd(), j) = settings_.dirichletBcRight[0];
        discretization_->f(discretization_->uIEnd(), j) = settings_.dirichletBcRight[0];
    }

    // v and g boundary values
    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
        // top boundary
        discretization_->v(i, discretization_->vJEnd()) = settings_.dirichletBcTop[1];
        discretization_->g(i, discretization_->vJEnd()) = settings_.dirichletBcTop[1];

        // bottom boundary
        discretization_->v(i, discretization_->vJBegin() - 1) = settings_.dirichletBcBottom[1];
        discretization_->g(i, discretization_->vJBegin() - 1) = settings_.dirichletBcBottom[1];
    }
}

void Computation::applyBoundaryValues() {
    // u boundary values
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
        // top boundary
        discretization_->u(i, discretization_->uJEnd()) = 2 * settings_.dirichletBcTop[0] - discretization_->u(i, discretization_->uJEnd() - 1);
        
        // bottom boundary
        discretization_->u(i, discretization_->uJBegin() - 1) = 2 * settings_.dirichletBcBottom[0] - discretization_->u(i, discretization_->uJBegin());
    }

    // v boundary values
    for (int j = (discretization_->vJBegin() - 1); j < (discretization_->vJEnd() + 1); j++) {
        // left boundary
        discretization_->v(discretization_->vIBegin() - 1, j) = 2 * settings_.dirichletBcLeft[1] -  discretization_->v(discretization_->vIBegin(), j);
        
        // right boundary
        discretization_->v(discretization_->vIEnd(), j) = 2 * settings_.dirichletBcRight[1] -  discretization_->v(discretization_->vIEnd() - 1, j);
    }

    if (settings_.computeHeat) {
        // Dirichlet at top and bottom boundary
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
            discretization_->t(i, discretization_->pJBegin()-1) = 2*settings_.dirichletBottomTemp - discretization_->t(i, discretization_->pJBegin());
            discretization_->t(i, discretization_->pJEnd()) = 2*settings_.dirichletTopTemp - discretization_->t(i, discretization_->pJEnd()-1);
        }

        // Neumann at left and right boundary
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
            discretization_->t(discretization_->pIBegin()-1, j) = discretization_->t(discretization_->pIBegin(), j);
            discretization_->t(discretization_->pIEnd(), j) = discretization_->t(discretization_->pIEnd()-1, j);
        }
    }
}

void Computation::computePreliminaryVelocities() {
    // compute f
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
        for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
            const double f_diffusion_term = (discretization_->computeD2uDx2(i,j) + discretization_->computeD2uDy2(i,j)) / settings_.re;
            const double f_convection_term = discretization_->computeDu2Dx(i,j) + discretization_->computeDuvDy(i,j);
            discretization_->f(i,j) = discretization_->u(i,j) + dt_*(f_diffusion_term - f_convection_term + settings_.g[0]);
            if (settings_.computeHeat) {
                discretization_->f(i,j) = discretization_->f(i,j) - dt_ * settings_.beta * settings_.g[0] * ((discretization_->t(i,j) + discretization_->t(i+1,j)) / 2.0 - settings_.initialTemp);
            }
            
        }
    }

    // compute g
    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
        for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++) {
            const double g_diffusion_term = (discretization_->computeD2vDx2(i,j) + discretization_->computeD2vDy2(i,j)) / settings_.re;
            const double g_convection_term = discretization_->computeDuvDx(i,j) + discretization_->computeDv2Dy(i,j);
            discretization_->g(i,j) = discretization_->v(i,j) + dt_*(g_diffusion_term - g_convection_term + settings_.g[1]);
            if (settings_.computeHeat) {
                discretization_->g(i,j) = discretization_->g(i,j) - dt_ * settings_.beta * settings_.g[1] * ((discretization_->t(i,j) + discretization_->t(i,j + 1)) / 2.0 - settings_.initialTemp);
            }
        }    
    }
}

void Computation::applyPreliminaryBoundaryValues() {
    // set f at bottom and top boundaries
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
        discretization_->f(i, discretization_->uJEnd()) = discretization_->u(i, discretization_->uJEnd());
        discretization_->f(i, discretization_->uJBegin() - 1) = discretization_->u(i, discretization_->uJBegin() - 1);
    }

    // set g at left and right boundaries
    for (int j = (discretization_->vJBegin() - 1); j < (discretization_->vJEnd() + 1); j++) {
        discretization_->g(discretization_->vIBegin() - 1, j) = discretization_->v(discretization_->vIBegin() - 1, j);
        discretization_->g(discretization_->vIEnd(), j) = discretization_->v(discretization_->vIEnd(), j);
    }
}

void Computation::computeRightHandSide() {
    for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
            const double F_x = (discretization_->f(i, j) - discretization_->f(i - 1, j)) / discretization_->dx();
            const double G_y = (discretization_->g(i, j) - discretization_->g(i, j - 1)) / discretization_->dy();

            discretization_->rhs(i, j) = (F_x + G_y) / dt_;
        }
    }
}

void Computation::computePressure() {
    pressureSolver_->solve();
}

void Computation::computeVelocities() {
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

void Computation::computeTemperature() {
    for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++) {
            discretization_->t(i, j) = discretization_->t(i, j) + dt_ * (1/settings_.re * 1/settings_.pr * (discretization_->computeD2TDx2(i, j) + discretization_->computeD2TDy2(i, j)) - discretization_->computeDuTDx(i, j) - discretization_->computeDvTDy(i, j));
        }
    }
}