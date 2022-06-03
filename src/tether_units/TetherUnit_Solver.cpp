#include "TetherUnit_Solver.hpp"

Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");

TetherUnit_Solver::TetherUnit_Solver(IntegrationInterface* integrator_, IntegrationInterface* integratorsStep_, double mass_distribution_, double tether_length_, unsigned int num_integrationSteps_)

{

    mass_distribution = mass_distribution_; 
    tether_length = tether_length_; 
    num_integrationSteps = num_integrationSteps_;

    integrator = integrator_;

    proximalStates.resize(numStates, 1); 
    proximalStates << 0, 0, 0, 1, 0, 0, 0, mass_distribution * tether_length * g, 0, 0, 0, 0, 0, 0, 0, 0.05, 0;
    distalStates.resize(numStates, 1); 
    fullStates.resize(numStates, num_integrationSteps);

}

TetherUnit_Solver::~TetherUnit_Solver() 

{


}

Eigen::MatrixXd TetherUnit_Solver::integrateDistalStates() 

{

    return integrator->integrate(proximalStates);

}

Eigen::MatrixXd TetherUnit_Solver::integrateFullStates() 

{

    
}


void TetherUnit_Solver::setTau(double tau)

{

    proximalStates(13, 1) = tau;

}

void TetherUnit_Solver::setInitialConditions(Eigen::MatrixXd initialConditions) 

{

    assertm(proximalStates.rows() == initialConditions.rows(), "Initial Conditions must be the correct size!");
    assertm(proximalStates.cols() == initialConditions.cols(), "Initial Conditions must be the correct size!");
   
    proximalStates = initialConditions; 

}

Eigen::MatrixXd TetherUnit_Solver::getProximalStates() {}

Eigen::MatrixXd TetherUnit_Solver::getFullStates() {} 

void TetherUnit_Solver::SolveJacobians() {}