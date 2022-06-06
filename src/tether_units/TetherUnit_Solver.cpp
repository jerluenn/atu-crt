#include "TetherUnit_Solver.hpp"

Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");

TetherUnit_Solver::TetherUnit_Solver(IntegrationInterface* integrator_, IntegrationInterface* integratorStep_, double mass_distribution_, double tether_length_, unsigned int num_integrationSteps_)

{

    mass_distribution = mass_distribution_; 
    tether_length = tether_length_; 
    num_integrationSteps = num_integrationSteps_;

    integrator = integrator_;
    integratorStep = integratorStep_;
    num_p = 3; 
    num_eta = 4;
    num_int_wrench = 6; 
    num_alpha = 1; 
    num_curvature = 1; 
    num_kappa = 1;
    num_tau = 1;

    E_yu.resize(6, 6); 
    B_yu.resize(6, 6); 
    B_yu.setZero();
    B_w_tip.resize(6, 6); 
    E_w_tip.resize(6, 6); 
    J_w_tip.resize(6, 6); 
    yu_dot_w_tip.resize(6, 6);

    dummyStates.resize(numStates, 1);
    tipWrench << 0, 0, 0, 0, 0, 0;
    proximalStates.resize(numStates, 1); 
    proximalStates << 0, 0, 0, 1, 0, 0, 0, mass_distribution * tether_length * g, 0, 0, 0, 0.183972, 0, 0, 0, 0.05, 0;
    distalStates.resize(numStates, 1); 
    distalStates.setZero();
    fullStates.resize(numStates, num_integrationSteps);

}

TetherUnit_Solver::~TetherUnit_Solver() 

{


}

Eigen::MatrixXd TetherUnit_Solver::integrateDistalStates() 

{

    distalStates = integrator->integrate(proximalStates);

    return distalStates;

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

Eigen::MatrixXd TetherUnit_Solver::integrateWithIncrement(unsigned int index) 

{

    dummyStates = proximalStates;
    dummyStates(index) += EPS; 

    return integrator->integrate(dummyStates);

}

Eigen::Matrix<double, 6, 1> TetherUnit_Solver::getBoundaryConditions_IncreaseWrenchTip(unsigned int index) 

{

    Eigen::Matrix<double, 6, 1> boundaryConditions; 
    Eigen::Matrix<double, 6, 1> internalTipWrench;
    Eigen::Matrix<double, 6, 1> tipWrench_plus;
    tipWrench_plus = tipWrench;
    tipWrench_plus(index) += EPS;

    internalTipWrench = distalStates.block<6, 1>(num_eta + num_p , 0);

    boundaryConditions = internalTipWrench - tipWrench_plus;

    return boundaryConditions;

}

Eigen::Matrix<double, 6, 1> TetherUnit_Solver::getBoundaryConditions(Eigen::MatrixXd distalConditions_) 

{

    Eigen::Matrix<double, 6, 1> boundaryConditions; 
    Eigen::Matrix<double, 6, 1> internalTipWrench;

    internalTipWrench = distalConditions_.block<6, 1>(num_eta + num_p , 0);

    boundaryConditions = internalTipWrench - tipWrench;

    return boundaryConditions;

}

Eigen::Matrix<double, 6, 1> TetherUnit_Solver::getBoundaryConditions() 

{

    Eigen::Matrix<double, 6, 1> boundaryConditions; 
    Eigen::Matrix<double, 6, 1> internalTipWrench;

    internalTipWrench = distalStates.block<6, 1>(num_eta + num_p , 0);

    boundaryConditions = internalTipWrench - tipWrench;

    return boundaryConditions;

}

Eigen::Matrix<double, 6, 1> TetherUnit_Solver::setTipWrench(Eigen::Matrix<double, 6, 1> p_FM) 

{

    tipWrench = p_FM;

}

Eigen::Matrix<double, 6, 1> TetherUnit_Solver::getTipWrench() 

{

    return tipWrench;
    
}

Eigen::MatrixXd TetherUnit_Solver::getProximalStates() {

    return proximalStates;

}

Eigen::MatrixXd TetherUnit_Solver::getFullStates() {} 

void TetherUnit_Solver::simulateStep(Eigen::Matrix<double, 6, 1> tip_wrench) 

{

    Eigen::Matrix<double, 6, 1> db; 

    solveJacobians(); 

    // std::cout << "yu_dot_w_tip: " << yu_dot_w_tip << "\n\n";
    // std::cout << "tip wrench: " << tip_wrench << "\n\n";

    db = yu_dot_w_tip * tip_wrench; 

    // std::cout << "db * dt: " << db * dt << "\n\n";
    // std::cout << "db: " << db << "\n\n";

    proximalStates.block<6, 1>(7, 0) += db * dt; 
    tipWrench += tip_wrench * dt ;

    integrateDistalStates();


}

void TetherUnit_Solver::solveJacobians() 

{

    Eigen::MatrixXd distalStates_Increment; 
    distalStates_Increment.resize(numStates, 1);

    Eigen::Matrix<double, 3, 1> position, positionPlus;
    position = distalStates.block<3, 1>(0, 0);

    Eigen::Matrix<double, 4, 1> eta; 
    eta = distalStates.block<4, 1>(3, 0);
    Eigen::Matrix3d rotation, rotationPlus; 
    rotation = MathUtils::quat2Rot(eta);

    Eigen::Matrix<double, 6, 1> BC, BC_Plus; 

    BC = getBoundaryConditions();

    Eigen::MatrixXd linear_twist; 
    linear_twist.resize(3, 1); 
    Eigen::MatrixXd rotation_twist; 
    rotation_twist.resize(3, 3); 
    Eigen::MatrixXd boundary_twist; 
    boundary_twist.resize(6, 1); 

    Eigen::Matrix4d se3; 

    //  std::cout << "BC: " << BC << "\n\n";

    R_dynamic.resize(9, 1);    

    for (unsigned int k = 0; k < 6; ++k) 
    
    {

        distalStates_Increment = integrateWithIncrement(k + num_p + num_eta); 
        positionPlus = distalStates_Increment.block<3, 1>(0, 0);
        eta = distalStates_Increment.block<4, 1>(3, 0);
        rotationPlus = MathUtils::quat2Rot(eta);
        BC_Plus = getBoundaryConditions(distalStates_Increment);

        // std::cout << "BC_Plus: " << BC_Plus << "\n\n"; 
       

        linear_twist = MathUtils::forwardFiniteDifferences(position, positionPlus, EPS); 
        rotation_twist = MathUtils::forwardFiniteDifferences(rotation, rotationPlus, EPS);         
        boundary_twist = MathUtils::forwardFiniteDifferences(BC, BC_Plus, EPS);

        se3.block<3, 3>(0, 0) = rotation_twist*rotation.transpose();
        se3.block<3, 1>(0, 3) = linear_twist;

        E_yu.col(k) = MathUtils::se3toVec(se3);
        B_yu.col(k) = boundary_twist; 

        // BC_Plus = getBoundaryConditions_IncreaseWrenchTip(k); 
        // boundary_twist = MathUtils::forwardFiniteDifferences(BC, BC_Plus, EPS);
        // B_w_tip.col(k) = boundary_twist;

    }

    // std::cout << "B_yu: " << B_yu << "\n\n";
    // std::cout << "B_w_tip: " << B_w_tip << "\n\n";

    Eigen::Matrix<double, 6, 6> B_w_tip_test; 
    B_w_tip_test.setIdentity(); 
    B_w_tip_test *= -1; 

    // std::cout << "B_w_tip_test: " << B_w_tip_test << "\n\n";

    // std::cout << "pinv(B_yu)*B_yu" << B_yu.completeOrthogonalDecomposition().pseudoInverse()*B_yu << "\n\n";

    yu_dot_w_tip = -B_yu.completeOrthogonalDecomposition().pseudoInverse()*B_w_tip_test;

}