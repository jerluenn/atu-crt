#include "TetherUnit_Solver.hpp"

Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");

TetherUnit_Solver::TetherUnit_Solver(IntegrationInterface* integrator_, IntegrationInterface* integratorStep_, double mass_distribution_, double tether_length_, unsigned int num_integrationSteps_, unsigned int num_tetherParts_, double lambdaLyapunov_, double lambdaDLS_, double Kp_, Eigen::MatrixXd initialSolution)

{

    mass_distribution = mass_distribution_; 
    tether_length = tether_length_; 
    num_integrationSteps = num_integrationSteps_;
    num_tetherParts = num_tetherParts_;
    total_tether_length = tether_length * num_tetherParts;

    integrator = integrator_;
    integratorStep = integratorStep_;
    num_p = 3; 
    num_eta = 4;
    num_int_wrench = 6; 
    num_alpha = 1; 
    num_curvature = 1; 
    num_kappa = 1;
    num_tau = 1;

    Kp = Kp_;
    lambdaLyapunov = lambdaLyapunov_;
    lambdaDLS = lambdaDLS_;
    initialise();
    integrateDistalStates();
    solveJacobians();
    poseError << 1, 1, 1, 1, 1, 1, 1;

}

TetherUnit_Solver::~TetherUnit_Solver() 

{


}

void TetherUnit_Solver::initialise(Eigen::MatrixXd initialSolution) 

{

    Eigen::MatrixXd zeros6x6(6, 6), zeros7x6(7, 6), zerosStates(numStates, 1); 
    zeros6x6.setZero(); 
    zeros7x6.setZero();
    zeroStates.setZero();
    Eigen::MatrixXd E_yu_tmp(6,6), B_yu_tmp(6,6), B_w_tip_tmp(6,6), E_w_tip_tmp(6,6);
    Eigen::MatrixXd J_w_tip_tmp(6,6), J_w_tip_eta_tmp(7,6), w_tip_dot_tmp(6,6);
    Eigen::MatrixXd yu_dot_w_tip_tmp(6,6), distalStates_tmp(numStates, 1), fullStates(num_integrationSteps*num_tetherParts + 1, numStates); 
    Eigen::MatrixXd proximalStates_tmp(numStates, 1)

    tipWrench << 0, 0, 0, 0, 0, 0;
    I_7x7.setIdentity();
    assertm(proximalStates.rows() == initialSolution.rows(), "Proximal states and initial solution must be equal in size.");
    assertm(Kp_ > 0, "Kp must be positive.");
    assertm(lambdaDLS_ > 0, "lambdaDLS must be positive.");
    assertm(lambdaLyapunov_ > 0, "lambdaLyapunov must be positive.");
    
    eta_dot.resize(4, 6);
    dummyStates.resize(numStates, 1);

    for (int i = 0; i < num_tetherParts; ++i) 
    
    {

        E_yu.push_back(zeros6x6);
        B_yu.push_back(zeros6x6);
        B_w_tip.push_back(zeros6x6);
        E_w_tip.push_back(zeros6x6); 
        J_w_tip.push_back(zeros6x6); 
        J_w_tip_eta.push_back(zeros7x6);
        w_tip_dot.push_back(zeros6x6);
        yu_dot_w_tip.push_back(zeros6x6);
        distalStates.push_back(zeroStates);
        proximalStates.push_back(zeroStates);

    }

    proximalStates[0] << initialSolution; 

}

void TetherUnit_Solver::integrateDistalStates(unsigned int stage_num) 

{

    distalStates[stage_num] = integrator->integrate(proximalStates[stage_num]);

}

Eigen::Matrix<double, 7, 1> TetherUnit_Solver::getDistalPose(unsigned int stage_num) 

{

    return distalStates[stage_num].block<7, 1>(0, 0);

}

Eigen::MatrixXd TetherUnit_Solver::getDistalStates(unsigned int stage_num) 

{

    return distalStates[stage_num];

}

Eigen::Matrix<double, 7, 6> TetherUnit_Solver::getJacobianEta_wrt_tip(unsigned int stage_num) 

{

    return J_w_tip_eta[stage_num];

}

unsigned int TetherUnit_Solver::getNumIntegrationSteps() 

{

    return num_integrationSteps;

}

void TetherUnit_Solver::integrateFullStates() 

{

    // Meant for plotting entire arm.

    Eigen::MatrixXd states_i(numStates, 1);
    states_i = proximalStates;
    fullStates.block<1, 17>(0, 0) = states_i.transpose();  

    for (unsigned int i = 0; i < num_integrationSteps; ++i) 
    
    {

        states_i = integratorStep->integrate(states_i);
        fullStates.block<1, 17>(i + 1, 0) = states_i.transpose();

    }

}

void TetherUnit_Solver::setTau(double tau)

{

    proximalStates[0](13, 1) = tau;

}

void TetherUnit_Solver::setInitialConditions(Eigen::MatrixXd initialConditions) 

{

    assertm(proximalStates.rows() == initialConditions.rows(), "Initial Conditions must be the correct size!");
    assertm(proximalStates.cols() == initialConditions.cols(), "Initial Conditions must be the correct size!");
   
    proximalStates[0] = initialConditions; 

}

void TetherUnit_Solver::setInitialConditions(Eigen::Matrix<double, 6, 1> initialConditions) 

{
   
    proximalStates[0].block<6, 1>(7, 0) = initialConditions; 

}

Eigen::MatrixXd TetherUnit_Solver::integrateWithIncrement(unsigned int index, unsigned int stage_num) 

{

    dummyStates = proximalStates[stage_num];
    dummyStates(index) += EPS; 

    return integrator->integrate(dummyStates);

}

Eigen::Matrix<double, 6, 1> TetherUnit_Solver::getBoundaryConditions_IncreaseWrenchTip(unsigned int index, unsigned int stage_num) 

{

    Eigen::Matrix<double, 6, 1> boundaryConditions; 
    Eigen::Matrix<double, 6, 1> internalTipWrench;
    Eigen::Matrix<double, 6, 1> tipWrench_plus;
    tipWrench_plus = tipWrench;
    tipWrench_plus(index) += EPS;

    internalTipWrench = distalStates[stage_num].block<6, 1>(num_eta + num_p , 0);

    boundaryConditions = internalTipWrench - tipWrench_plus;

    return boundaryConditions;

}

Eigen::Matrix<double, 6, 1> TetherUnit_Solver::getBoundaryConditions(Eigen::MatrixXd distalConditions_, unsigned int stage_num) 

{

    Eigen::Matrix<double, 6, 1> boundaryConditions; 
    Eigen::Matrix<double, 6, 1> internalTipWrench;
    internalTipWrench = distalConditions_.block<6, 1>(num_eta + num_p , 0);

    if (stage_num == num_tetherParts) 
    
    {

        boundaryConditions = internalTipWrench - tipWrench;        

    } 

    else 
    
    {

        boundaryConditions = internalTipWrench - proximalStates[stage_num + 1].block<6, 1>(num_eta + num_p , 0);

    }



    return boundaryConditions;

}


Eigen::Matrix<double, 7, 1> TetherUnit_Solver::getPoseError() 

{

    return poseError;

}

void TetherUnit_Solver::solveReactionForcesStep(Eigen::MatrixXd poseDesired) 

{

    assertm(poseDesired.rows() == 7, "poseDesired must have 7 rows.");
    assertm(poseDesired.cols() == 1, "poseDesired must have 1 col.");

    poseError = poseDesired - getDistalPose();
    tipWrenchInput = Kp * J_w_tip_eta.transpose() * (J_w_tip_eta * J_w_tip_eta.transpose() + pow(lambdaDLS, 2)*I_7x7).inverse() * poseError;
    simulateStep(tipWrenchInput);

}

Eigen::Matrix<double, 6, 1> TetherUnit_Solver::getBoundaryConditions(unsigned int stage_num) 

{

    Eigen::Matrix<double, 6, 1> boundaryConditions; 
    Eigen::Matrix<double, 6, 1> internalTipWrench;

    internalTipWrench = distalStates[stage_num].block<6, 1>(num_eta + num_p , 0);

    boundaryConditions = internalTipWrench - tipWrench;

    return boundaryConditions;

}

void TetherUnit_Solver::setTipWrench(Eigen::Matrix<double, 6, 1> p_FM) 

{

    tipWrench = p_FM;

}

void TetherUnit_Solver::saveData(std::string fileName, Eigen::MatrixXd matrix)
{
    //https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
 
    std::ofstream file(fileName);
    if (file.is_open())
    {
        file << matrix.format(CSVFormat);
        file.close();
    }
}

Eigen::Matrix<double, 6, 1> TetherUnit_Solver::getTipWrench() 

{

    return tipWrench;
    
}

Eigen::MatrixXd TetherUnit_Solver::getProximalStates(unsigned int stage_num) {

    return proximalStates[stage_num];

}

Eigen::MatrixXd TetherUnit_Solver::getFullStates(std::string fileName) 

{

    saveData(fileName, fullStates);

    return fullStates;

} 

Eigen::MatrixXd TetherUnit_Solver::getFullStates() 

{

    return fullStates;

} 

void TetherUnit_Solver::simulateStep(Eigen::Matrix<double, 6, 1> tip_wrench) 

{

    Eigen::Matrix<double, 6, 1> d_yu;  

    // Need to solve Jacobians here iteratively, and then finally find out d_yu by 
    // solving using Jacobians backwards. Use Lyapunov step at every step to make sure 
    // that error remains very close to zero. 

    solveJacobians(); 

    d_yu = yu_dot_w_tip * tip_wrench; 
    // d_yu -= lambdaLyapunov*B_yu.completeOrthogonalDecomposition().pseudoInverse()*(getBoundaryConditions() - (1/lambdaLyapunov)*tip_wrench*dt);

    proximalStates.block<6, 1>(7, 0) += d_yu * dt; 
    tipWrench += tip_wrench * dt ;

    integrateDistalStates();

    d_yu = -lambdaLyapunov*B_yu.completeOrthogonalDecomposition().pseudoInverse()*(getBoundaryConditions() - (1/lambdaLyapunov)*tip_wrench*dt);

    proximalStates.block<6, 1>(7, 0) += d_yu * dt; 
    integrateDistalStates();


}


void TetherUnit_Solver::solveJacobians(unsigned int stage_num) 

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

    Eigen::MatrixXd conv(4, 3);
    Eigen::Matrix3d I_3x3; 
    I_3x3.setIdentity();

    Eigen::Matrix<double, 6, 6> I_6x6; 
    I_6x6.setIdentity();

    BC = getBoundaryConditions(stage_num);

    Eigen::MatrixXd linear_twist; 
    linear_twist.resize(3, 1); 
    Eigen::MatrixXd rotation_twist; 
    rotation_twist.resize(3, 3); 
    Eigen::MatrixXd boundary_twist; 
    boundary_twist.resize(6, 1); 

    Eigen::Matrix4d se3; 


    R_dynamic.resize(9, 1);    

    for (unsigned int k = 0; k < 6; ++k) 
    
    {

        distalStates_Increment = integrateWithIncrement(k + num_p + num_eta, stage_num); 
        positionPlus = distalStates_Increment.block<3, 1>(0, 0);
        eta = distalStates_Increment.block<4, 1>(3, 0);
        rotationPlus = MathUtils::quat2Rot(eta);
        BC_Plus = getBoundaryConditions(distalStates_Increment, stage_num);

        linear_twist = MathUtils::forwardFiniteDifferences(position, positionPlus, EPS); 
        rotation_twist = MathUtils::forwardFiniteDifferences(rotation, rotationPlus, EPS);         
        boundary_twist = MathUtils::forwardFiniteDifferences(BC, BC_Plus, EPS);

        se3.block<3, 3>(0, 0) = rotation_twist*rotation.transpose();
        se3.block<3, 1>(0, 3) = linear_twist;

        E_yu[stage_num].col(k) = MathUtils::se3toVec(se3);
        B_yu[stage_num].col(k) = boundary_twist; 


    }

    yu_dot_w_tip[stage_num] = -B_yu[stage_num].completeOrthogonalDecomposition().pseudoInverse()*B_w_tip[stage_num];

    J_w_tip[stage_num] = E_yu * yu_dot_w_tip[stage_num];

    w_tip_dot[stage_num] = (J_w_tip[stage_num].transpose()*J_w_tip[stage_num] + 0.1*I_6x6).inverse()*J_w_tip[stage_num].transpose();

    conv.block<1, 3>(0, 0) = -eta.segment<3>(1).transpose();
    conv.block<3, 3>(1, 0) = eta(0) * I_3x3 + MathUtils::skew_m(eta.segment<3>(1));
    eta_dot = conv*J_w_tip[stage_num].block<3, 6>(3, 0);

    J_w_tip_eta[stage_num] << J_w_tip[stage_num].block<3, 6>(0,0), eta_dot;


}