#include <TetherUnit_Solver.hpp>
#include <IntegratorInterface.hpp>
#include <TetherUnit_Interface.hpp>
#include "acados_sim_solver_tetherunit_integrator.h"
#include "acados_sim_solver_tetherunit_stepIntegrator.h"

int main () 

{

    /* 
    Here, we get the integration solvers from the acados include directories. 
    */


    sim_solver_capsule *capsule = tetherunit_integrator_acados_sim_solver_create_capsule(); 
    tetherunit_integrator_acados_sim_create(capsule); 

    sim_solver_capsule *capsule_step = tetherunit_stepIntegrator_acados_sim_solver_create_capsule(); 
    tetherunit_stepIntegrator_acados_sim_create(capsule_step); 

    IntegrationInterface i1(capsule), i2(capsule_step);

    TetherUnit_Solver TSolver(&i1, &i2, 0.035, 2.3, 50); 

    /* 
    Here, we set some inputs to test the Jacobians solved by TetherUnit_Solver. 
     */

    Eigen::Matrix<double, 6, 1> tipWrench; 
    Eigen::Matrix<double, 6, 1> twist;
    Eigen::Matrix<double, 7, 1> poseDesired, poseCurrent, poseError; 
    Eigen::Matrix<double, 7, 6> JacobianEta; 
    Eigen::Matrix<double, 7, 7> I_7x7; 
    I_7x7.setIdentity();
    double lambda = 1;
    // twist << -0.1, 0.0, 0.1, 0.0, -0.8, 0.0;
    tipWrench << -0.5, 0.0, 0.8, 0.0, -0.085, 0.0; 
    poseDesired << -1.4, 0, 1.6, 1, 0, 0, 0 ;
    TSolver.integrateDistalStates(); 
    poseCurrent = TSolver.getDistalPose();
    poseError = poseDesired - poseCurrent;
    std::cout << poseError.norm() << "\n";
    TSolver.solveJacobians(); 
    JacobianEta = TSolver.getJacobianEta_wrt_tip();
    

    std::cout.precision(10);

    std::cout << "Distal states: " << TSolver.getDistalStates() << "\n\n"; 

    std::cout << "Boundary conditions" << TSolver.getBoundaryConditions() << "\n\n";

    TSolver.timer.tic();

    while (poseError.norm() > 5e-3)
    
    {

        // TSolver.updateTipWrench(twist);
        // std::cout << "Boundary conditions" << TSolver.getBoundaryConditions() << "\n\n";
        tipWrench = 5*JacobianEta.transpose() * (JacobianEta * JacobianEta.transpose() + pow(lambda, 2)*I_7x7).inverse() * poseError;
        TSolver.simulateStep(tipWrench);
        poseCurrent = TSolver.getDistalPose();
        poseError = poseDesired - poseCurrent;
        std::cout << "poseError norm: " << poseError.norm() << "\n\n";
        std::cout << "Boundary conditions norm: " << TSolver.getBoundaryConditions().norm() << "\n\n";

    }

    TSolver.integrateFullStates();
    TSolver.getFullStates("test.txt");


    TSolver.timer.toc();
    

    std::cout << "Proximal states: " << TSolver.getProximalStates() << "\n\n";

    std::cout << "Distal states: " << TSolver.getDistalStates() << "\n\n"; 

    std::cout << "Boundary conditions" << TSolver.getBoundaryConditions() << "\n\n";

    std::cout << "Tip Wrench: " << TSolver.getTipWrench() << "\n\n";
    
    TetherUnit_Interface T_Object(&TSolver);
    T_Object.solveBoundaryValueProblem(true);

    return 0; 

}