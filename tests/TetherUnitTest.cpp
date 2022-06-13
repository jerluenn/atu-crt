#include <TetherUnit_Solver.hpp>
#include <IntegratorInterface.hpp>
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

    TetherUnit_Solver tetherobject(&i1, &i2, 0.035, 2.3, 50); 

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
    poseDesired << -1.6, -0.4, 1.2, 0.9914, 0.1305, 0, 0 ;
    tetherobject.integrateDistalStates(); 
    poseCurrent = tetherobject.getDistalPose();
    poseError = poseDesired - poseCurrent;
    std::cout << poseError.norm() << "\n";
    tetherobject.solveJacobians(); 
    JacobianEta = tetherobject.getJacobianEta_wrt_tip();
    

    std::cout.precision(10);

    std::cout << "Distal states: " << tetherobject.getDistalStates() << "\n\n"; 

    std::cout << "Boundary conditions" << tetherobject.getBoundaryConditions() << "\n\n";

    tetherobject.timer.tic();

    while (poseError.norm() > 5e-3)
    
    {

        // tetherobject.updateTipWrench(twist);
        // std::cout << "Boundary conditions" << tetherobject.getBoundaryConditions() << "\n\n";
        tipWrench = 5*JacobianEta.transpose() * (JacobianEta * JacobianEta.transpose() + pow(lambda, 2)*I_7x7).inverse() * poseError;
        tetherobject.simulateStep(tipWrench);
        poseCurrent = tetherobject.getDistalPose();
        poseError = poseDesired - poseCurrent;
        std::cout << "poseError norm: " << poseError.norm() << "\n\n";
        std::cout << "Boundary conditions norm: " << tetherobject.getBoundaryConditions().norm() << "\n\n";

    }

    tetherobject.integrateFullStates();
    tetherobject.getFullStates("test.txt");


    tetherobject.timer.toc();
    

    std::cout << "Proximal states: " << tetherobject.getProximalStates() << "\n\n";

    std::cout << "Distal states: " << tetherobject.getDistalStates() << "\n\n"; 

    std::cout << "Boundary conditions" << tetherobject.getBoundaryConditions() << "\n\n";

    std::cout << "Tip Wrench: " << tetherobject.getTipWrench() << "\n\n";
    

    // IpoptSolver ipopt;

    // ipopt.solve();

    return 0; 

}