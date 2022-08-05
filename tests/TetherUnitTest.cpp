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

    /* 
    Here, we set some inputs to test the Jacobians solved by TetherUnit_Solver. 
     */
    Eigen::Matrix<double, 7, 1> poseDesired;
    poseDesired << 2.8, 0.0, 0.4, 1, 0, 0, 0 ;
    double mass_distribution = 0.035; 
    double tether_length = 3.1;
    double g = 9.81;
    double w_t = 2500; 
    double Kp = 0.50; 
    double lambdaDLS = 5; 
    double lambdaLyap = 50.0; 
    unsigned int numIntegrationSteps = 50; 

    Eigen::MatrixXd proximalStates(17, 1);
    proximalStates << 0, 0, 0, 1, 0, 0, 0, mass_distribution * tether_length * g, 0, 0, 0, 0.21544033, 0, 0, 0, 0.05, 0;

    TetherUnit_Solver TSolver(&i1, &i2, mass_distribution, tether_length, numIntegrationSteps, 
            lambdaLyap, lambdaDLS, w_t, Kp, proximalStates); 

    std::cout.precision(10);

    std::cout << "Boundary conditions" << TSolver.getBoundaryConditions() << "\n\n";

    TSolver.timer.tic();

    // while (TSolver.getPoseError().norm() > 5e-3)

    Eigen::Matrix<double, 6, 1> tipWrench; 
    tipWrench.setZero(); 
    tipWrench(0, 0) = -0.035; 
    tipWrench(2, 0) = -0.0005; 
    tipWrench(4, 0) = 0.0015;

    
    while (TSolver.getPoseError().norm() > 1e-2)
    
    {

        TSolver.timer.tic();
        // TSolver.solveReactionForcesStep(poseDesired);
        TSolver.simulateStep(tipWrench);
        // std::cout << "Jac: \n" << TSolver.getJacobianEta_wrt_tip() << "\n\n";
        std::cout << "poseError norm: " << TSolver.getPoseError().norm() << "\n\n";
        std::cout << "Proximal states: " << TSolver.getProximalStates().transpose() << "\n\n";
        std::cout << "Distal states: " << TSolver.getDistalStates().transpose() << "\n\n"; 
        // std::cout << "Tip Wrench: " << TSolver.getTipWrench() << "\n\n";
        if (std::isnan(TSolver.getPoseError().norm()))
        {

            std::cout << "SINGULARITY REACHED. " << "\n\n";
            break; 

        }
        std::cout << "Boundary conditions norm: " << TSolver.getBoundaryConditions().norm() << "\n\n";
        TSolver.timer.toc();

    }

    // tipWrench.setZero(); 
    // // tipWrench(0, 0) = -0.030; 
    // // tipWrench(2, 0) = 0.001; 
    // tipWrench(4, 0) = 0.00835;


    // for (int i = 0; i < 2000; ++i)
    
    // {

    //     TSolver.timer.tic();
    //     // TSolver.solveReactionForcesStep(poseDesired);
    //     TSolver.simulateStep(tipWrench);
    //     // std::cout << "Jac: \n" << TSolver.getJacobianEta_wrt_tip() << "\n\n";
    //     std::cout << "poseError norm: " << TSolver.getPoseError().norm() << "\n\n";
    //     std::cout << "Proximal states: " << TSolver.getProximalStates().transpose() << "\n\n";
    //     std::cout << "Distal states: " << TSolver.getDistalStates().transpose() << "\n\n"; 
    //     std::cout << "Tip Wrench: " << TSolver.getTipWrench() << "\n\n";
    //     if (std::isnan(TSolver.getPoseError().norm()))
    //     {

    //         std::cout << "SINGULARITY REACHED. " << "\n\n";
    //         break; 

    //     }
    //     std::cout << "Boundary conditions norm: " << TSolver.getBoundaryConditions().norm() << "\n\n";
    //     TSolver.timer.toc();

    // }

    TSolver.integrateFullStates();
    TSolver.getFullStates("test.txt");

    std::cout << "Proximal states: " << TSolver.getProximalStates() << "\n\n";
    std::cout << "Distal states: " << TSolver.getDistalStates() << "\n\n"; 
    std::cout << "Boundary conditions" << TSolver.getBoundaryConditions() << "\n\n";
    std::cout << "Tip Wrench: " << TSolver.getTipWrench() << "\n\n";


    free(capsule);
    free(capsule_step);

    return 0; 

}