#include <TetherUnit_Solver.hpp>
#include <IntegratorInterface.hpp>
#include "acados_sim_solver_tetherunit_integrator.h"
#include "acados_sim_solver_tetherunit_stepIntegrator.h"

int main () 

{

    sim_solver_capsule *capsule = tetherunit_integrator_acados_sim_solver_create_capsule(); 
    tetherunit_integrator_acados_sim_create(capsule); 

    sim_solver_capsule *capsule_step = tetherunit_stepIntegrator_acados_sim_solver_create_capsule(); 
    tetherunit_stepIntegrator_acados_sim_create(capsule_step); 

    IntegrationInterface i1(capsule), i2(capsule_step);

    TetherUnit_Solver tetherobject(&i1, &i2, 0.035, 2.3, 50); 

    Eigen::Matrix<double, 6, 1> tipWrench; 
    Eigen::Matrix<double, 6, 1> twist;

    // twist << -0.1, 0.0, 0.1, 0.0, -0.8, 0.0;

    // tipWrench << 0.02, 0.0, 0.0, 0.0, -0.01, 0.0;
    tipWrench << -0.05, 0.008, 0.15, 0.0, -0.009, 0.0; 
    

    

    std::cout.precision(10);

    std::cout << "Distal states: " << tetherobject.integrateDistalStates() << "\n\n"; 

    std::cout << "Boundary conditions" << tetherobject.getBoundaryConditions() << "\n\n";

    tetherobject.timer.tic();

    for (int i = 0; i < 200; ++i) 
    
    {

        // tetherobject.timer.tic();

        // tetherobject.updateTipWrench(twist);
        // std::cout << "Boundary conditions" << tetherobject.getBoundaryConditions() << "\n\n";
        tetherobject.simulateStep(tipWrench);

        // tetherobject.timer.toc();
        

    }

    std::cout << "Boundary conditions" << tetherobject.getBoundaryConditions() << "\n\n";


    for (int i = 0; i < 200; ++i) 
    
    {

        // tetherobject.timer.tic();

        // tetherobject.updateTipWrench(twist);
        // std::cout << "Boundary conditions" << tetherobject.getBoundaryConditions() << "\n\n";
        tetherobject.simulateStep(-tipWrench);

        // tetherobject.timer.toc();
        

    }

    std::cout << "Boundary conditions" << tetherobject.getBoundaryConditions() << "\n\n";


    for (int i = 0; i < 200; ++i) 
    
    {

        // tetherobject.timer.tic();

        // tetherobject.updateTipWrench(twist);
        // std::cout << "Boundary conditions" << tetherobject.getBoundaryConditions() << "\n\n";
        tetherobject.simulateStep(tipWrench);

        // tetherobject.timer.toc();
        

    }

    tetherobject.timer.toc();
    

    std::cout << "Proximal states: " << tetherobject.getProximalStates() << "\n\n";

    std::cout << "Distal states: " << tetherobject.integrateDistalStates() << "\n\n"; 

    std::cout << "Boundary conditions" << tetherobject.getBoundaryConditions() << "\n\n";

    std::cout << "Tip Wrench: " << tetherobject.getTipWrench() << "\n\n";
    

    // IpoptSolver ipopt;

    // ipopt.solve();

    return 0; 

}