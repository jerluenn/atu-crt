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

    // tipWrench << 0.02, 0.0, 0.0, 0.0, -0.01, 0.0;
    tipWrench << 0.01, 0.0, 0.0, 0.0, -0.005, 0.0; 

    tetherobject.timer.tic();

    std::cout.precision(10);

    std::cout << "Distal states: " << tetherobject.integrateDistalStates() << "\n\n"; 

    std::cout << "Boundary conditions" << tetherobject.getBoundaryConditions() << "\n\n";

    for (int i = 0; i < 2200; ++i) 
    
    {

        tetherobject.simulateStep(tipWrench);

    }

    tetherobject.timer.toc();

    std::cout << "Proximal states: " << tetherobject.getProximalStates() << "\n\n";

    std::cout << "Distal states: " << tetherobject.integrateDistalStates() << "\n\n"; 

    std::cout << "Boundary conditions" << tetherobject.getBoundaryConditions() << "\n\n";
    

    // IpoptSolver ipopt;

    // ipopt.solve();

    return 0; 

}