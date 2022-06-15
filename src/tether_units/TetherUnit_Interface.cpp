#include "TetherUnit_Interface.hpp"

using namespace ifopt;

TetherUnit_Interface::TetherUnit_Interface(TetherUnit_Solver* TetherUnit_) 

{

    TU_Solver = TetherUnit_; 

}

void TetherUnit_Interface::solveBoundaryValueProblem(bool print_level) 

{

    Problem nlp; 

    nlp.AddVariableSet  (std::make_shared<ExVariables>());
    nlp.AddConstraintSet(std::make_shared<ExConstraint>("Constraints", TU_Solver));
    nlp.AddCostSet      (std::make_shared<ExCost>("Cost", TU_Solver));
    // nlp.PrintCurrent();

    IpoptSolver ipopt;
    ipopt.SetOption("linear_solver", "mumps");
    ipopt.SetOption("jacobian_approximation", "exact");

    ipopt.Solve(nlp);
    // Eigen::VectorXd x = nlp.GetOptVariables()->GetValues();
    // std::cout << x.transpose() << std::endl;

    if (print_level == true) 
    
    {

        std::cout << "Sol: " << "TestString" << "\n";

    }

}

TetherUnit_Interface::~TetherUnit_Interface() 

{


};