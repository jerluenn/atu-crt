#ifndef TETHERUNIT_INTERFACE_H
#define TETHERUNIT_INTERFACE_H 

#include <TetherUnit_Solver.hpp>
#include "IPOPTTetherUnit.hpp"

class TetherUnit_Interface

{

    public: 

        TetherUnit_Interface(TetherUnit_Solver* TetherUnit_); 
        void solveBoundaryValueProblem(bool print_level);
        virtual ~TetherUnit_Interface();

    private: 

        TetherUnit_Solver* TU_Solver; 

};

#endif