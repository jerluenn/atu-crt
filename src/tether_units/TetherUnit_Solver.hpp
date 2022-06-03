#ifndef TETHERUNIT_H
#define TETHERUNIT_H 

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <IntegratorInterface.hpp>
#include <MathUtils.hpp>
#include <vector>
#include <cassert>
#include <chrono>
#include <map>
#include <ifopt/ipopt_solver.h>
#include "IPOPTTetherUnit.hpp"

#define assertm(exp, msg) assert(((void)msg, exp))

class TetherUnit_Solver 

{

    public: 

        TetherUnit_Solver(IntegrationInterface* integrator_, IntegrationInterface* integratorsStep_, double mass_distribution_, double tether_length_, unsigned int num_integrationSteps_);
        TetherUnit_Solver(); 
        virtual ~TetherUnit_Solver(); 
        Eigen::MatrixXd getProximalStates();
        Eigen::MatrixXd getFullStates();
        void setInitialConditions(Eigen::MatrixXd initialConditions); 
        Eigen::MatrixXd integrateDistalStates(); 
        Eigen::MatrixXd integrateFullStates();
        void setTau(double tau); 
        void SolveJacobians(); 

    private: 

        Eigen::MatrixXd fullStates; 
        Eigen::MatrixXd Jacobians_ChangeName; 
        Eigen::MatrixXd distalStates;
        Eigen::MatrixXd proximalStates;
        double mass_distribution; 
        double tether_length;
        const unsigned int numStates = 17; 
        unsigned int num_integrationSteps;
        constexpr static double g = 9.81; 
        IntegrationInterface* integrator;
        IntegrationInterface* integratorStep;

};

#endif