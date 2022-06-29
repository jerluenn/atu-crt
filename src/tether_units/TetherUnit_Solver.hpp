#ifndef TETHERUNIT_H
#define TETHERUNIT_H 

#include <iostream>
#include <stdlib.h>
#include <fstream>
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

#define assertm(exp, msg) assert(((void)msg, exp))

class TetherUnit_Solver 

{

    public: 
        TetherUnit_Solver(IntegrationInterface* integrator_, IntegrationInterface* integratorStep_, double mass_distribution_, double tether_length_, unsigned int num_integrationSteps_, double lambdaLyapunov_, double lambdaDLS_, double Kp_, Eigen::MatrixXd initialSolution);
        TetherUnit_Solver(); 
        virtual ~TetherUnit_Solver(); 
        Eigen::MatrixXd getProximalStates();
        Eigen::MatrixXd getFullStates(std::string fileName);
        Eigen::MatrixXd getFullStates();
        Eigen::MatrixXd getDistalStates();
        unsigned int getNumIntegrationSteps();
        Eigen::Matrix<double, 7, 1> getDistalPose();
        Eigen::Matrix<double, 7, 6> getJacobianEta_wrt_tip(); 
        void integrateDistalStates();
        void integrateFullStates();
        void setInitialConditions(Eigen::MatrixXd initialConditions); 
        void setInitialConditions(Eigen::Matrix<double, 6, 1> initialConditions);
        Eigen::Matrix<double, 6, 1> getTipWrench();

        Eigen::Matrix<double, 6, 1> getBoundaryConditions();
        Eigen::Matrix<double, 6, 1> getBoundaryConditions_IncreaseWrenchTip(unsigned int index);
        Eigen::Matrix<double, 6, 1> getBoundaryConditions(Eigen::MatrixXd distalConditions_);
        void setTipWrench(Eigen::Matrix<double, 6, 1> p_FM);
        void updateTipWrench(Eigen::Matrix<double, 6, 1> twist);
        void simulateStep(Eigen::Matrix<double, 6, 1> tip_wrench);
        void setTau(double tau); 
        void solveJacobians(); 
        void solveReactionForcesStep(Eigen::MatrixXd poseDesired);
        Eigen::Matrix<double, 7, 1> getPoseError();

        MathUtils::Timer timer;

    private: 

        Eigen::MatrixXd fullStates; 
        Eigen::MatrixXd E_yu, B_yu, E_w_tip, B_w_tip, J_w_tip, yu_dot_w_tip, eta_dot, w_tip_dot, J_w_tip_eta; 
        Eigen::MatrixXd distalStates;
        Eigen::MatrixXd proximalStates;
        Eigen::MatrixXd dummyStates;
        Eigen::Matrix<double, 7, 1> poseError;
        Eigen::Matrix<double, 7, 7> I_7x7; 
        Eigen::MatrixXd J_test;
        void saveData(std::string fileName, Eigen::MatrixXd matrix);
        Eigen::MatrixXd integrateWithIncrement(unsigned int index);
        Eigen::Matrix<double, 6, 1> tipWrench, tipWrenchInput; // external tip wrench.
        double mass_distribution; 
        double tether_length;
        const unsigned int numStates = 17; 
        unsigned int num_p, num_eta, num_int_wrench, num_alpha, num_curvature, num_kappa, num_tau;
        Eigen::MatrixXd R_dynamic; 
        unsigned int num_integrationSteps;
        constexpr static double g = 9.81; 
        constexpr static double EPS = 1e-10;
        IntegrationInterface* integrator;
        IntegrationInterface* integratorStep;
        double lambdaLyapunov;
        double lambdaDLS;
        double Kp;
        double dt = 0.01;

};

#endif