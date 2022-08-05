#include <TetherUnit_Solver.hpp>
#include <IntegratorInterface.hpp>
#include "acados_sim_solver_tetherunit_integrator.h"
#include "acados_sim_solver_tetherunit_stepIntegrator.h"
#include "ros/ros.h"
#include "std_msgs/String.h"
#include "std_msgs/Float64MultiArray.h"
#include "geometry_msgs/PoseStamped.h"

class Single_ATU 

{

    public: 

        Single_ATU(ros::NodeHandle* nh, TetherUnit_Solver* TU_Solver)
        {

            ROS_INFO("Opening node...");
            this->TU_Solver = TU_Solver;
            sub_relative_pose = nh->subscribe("/relative_pose", 1, &Single_ATU::relative_pose_callback, this);
            sub_loadcell_filtered = nh->subscribe("/loadcell_filtered", 1, &Single_ATU::loadcell_callback, this);
            std::cout << "Boundary conditions" << TU_Solver->getBoundaryConditions() << "\n\n";
            poseDesired << -2.8, 0.0, 1.2, 1.0, 0.0, 0.0, 0.0;

            ros::Rate rate(500); 
            
            while (ros::ok()) 
            {

                TU_Solver->solveReactionForcesStep(poseDesired);
                std::cout << "Distal Pose: " << TU_Solver->getDistalStates() << "\n\n";
                std::cout << "Desired POse: " << poseDesired << "\n\n";
                std::cout << "Boundary conditions norm: " << TU_Solver->getBoundaryConditions().norm() << "\n\n";
                rate.sleep();
                ros::spinOnce();

            }
            

        }

        void relative_pose_callback(const geometry_msgs::PoseStamped::ConstPtr& msg)
        {
        
            // poseDesired << msg->pose.position.x, msg->pose.position.y, msg->pose.position.z, msg->pose.orientation.w, msg->pose.orientation.x,
            // msg->pose.orientation.y, msg->pose.orientation.z;  

        }

        void loadcell_callback(const std_msgs::Float64MultiArray& msg) 
        {

            TU_Solver->setTau(msg.data.data()[2]);

        }

        

    private: 

        ros::Subscriber sub_relative_pose;
        ros::Subscriber sub_loadcell_filtered;
        TetherUnit_Solver* TU_Solver; 
        Eigen::Matrix<double, 7, 1> poseDesired;

};

int main (int argc, char **argv) 

{

    ros::init(argc, argv, "single_atu");
    ros::NodeHandle n; 


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
    poseDesired << -1.2, 0.5, 0.9, 1, 0, 0, 0 ;
    double mass_distribution = 0.035; 
    double tether_length = 3.1;
    double g = 9.81;
    double w_t = 0.5; 
    double Kp = 0.25; 
    double lambdaDLS = 6.0; 
    double lambdaLyap = 50.0; 
    unsigned int numIntegrationSteps = 50; 

    Eigen::MatrixXd proximalStates(17, 1);
    proximalStates << 0, 0, 0, 1, 0, 0, 0, mass_distribution * tether_length * g, 0, 0, 0, 0.21544033, 0, 0, 0, 0.05, 0;

    TetherUnit_Solver TSolver(&i1, &i2, mass_distribution, tether_length, numIntegrationSteps, 
            lambdaLyap, lambdaDLS, w_t, Kp, proximalStates); 


    Single_ATU(&n, &TSolver);

    return 0; 

}