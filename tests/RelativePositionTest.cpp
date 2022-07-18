#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "ros/ros.h"
#include "geometry_msgs/PoseStamped.h"

typedef struct Pose

{

    Eigen::Quaterniond eta; 
    Eigen::Vector3d p; 

}  Pose;

class SubscriberNode 

{

    public: 

        SubscriberNode(ros::NodeHandle *nh) 
        {

            ROS_INFO("Opening node...");
            sub_p1 = nh->subscribe("vrpn_client_node/ground_station/pose", 1, &SubscriberNode::p1_callback, this);
            sub_p2 = nh->subscribe("vrpn_client_node/tension_mount/pose", 1, &SubscriberNode::p2_callback, this);    
            pub_p = nh->advertise<geometry_msgs::PoseStamped>("/relative_pose", 10);        
            initialise();
            ros::Rate rate(100); 
            
            while (ros::ok()) 
            {

                pub_relativePose();
                rate.sleep();
                ros::spinOnce();

            }
            

        }

        void initialise() 
        {


        }

        void p1_callback(const geometry_msgs::PoseStamped::ConstPtr& msg )
        {

            p1.p << msg->pose.position.x, msg->pose.position.y, msg->pose.position.z;
            p1.eta.w() = msg->pose.orientation.w;
            p1.eta.vec() << msg->pose.orientation.x, msg->pose.orientation.y, msg->pose.orientation.z;

        }

        void p2_callback(const geometry_msgs::PoseStamped::ConstPtr& msg )
        {

            p2.p << msg->pose.position.x, msg->pose.position.y, msg->pose.position.z;
            p2.eta.w() = msg->pose.orientation.w;
            p2.eta.vec() << msg->pose.orientation.x, msg->pose.orientation.y, msg->pose.orientation.z;

        }


        void pub_relativePose() 
        {

            Eigen::Quaterniond q_tmp(0, 0, 0, 0);
            Eigen::Quaterniond relativeRotation;
            p1.eta.normalize();
            p2.eta.normalize();
            relativeRotation = p1.eta.inverse() * p2.eta;
            q_tmp.vec() = p2.p - p1.p; 
            q_tmp = p1.eta.inverse() * q_tmp * p1.eta;
            relativePose.pose.position.x = q_tmp.vec()(0);
            relativePose.pose.position.y = q_tmp.vec()(1);
            relativePose.pose.position.z = q_tmp.vec()(2);
            relativePose.pose.orientation.w = relativeRotation.w();
            relativePose.pose.orientation.x = relativeRotation.vec()(0);
            relativePose.pose.orientation.y = relativeRotation.vec()(1);
            relativePose.pose.orientation.z = relativeRotation.vec()(2);
            pub_p.publish(relativePose);    
       

        }

    private: 

        ros::Subscriber sub_p1; 
        ros::Subscriber sub_p2; 
        ros::Publisher pub_p; 
        geometry_msgs::PoseStamped relativePose; // pose {1,2}
        Pose p1; 
        Pose p2; 
        Pose relativePoseEig;


};

int main (int argc, char** argv) 

{

    ros::init(argc, argv, "relative_pose_test");
    ros::NodeHandle n; 
    SubscriberNode s(&n);

    return 0; 

}