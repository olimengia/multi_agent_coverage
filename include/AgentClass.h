/*
    Header file for AgentClass.cpp
*/
#ifndef AGENT_CLASS_H
#define AGENT_CLASS_H

#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

class Agent {
    // Pose
    double x;
    double y;
    double theta;
    // Bounds
    double vlb;  // forward velocity lower bound
    double vub;  // forward velocity upper bound
    double wlb;  // angular velocity lower bound
    double wub;  // angular velocity upper bound
  
  public:
    Vector3d get_pose(void);
    Vector2d get_forward_velocity_bound(void);
    Vector2d get_angular_velocity_bound(void);
    void set_x(double x_new);
    void set_y(double y_new);
    void set_theta(double theta_new);
    void set_v_lower(double vlb_new);
    void set_v_upper(double vub_new);
    void set_w_lower(double wlb_new);
    void set_w_upper(double wub_new);
};
#endif