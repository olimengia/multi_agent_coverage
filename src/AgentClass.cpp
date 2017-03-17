/*
    AgentClass: a class for the agent (robot) object
*/

#include <AgentClass.h>

// Return pose 
Vector3d Agent::get_pose(void) {
    Vector3d pose;
    pose(0) = x;
    pose(1) = y;
    pose(2) = theta;
	return pose;
}

// Return forward velocity bounds. 1st element is velocity lower bound
// 2nd element is velocity upper bound
Vector2d Agent::get_forward_velocity_bound() {
    Vector2d m;
    m(0) = vlb;
    m(1) = vub;
    return m;
}

// Return angular velocity bounds 1st element is angular velocity lower bound
// 2nd element is angular velocity upper bound
Vector2d Agent::get_angular_velocity_bound() {
    Vector2d m;
    m(0) = wlb;
    m(1) = wub;
    return m;
}

// Set x pose
void Agent::set_x(double x_new) {
    x = x_new;
}

// Set y pose
void Agent::set_y(double y_new) {
    y = y_new;
}

// Set theta pose
void Agent::set_theta(double theta_new) {
    theta = theta_new;
}

void Agent::set_v_lower(double vlb_new) {
    vlb = vlb_new;
}

void Agent::set_v_upper(double vub_new) {
    vub = vub_new;
}

void Agent::set_w_lower(double wlb_new) {
    wlb = wlb_new;
}

void Agent::set_w_upper(double wub_new) {
    wub = wub_new;
}