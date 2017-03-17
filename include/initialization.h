/*
 * File Header:
     Initializ the coverage problem. Domain, speed limits, number of agents...
 */

#ifndef INITIALIZATION
#define INITIALIZATION

#include <Eigen/Dense>
#include <math.h>

using namespace Eigen;

// Agent locations
struct Pose {
    Vector3d x;
    Vector3d y;
    VectorXd theta;
};

struct DomainBounds {
    double xmin;
    double xmax;
    double ymin;
    double ymax;
};

struct ERG {
    double s; // parameter of soblev norm
    int Nkx;
    int Nky;
    // used in calculation CK and Bj in SMC update
    MatrixXd KX;
    MatrixXd KY;
    MatrixXd LK;
    MatrixXd HK;  // normalizer of fourier basis functions, which will be assigned in GetFourier Coeff and
          // used in multiple places
    MatrixXd muk;
};

struct OPT {
	Vector2d L;    // [Lx; Ly]
	int nagents;  //# of agents 
	VectorXd vlb;  // forward velocity lower bound
	VectorXd vub;  // forward velocity upper bound
	VectorXd wlb;  // angular velocity lower bound
	VectorXd wub;  // angular velocity upper bound
};

#endif