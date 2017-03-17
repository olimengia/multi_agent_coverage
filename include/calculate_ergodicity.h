#ifndef CALCULATE_ERGODICITY
#define CALCULATE_ERGODICITY

#include <Eigen/Dense>

using namespace Eigen;

double calculate_ergodicity(MatrixXd &ck, MatrixXd &LK, MatrixXd &muk);

#endif