#ifndef GENERATE_UTILITY_MAP
#define GENERATE_UTILITY_MAP

#include <Eigen/Dense>

using namespace Eigen;

MatrixXd mvnpdf(int rows, double x_interval, int cols, double y_interval, Vector2d &mean, Matrix2d &covar);

MatrixXd generate_utility_map(double xmin, double xmax, double xdel, double ymin, double ymax, double ydel, int addnoise);

#endif