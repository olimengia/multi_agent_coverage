#ifndef GET_FOURIER_COEFF
#define GET_FOURIER_COEFF

#include <Eigen/Dense>

using namespace Eigen;

struct FourierReturn {
	MatrixXd muk;
	MatrixXd HK;
};

FourierReturn *get_fourier_coeff(MatrixXd &mu, MatrixXd &X, MatrixXd &Y, int Nkx, int Nky, double Lx, double Ly);

#endif