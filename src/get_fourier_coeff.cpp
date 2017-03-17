#include <get_fourier_coeff.h>
#include <math.h>
#include <iostream>

FourierReturn *get_fourier_coeff(MatrixXd &mu, MatrixXd &X, MatrixXd &Y, int Nkx, int Nky, double Lx, double Ly) {
	// Initializing Fourier Coefficients of prior
	MatrixXd muk = MatrixXd::Zero(Nkx, Nky);

	VectorXd temp1(Nky);
	temp1 << 1, sqrt(0.5) * VectorXd::Ones(Nky-1);

	RowVectorXd temp2(Nkx);
	temp2 << 1, sqrt(0.5) * RowVectorXd::Ones(Nkx-1);

	MatrixXd HK = sqrt(Lx*Ly) * temp1 * temp2;

	double const1 = Nkx * M_PI / Lx;
	double const2 = Nky * M_PI / Ly;

	// Our implementation
	for (int i = 0; i < Nkx; i++) {
		for (int j = 0; j < Nky; j++) {
			/*std::cout << "X size is : " << X.rows() << " x " << X.cols() << std::endl;
			std::cout << "Y size is : " << Y.rows() << " y " << Y.cols() << std::endl;*/
			MatrixXd temp1 = ((const1 * X).array().cos()).matrix();
			MatrixXd temp2 = ((const2 * Y).array().cos()).matrix();

			/*std::cout << "mu size is : " << mu.rows() << " x " << mu.cols() << std::endl;
			std::cout << "temp1 size is : " << temp1.rows() << " x " << temp1.cols() << std::endl;
			std::cout << "temp2 size is : " << temp2.rows() << " x " << temp2.cols() << std::endl;*/
			muk(i, j) = ((mu.cwiseQuotient(temp1)).cwiseQuotient(temp2)).sum() / HK(j, i);
		}
	}

	FourierReturn *result = (FourierReturn*)calloc(1, sizeof(FourierReturn));
	result->muk = muk;
	result->HK = HK;

	return result;
}