#include <generate_utility_map.h>
#include <fstream>
#include <iostream>
#include <math.h>


MatrixXd mvnpdf(int cols, double x_interval, int rows, double y_interval, Vector2d &mean, Matrix2d &covar) {
	MatrixXd result(rows, cols);

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			Vector2d loc;
			//std::cout << "mean: " << mean << std::endl;
			loc << i * y_interval, j * x_interval;
			//std::cout << "loc: " << loc << std::endl;
			//std::cout << "covar: " << covar << std::endl;

			result(i, j) = exp(-0.5 * (loc - mean).transpose() * covar.inverse() * (loc - mean)) /
							   sqrt(pow(2*M_PI, mean.size()) * covar.determinant());
			
			/*std::cout << "first part: " << exp(-0.5 * (loc - mean).transpose() * covar.inverse() * (loc - mean)) << std::endl;
			std::cout << "power: " << pow(2*M_PI, mean.size()) << std::endl;
			std::cout << "determinant: " << covar.determinant() << std::endl;
			std::cout << "second part: " << sqrt(pow(2*M_PI, mean.size()) * covar.determinant()) << std::endl;*/
			// Set anything negative to 0
			if (result(i, j) < 0)
				result(i, j) = 0;
			
		}
	}
	return result;
}

MatrixXd generate_utility_map(double xmin, double xmax, double xdel, double ymin, double ymax, double ydel, int addnoise) {
	int x_range = (int)((xmax-xdel-xmin)/xdel + 1);
	int y_range = (int)((ymax-ydel-ymin)/ydel + 1);

	double n;
	if (addnoise)
		n = .05;
	else
		n = 0;

	Vector2d m;
	m << 100, 200;
	Matrix2d s = 150 * (Matrix2d::Identity(2,2));
	Vector2d m2;
	m2 << 220, 200;
	Matrix2d s2 = 800 * (Matrix2d::Identity(2,2));
	Vector2d m3;
	m3 << 120, 120;
	Matrix2d s3 = 600 * (Matrix2d::Identity(2,2));

	MatrixXd G1 = mvnpdf(x_range, xdel, y_range, ydel, m, s);
	MatrixXd G2 = mvnpdf(x_range, xdel, y_range, ydel, m2, s2);
	MatrixXd G3 = mvnpdf(x_range, xdel, y_range, ydel, m3, s3);
	MatrixXd G = G1+G2+G3;

	// Normalization 
	G = G / G.maxCoeff();

	return G;
}