#include <calculate_ergodicity.h>

double calculate_ergodicity(MatrixXd &ck, MatrixXd &LK, MatrixXd &muk) {
	MatrixXd temp = LK.cwiseProduct(ck-muk);
	return (temp.array().square()).sum();
}