#ifndef SMC_UPDATE
#define SMC_UPDATE

#include <Eigen/Dense>
#include <initialization.h>

using namespace Eigen;

struct SMCUpdateReturn {
	VectorXd x;
	VectorXd y;
	VectorXd theta;
};

SMCUpdateReturn *smc_update(Pose &pose, OPT &opt, DomainBounds &domain_bounds, ERG &erg, MatrixXd &Ck, double time, double dt);

#endif