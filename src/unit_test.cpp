#include <initialization.h>
#include <AgentClass.h>
#include <calculate_ergodicity.h>
#include <generate_utility_map.h>
#include <get_fourier_coeff.h>
#include <smc_update.h>
#include <math.h>

using namespace Eigen;

// Domain bounds
#define DOMAIN_XMIN 0.0
#define DOMAIN_XMAX 300.0
#define DOMAIN_YMIN 0.0
#define DOMAIN_YMAX 300.0
#define NUM_AGENTS 3
#define NSTEPS 5000
#define DT 0.1
#define PI 3.14159265358979323846

int main(void){
	double Lx = DOMAIN_XMAX - DOMAIN_XMIN;
	double Ly = DOMAIN_YMAX - DOMAIN_YMIN;

	double vlb_init = 0.1;
	double vub_init = 5;
	double wlb_init = -0.5;
	double wub_init = 0.5;

	Vector3d x;
	x(0) = 100;
	x(1) = 200;
	x(2) = 50;
	//x << 100 << 200 << 50;
	Vector3d y;
	y(0) = 50;
	y(1) = 75;
	y(2) = 250;
	//y << 50 << 75 << 250;
	Vector3d thetas;
	thetas(0) = PI;
	thetas(1) = PI;
	thetas(2) = PI;
	//thetas << M_PI << M_PI << M_PI;

	// Agents params: Velocity Bounds, number of agents...

	/********* Initialization ends *************/
	// Initialize the agent locations and limits
	Agent agents[NUM_AGENTS];

	for (int i = 0; i < NUM_AGENTS; i++) {
		agents[i].set_x(x(i));
		agents[i].set_y(y(i));
		agents[i].set_theta(thetas(i));
		agents[i].set_v_lower(vlb_init);
		agents[i].set_v_upper(vub_init);
		agents[i].set_w_lower(wlb_init);
		agents[i].set_w_upper(wub_init);
	}

	// Initialize ergodicity params
	double s = 1.5;
	int Nkx = 50;
	int Nky = 50;
	MatrixXi KX = (VectorXi::LinSpaced(Nkx, 0, Nkx - 1)) * RowVectorXi::Ones(Nky);   // Nkx x 1
	MatrixXi KY = VectorXi::Ones(Nky) * RowVectorXi::LinSpaced(Nky, 0, Nky - 1);  // Nky x 1
	MatrixXd LK = (((KX.array().square() + KY.array().square() + 1.0).array()).pow(s)).inverse();

	MatrixXd HK;
	MatrixXd muk;

	return 0;
}