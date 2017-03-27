#include <smc_update.h>
#include <math.h>
#include <iostream>

SMCUpdateReturn *smc_update(Pose &pose, OPT &opt, DomainBounds &domain_bounds, ERG &erg, MatrixXd &Ck, double time, double dt) {
	
	double const1 = M_PI / opt.L(0);
	double const2 = M_PI / opt.L(1);

	VectorXd xrel_vec(opt.nagents);
	VectorXd yrel_vec(opt.nagents);
	// Calculating the fourier coefficients of time average statistics distribution
	for (int iagent = 0; iagent < opt.nagents; iagent++) {
		double xrel = pose.x(iagent) - domain_bounds.xmin;
		xrel_vec(iagent) = xrel;
		double yrel = pose.y(iagent) - domain_bounds.ymin;
		yrel_vec(iagent) = yrel;

		MatrixXd temp1 = erg.KX.array().cos().matrix()*const1 * xrel;
		MatrixXd temp2 = erg.KY.array().cos().matrix()*const2 * yrel;
		Ck = Ck + (temp1.cwiseProduct(temp2) * dt).cwiseQuotient(erg.HK.transpose());
	}

	MatrixXd const3 = Ck - opt.nagents*time*erg.muk;
	MatrixXd const4 = -(erg.KX) * M_PI / opt.L(0);
	MatrixXd const5 = -(erg.KY) * M_PI / opt.L(1);
	MatrixXd const6 = erg.KX * M_PI / opt.L(0);
	MatrixXd const7 = erg.KY * M_PI / opt.L(1);

	SMCUpdateReturn *result = (SMCUpdateReturn*)calloc(1, sizeof(SMCUpdateReturn));
	result->x = VectorXd::Zero(opt.nagents);
	result->y = VectorXd::Zero(opt.nagents);
	result->theta = VectorXd::Zero(opt.nagents);
	result->Ck = Ck;

	for (int iagent = 0; iagent < opt.nagents; iagent++) {
		MatrixXd temp1 = const4.cwiseProduct((const6 * xrel_vec(iagent)).array().sin().matrix());
		MatrixXd temp2 = (const7 * yrel_vec(iagent)).array().cos().matrix();
		double Bjx = ((erg.LK).cwiseQuotient((erg.HK).transpose()).cwiseProduct(const3).cwiseProduct(temp1).cwiseProduct(temp2)).sum();

		MatrixXd temp3 = const5.cwiseProduct((const6 * xrel_vec(iagent)).array().cos().matrix());
		MatrixXd temp4 = (const7 * yrel_vec(iagent)).array().sin().matrix();
		double Bjy = ((erg.LK).cwiseQuotient((erg.HK).transpose()).cwiseProduct(const3).cwiseProduct(temp3).cwiseProduct(temp4)).sum();

		double gamma_v = Bjx*cos(pose.theta(iagent)) + Bjy*sin(pose.theta(iagent));
		double gamma_w = -Bjx*sin(pose.theta(iagent)) + Bjy*cos(pose.theta(iagent));

		// Updating agent location based on SMC feedback control law
		double v;
		double w;
		if (gamma_v >= 0)
			v = opt.vlb(iagent);
		else
			v = opt.vub(iagent);

		if (gamma_w >= 0)
			w = opt.wlb(iagent);
		else 
			w = opt.wub(iagent);


		// Velocity motion model 
		if (abs(w) < pow(1, -10)) {
			result->x(iagent) = pose.x(iagent) + v*dt*cos(pose.theta(iagent));
			result->y(iagent) = pose.y(iagent) + v*dt*sin(pose.theta(iagent));
		}
		else {
			result->x(iagent) = pose.x(iagent) + v/w*(sin(pose.theta(iagent) + w*dt) - sin(pose.theta(iagent)));
			result->y(iagent) = pose.y(iagent) + v/w*(cos(pose.theta(iagent)) - cos(pose.theta(iagent)+w*dt));
		}
			result->theta(iagent) = pose.theta(iagent) + w*dt;
	}
	return result;
}