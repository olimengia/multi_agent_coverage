/*
    initialization.cpp:
        Initialzie the enviroment for the agents to perform in
*/

#include <initialization.h>
#include <AgentClass.h>
#include <calculate_ergodicity.h>
#include <generate_utility_map.h>
#include <get_fourier_coeff.h>
#include <smc_update.h>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace Eigen;

// Domain bounds
#define DOMAIN_XMIN 0.0
#define DOMAIN_XMAX 300.0
#define DOMAIN_YMIN 0.0
#define DOMAIN_YMAX 300.0
#define NUM_AGENTS 3
#define NSTEPS 3
#define DT 0.1

int main(void) {

    double Lx = DOMAIN_XMAX - DOMAIN_XMIN;
    double Ly = DOMAIN_YMAX - DOMAIN_YMIN;

    VectorXd vlb_init = VectorXd::Ones(NUM_AGENTS) * 0.1;
    VectorXd vub_init = VectorXd::Ones(NUM_AGENTS) * 5;
    VectorXd wlb_init = VectorXd::Ones(NUM_AGENTS) * -0.5;
    VectorXd wub_init = VectorXd::Ones(NUM_AGENTS) * 0.5;

    Vector3d x;
    x << 100, 200, 50;
    Vector3d y;
    y << 50, 75, 250;
    Vector3d thetas;
    thetas << M_PI, M_PI, M_PI;

    // Agents params: Velocity Bounds, number of agents...

    /********* Initialization ends *************/
    // Initialize the agent locations and limits
    /*Agent agents[NUM_AGENTS];

    for (int i = 0; i < NUM_AGENTS; i++) {
        agents[i].set_x(x(i));
        agents[i].set_y(y(i));
        agents[i].set_theta(thetas(i));
        agents[i].set_v_lower(vlb_init);
        agents[i].set_v_upper(vub_init);
        agents[i].set_w_lower(wlb_init);
        agents[i].set_w_upper(wub_init);
    }*/

    // Initialize ergodicity params
    double s = 1.5;
    int Nkx = 50;
    int Nky = 50;
    MatrixXd KX = (VectorXd::LinSpaced(Nkx, 0, Nkx-1)) * RowVectorXd::Ones(Nky);   // Nkx x 1
    MatrixXd KY = VectorXd::Ones(Nky) * RowVectorXd::LinSpaced(Nky, 0, Nky-1);  // Nky x 1
    MatrixXd LK = (((KX.array().square() + KY.array().square() + 1.0).array()).pow(s)).matrix().inverse();

    /********* Initialization ends *************/
    std::cout << "Initialization complete" << std::endl;

    int add_noise = 0;
    double xdel = 1;    // Resolution of x and y 
    double ydel = 1;

    MatrixXd mu = generate_utility_map(DOMAIN_XMIN, DOMAIN_XMAX, xdel, DOMAIN_YMIN, DOMAIN_YMAX, ydel, add_noise);
    
    std::ofstream outfile;

    // Normalize information distribution
    mu = mu / mu.sum();
    /*outfile.open("mu.txt");
    outfile << "normalized mu is: "<< std::endl;
    outfile << mu << std::endl;
    outfile.close();*/
    MatrixXd X = RowVectorXd::LinSpaced(mu.cols(), DOMAIN_XMIN, DOMAIN_XMAX).replicate(mu.rows(),1);
    MatrixXd Y = VectorXd::LinSpaced(mu.rows(), DOMAIN_YMIN, DOMAIN_YMAX).replicate(1,mu.cols());

    FourierReturn *fourier_coeff = get_fourier_coeff(mu, X, Y, Nkx, Nky, Lx, Ly);

    // mu is used for plotting

    /****** Simulation starts ******/
    std::cout << "Simulation starts" << std::endl;    

    // Initialize Fourier coefficients of coverage distribution 
    MatrixXd Ck = MatrixXd::Zero(Nkx, Nky);
    VectorXd ergodicity_metric = VectorXd::Zero(NSTEPS, 1);

    // Plot positions of agents initially 

    Matrix<double, NSTEPS, NUM_AGENTS> *positions = new Matrix<double, NSTEPS, NUM_AGENTS>[3];

    // Pack L and number of agents into struct to pass into smc_update
    Vector2d L;
    L << Lx, Ly;
    OPT opt;
    opt.L = L;
    opt.nagents = NUM_AGENTS;
    opt.vlb = vlb_init;
    opt.vub = vub_init;
    opt.wlb = wlb_init;
    opt.wub = wub_init;

    // Pack pose and opt into struct to pass into smc_update
    Pose pose;
    pose.x = x;
    pose.y = y;
    pose.theta = thetas;

    // Pack domain bounds into struct to pass into smc_update
    DomainBounds domain_bounds;
    domain_bounds.xmin = DOMAIN_XMIN;
    domain_bounds.xmax = DOMAIN_XMAX;
    domain_bounds.ymin = DOMAIN_YMIN;
    domain_bounds.ymax = DOMAIN_YMAX;

    // Pack ergodicity parameters into struct to pass into smc_update
    ERG erg;
    erg.KX = KX;
    erg.KY = KY;
    erg.LK = LK;
    erg.HK = fourier_coeff->HK;
    erg.muk = fourier_coeff->muk;


    // Hadi said "alway start from 0"
    for (int n = 0; n < NSTEPS; n++) {
        double time = n * DT;
        SMCUpdateReturn *smc_result = smc_update(pose, opt, domain_bounds, erg, Ck, time, DT);
        for (int iagent = 0; iagent < NUM_AGENTS; iagent++) {
            (positions[0])(n, iagent) = (smc_result->x)(iagent);
            (positions[1])(n, iagent) = (smc_result->y)(iagent);
        }
        ergodicity_metric(n) = calculate_ergodicity(Ck, LK, fourier_coeff->muk);
    }
    
    return 0;
}
