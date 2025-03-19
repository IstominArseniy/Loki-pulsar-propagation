#pragma PROCESS__UNCTIONS
#include <Eigen/Dense>
using Eigen::Vector3d;


Vector3d vB(double R);
Vector3d vBsplit(double R);

double delta (double R);
double BetaB (double R);
double delta_derivative(double R);
double BetaB_derivative(double R);
double Q (double R);
double Lambda (double R);
double omegaB (double R);
double gammaU (double R);
double omegaW (double R);
double dtau (double R);
double psi_m (double R);
double r_perp(double R);
double phi_pc(double R);
double gFunc (double R);

Vector3d vR (double R);
Vector3d vMoment (double R);
Vector3d vUdr (double R);
double theta_kb (double R);
double find_initial_point(bool use_binary_search=true);
double Lambda_derivative(double R);
double RM_dencity(double R);

double approximate_solution_theta0(double R, int mode);
double approximate_solution_theta1(double R, int mode);
