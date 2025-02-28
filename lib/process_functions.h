#pragma PROCESS__UNCTIONS
using namespace std;

vector <double> vB(double R);
vector <double> vBsplit(double R);

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
vector <double> vR (double R);
double r_perp(double R);
double phi_pc(double R);
double gFunc (double R);

vector <double> vMoment (double R);
vector <double> vUdr (double R);
double theta_kb (double R);
double find_initial_point(bool use_binary_search=true);
double Lambda_derivative(double R);
double RM_dencity(double R);

double approximate_solution_theta0(double R, int mode);
double approximate_solution_theta1(double R, int mode);
