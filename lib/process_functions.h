#pragma PROCESS_FUNCTIONS
using namespace std;

vector <double> vB(double R);
vector <double> vBsplit(double R);

double delta (double R);
double BetaB (double R);
double Q (double R);
double Lambda (double R);

double dtau (double R);
double psi_m (double R);
vector <double> vR (double R);
double r_perp(double R);
double phi_pc(double R);
double gFunc (double R);

vector <double> vMoment (double R);
vector <double> vUdr (double R);
double theta_kb (double R);
double find_initial_point();
double dLambda_dl(double R);
double RM_dencity(double R);
double approximate_solution_theta0(double R, int mode);
double approximate_solution_theta1(double R, int mode);