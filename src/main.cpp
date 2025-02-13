#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <typeinfo>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>

using namespace std;

#include "../lib/constants.h"
#include "../lib/read_write.h"

#include "../lib/functions.h"
#include "../lib/process_functions.h"
#include "../lib/initialize.h"

#include "../lib/integrator.h"
#include "../lib/RHS.h"

/*
  Default libraries
*/
#include "../lib/diffeqsolver.h"

void displayVector (vector <double> a) {
    cout << endl << a[0] << endl << a[1] << endl << a[2] << endl;
}

int main(int argc, char* argv[]) {
  initialize(argc, argv);

  ofstream output0(Globals::out_path + "/" + Globals::RUN_ID + "_0.dat");
  ofstream output1(Globals::out_path + "/" + Globals::RUN_ID + "_1.dat");
  ofstream output2(Globals::out_path + "/" + Globals::RUN_ID + "_log.dat");
  ofstream output3(Globals::out_path + "/" + Globals::RUN_ID + "_RMs.dat");
  ofstream output4(Globals::out_path + "/" + Globals::RUN_ID + "PAs.dat");



  // SIMULATION STARTS HERE />

  double phi_t_start = read_from_file(Globals::input_name, "phi_start");
  double phi_t_end = read_from_file(Globals::input_name, "phi_end");
  double phi_t_step = read_from_file(Globals::input_name, "phi_step");
  for (double phi_t = phi_t_start; phi_t <= phi_t_end; phi_t += phi_t_step) { // Phase switch
    cout << "PHI: " << phi_t << endl;
    Globals::PHI0 = phi_t * constants::PI / 180.0;
    findInitPoints (Globals::PHI0);
    // cout << r_perp(0) << " " << phi_pc(0) << endl;
    double x1, x2, dep_vars[2];
    // attempt to avoid initial osc. region
    x1 = find_initial_point();
    // x1 = 100;
    cout << x1 << endl;
    // for(int i=800; i<1700; i+=50){
    //   // std::cout<<std::fabs(omegaB(i) / (Globals::gamma0 * gammaU(i) * omegaW(i)))<< " ";
    //   std::cout << omegaW(i) << " ";
    //   //std::cout<<Lambda(i) << " ";
    // }
    // x2 = 1.5 * Globals::RESCAPE;
    x2 = Globals::RLC;
    // Initial values />
    // if (Globals::mode == 0) { // X-mode
    //   dep_vars[0] = BetaB(x1) + delta(x1) + constants::PI / 2.0;
    //   dep_vars[1] = Arcsinh(1.0 / Q(x1)) / 2.0;
    // } else { // O-mode
    //   dep_vars[0] = BetaB(x1) + delta(x1);
    //   dep_vars[1] = Arcsinh(-1.0 / Q(x1)) / 2.0;
    // }
    dep_vars[0] = approximate_solution_theta0(x1, Globals::mode);
    dep_vars[1] = approximate_solution_theta1(x1, Globals::mode);
    // </ Initial values
    double PA = dep_vars[0] * 180 / constants::PI;
    //double RM = integrate(RM_dencity, Globals::R_em, Globals::RLC);
    //std::cout << RM_dencity(Globals::R_em) << " " << RM_dencity(Globals::RLC) << " " << RM << " " << Globals::R_em * 1e6 << std::endl;
    double tau = constants::PI * constants::R_star * integrate(dtau, x1, Globals::RLC) / (constants::c * Globals::omega);
    // testing absorbtion position ---------
    double dx = (Globals::RLC - x1) / 100;
    double tmp_x = x1;

    for(int i=0; i<100; i++){
      output2 << gFunc(tmp_x) << ", ";
      tmp_x += dx;
    }
    output2 << endl;
    output2 << tau << endl;
    // -------------------------------------
    double II0 = gFunc(0);
    double II = II0 * exp (-tau);
    double PA0_rad = dep_vars[0];
    double VV = II * tanh(2.0 * dep_vars[1]);
    output0 << phi_t << " " << II0 << " " << VV << " " << PA << endl;

    int nvar = 2, nok = 0, nbad = 0;
    double deps = 1e-6, h1 = 1.0e-3, hmin = 1.0e-15;
    odeint(dep_vars, nvar, x1, x2, deps, h1, hmin, nok, nbad, RHS);

    VV = II * tanh(2.0 * dep_vars[1]);
    PA = dep_vars[0] * 180.0 / constants::PI;
    // cout << "\tI: " << II << "\n\tV: " << VV << "\n\tPA: " << -PA << endl << endl;
    output1 << phi_t << " " << II << " " << VV << " " << PA << endl ;
    output3 << (dep_vars[0] - PA0_rad) / std::pow((constants::c / Globals::freqGHz / 1e9), 2) << " ";
    output4 << dep_vars[0] << " ";
    //std::cout << (dep_vars[0] - PA0_rad) << " " << dep_vars[0] * 180 / constants::PI << " " <<  PA0_rad * 180 / constants::PI << std::endl;
  }
  output0.close();
  output1.close();
  output3.close();
	return 0;
}
