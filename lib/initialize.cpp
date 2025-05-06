#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "read_write.h"
#include "constants.h"
#include "functions.h"
#include "initialize.h"
#include "process_functions.h"
using namespace std;
using Eigen::Vector3d;

// finding initial point of generation for dipole field b0 (along field line)

Vector3d b0 (double th, double ph) {
  Vector3d n0(3);
  n0(0) = sin(th) * cos(ph);
  n0(1) = sin(th) * sin(ph);
  n0(2) = cos(th);
  Vector3d m(3);
  m(0) = sin(Globals::alpha) * cos(Globals::PHI0);
  m(1) = sin(Globals::alpha) * sin(Globals::PHI0);
  m(2) = cos(Globals::alpha);
  return Bfield(n0 * Globals::R_em, m);
}

double func1 (double x, double y) {
  return b0(x, y).cross(Globals::o)(0);
}

double func2 (double x, double y) {
  return b0(x, y).cross(Globals::o)(1);
}

double DX (double (*func)(double, double), double x, double y) {
  double h = 0.00001;
  double fm2 = func(x - 2 * h, y);
  double fp2 = func(x + 2 * h, y);
  double fm1 = func(x - h, y);
  double fp1 = func(x + h, y);
  return (fm2 - 8 * fm1 + 8 * fp1 - fp2) / (12 * h);
}
double DY (double (*func)(double, double), double x, double y) {
  double h = 0.00001;
  double fm2 = func(x, y - 2 * h);
  double fp2 = func(x, y + 2 * h);
  double fm1 = func(x, y - h);
  double fp1 = func(x, y + h);
  return (fm2 - 8 * fm1 + 8 * fp1 - fp2) / (12 * h);
}

/// @brief Function, which finds emission point sperical coordinates and put them into corresponding global variables 
void setInitPoints () {
  double X = Globals::alpha, Y = Globals::PHI0;
  for(int i = 0; i < 50; i ++) {
      double f1x = DX(func1, X, Y);
      double f2x = DX(func2, X, Y);
      double f1y = DY(func1, X, Y);
      double f2y = DY(func2, X, Y);
      double f1 = func1(X, Y);
      double f2 = func2(X, Y);
      double dX = (f1y * f2 - f1 * f2y) / (f1x * f2y - f1y * f2x);
      double dY = (f1x * f2 - f1 * f2x) / (f1y * f2x - f2y * f1x);
      X += dX;
      Y += dY;
  }
  Globals::theta_em = X;
  Globals::phi_em = Y;
  
}

// initialize Globals >
double Globals::theta_em, Globals::phi_em,
      Globals::B12, Globals::B0,
      Globals::freqGHz, Globals::omega,
      Globals::Period, Globals::Omega,
      Globals::lambda, Globals::gamma0, Globals::f0,
      Globals::mode, Globals::fr, Globals::fphi,
      Globals::alpha_deg, Globals::beta_deg, Globals::alpha, Globals::beta, Globals::dzeta,
      Globals::PHI0,
      Globals::R_em, Globals::RLC, Globals::RESCAPE, Globals::ROMODE, Globals::L_SHIFT;
Vector3d Globals::vOmega;
Vector3d Globals::o;
interpolator2D Globals::density_interpolation;

string Globals::RUN_ID, Globals::input_name, Globals::out_path, Globals::density_filename;
// < initialize Globals

void define_Globals() {
  Globals::B12 = read_from_file(Globals::input_name, "B12"); // Surface B-field in 10^12 Gs
  Globals::Period = read_from_file(Globals::input_name, "Period"); // Rotation period in sec
  Globals::freqGHz = read_from_file(Globals::input_name, "freqGHz"); // Radiation frequency in GHz

  Globals::lambda = read_from_file(Globals::input_name, "lambda"); // Plasma multiplicity // normally 10000
  Globals::gamma0 = read_from_file(Globals::input_name, "gamma0"); // Plasma mean gamma factor // normally 50
  Globals::f0 = read_from_file(Globals::input_name, "f0"); // Polar Cap gap width // normally 0.5
  Globals::R_em = read_from_file(Globals::input_name, "R_em"); // Emission radius in star radii

  Globals::mode = read_from_file(Globals::input_name, "mode"); // 1 = O-mode & 0 = X-mode
  Globals::fr = read_from_file(Globals::input_name, "fr"); // Split-monopole parameter
  Globals::fphi = read_from_file(Globals::input_name, "fphi"); // Split-monopole parameter

  Globals::alpha_deg = read_from_file(Globals::input_name, "alpha_deg");
  Globals::beta_deg = read_from_file(Globals::input_name, "beta_deg");
  Globals::density_filename = read_from_file_str(Globals::input_name, "density_filename", "default");
  if(Globals::density_filename != "default"){
    Globals::density_interpolation.init(Globals::density_filename, 0, 1, 0, 2 * constants::PI);
  }
  Globals::alpha = Globals::alpha_deg * constants::PI / 180.0; // Inclination angle in radians
  Globals::beta = Globals::beta_deg * constants::PI / 180.0; // Line of sight angle in radians;
  Globals::dzeta = Globals::alpha - Globals::beta; // Minimum angle between the rotation axis and the line of sight

  Globals::B0 = Globals::B12 * 1.0e12; // Surface magnetic field in Gs
  Globals::Omega = 2.0 * constants::PI / Globals::Period; // Rotation frequency
  Globals::omega = 2.0 * constants::PI * Globals::freqGHz * 1.0e9; // Radiation circular frequency
  Globals::L_SHIFT = 1.0; // Shift to avoid zero k_B angle

  Globals::vOmega << 0, 0, Globals::Omega;
  Globals::o << std::sin(Globals::dzeta), 0, std::cos(Globals::dzeta);

  Globals::RLC = (constants::c / Globals::Omega) / constants::R_star;
  Globals::RESCAPE = 1.0e3 * pow(Globals::lambda / 1.0e4, 1.0/3.0) * pow(Globals::gamma0 / 100.0, -6.0/5.0) * pow(Globals::B0 / 1.0e12, 2.0/5.0) * pow(Globals::freqGHz, -2.0/5.0) * pow(Globals::Period, -1.0/5.0);
  Globals::ROMODE = 1.0e2 * pow(Globals::lambda / 1.0e4, 1.0/3.0) * pow(Globals::gamma0 / 100.0, 1.0/3.0) * pow(Globals::B0 / 1.0e12, 1.0/3.0) * pow(Globals::freqGHz, -2.0/3.0) * pow(Globals::Period, -1.0/3.0);
}

void initialize(int argc, char* argv[]) {
  read_in_out(Globals::input_name, Globals::out_path, argc, argv);

  define_Globals();

  Globals::RUN_ID = read_from_file_str(Globals::input_name, "run_id", "my_run");


  string MODE;
  if (Globals::mode == 0) MODE = "X-mode";
  else MODE = "O-mode";

  // Create output directory if doesn't exist
  struct stat st = {0};
  if (stat(Globals::out_path.c_str(), &st) == -1) {
    mkdir(Globals::out_path.c_str(), 0700);
  }

  ofstream outputData(Globals::out_path + "/" + Globals::RUN_ID + ".dat");
  outputData
      << "alpha = " << Globals::alpha_deg
      << "\nbeta = " << Globals::beta_deg
      << "\n\nPeriod = " << Globals::Period
      << "\nB12 = " << Globals::B12
      << "\nfGHz = " << Globals::freqGHz
      << "\n\nlambda = " << Globals::lambda
      << "\ngamma0 = " << Globals::gamma0
      << "\nf0 = " << Globals::f0
      << "\nR_em = " << Globals::R_em
      << "\n\n" + MODE
      << "\nfr = " << Globals::fr
      << "\nfphi = " << Globals::fphi
      << "\n\n\nR_LC = " << Globals::RLC
      << "\nR_escape = " << Globals::RESCAPE
      << "\nR_A = " << Globals::ROMODE;
  outputData.close();
}
