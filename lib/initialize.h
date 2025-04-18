#include "interpolator.h"

#pragma INITIALIZE
using Eigen::Vector3d;

namespace Globals {
  extern double theta_em, phi_em,
        B12, B0,
        freqGHz, omega,
        Period, Omega,
        lambda, gamma0, f0,
        mode, fr, fphi,
        alpha_deg, beta_deg, alpha, beta, dzeta,
        PHI0,
        R_em, RLC, RESCAPE, ROMODE, L_SHIFT;
  extern Vector3d vOmega;
  extern Vector3d o;
  extern std::string RUN_ID, input_name, out_path, density_filename;
  extern interpolator2D density_interpolation;
}

void setInitPoints ();

void define_Globals();

void initialize(int argc, char* argv[]);
