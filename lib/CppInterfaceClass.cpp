#include <map>
#include <vector>
#include <cmath>
#include <constants.h>
#include <functions.h>
#include <CppModelClass.h>

CppFixedHeightModel::CppFixedHeightModel(std::map<std::string, double> psr_dict, std::map<std::string, double> param_dict){
    std::map<std::string, double> psr_dict_ = psr_dict;
    std::map<std::string, double> param_dict_ = param_dict;
    fr = param_dict['fr'];
    fphi = param_dict['fphi'];
    // resolve parameters 
}

std::map<std::string, double> CppFixedHeightModel::find_ILVPA(double phi, int mode)
{
    calc = FixedHeightCalculator(phi, mode, psr_dict_, param_dict_)
    double l1 = calc.find_initial_point()
    double l2; // DEFINE!!!
    std::vector<double> theta_init = calc.find_approximate_KO_solution(l1);
    std::vector<double> theta_final = calc.solve_KO_equations(theta_init, l1, l2);
    double I = calc.find_intensity();
    double V = I * std::tanh(2.0 * theta_final[1]);
    double PA = theta_final[0] * 180.0 / constants::PI;
    std::map<std::string, double> > result;
    result['I'] = I;
    result['L'] = std::sqrt(I**2 - V**2);
    result['V'] = V;
    result['PA'] = PA;
    return result;
}

// Vector3d CppFixedHeightModel::Bfield(Vector3d vR, Vector3d m){
//   double Rdist = vR.norm();
//   Vector3d n = vR.normalized();
//   Vector3d Bdipole = 3 * m.dot(n) * n - m; // Dipole component
//   Vector3d Bwind; //Wind component
//   Bwind(0) = Rdist/Globals::RLC * Globals::fr * n(0) - std::pow(Rdist/Globals::RLC, 2)* 
//   Globals::fphi * Globals::fr * (-n(1));
//   Bwind(1) = Rdist/Globals::RLC * Globals::fr * n(1) - std::pow(Rdist/Globals::RLC, 2)* 
//   Globals::fphi * Globals::fr * n(0);
//   Bwind(2) = Rdist/Globals::RLC * Globals::fr * n(2);
//   return Bdipole + Bwind;
// }



FixedHeightCalculator::FixedHeightCalculator(double phi, int mode, std::map<std::string, double> psr_dict, std::map<std::string, double> param_dict){
    //----------INIT PSR PARAMS--------------------
    B12 = psr_dict['B12'];
    freqGhz = psr_dict['freqGhz'];
    omega_obs = 2 * constants::PI * freqGhz * 1e9;
    Period = psr_dict['Period'];
    Omega = 2 * constants::PI / Period;
    chi_deg = psr_dict['chi_deg'];
    beta_deg = psr_dict['beta_deg'];
    chi = chi_deg * constants::PI / 180;
    beta = beta_deg * constants::PI / 180;
    dzeta = chi + beta;
    Rs = psr_dict['Rs']
    RLC = (constants::c / Omega) / Rs
    //----------------------------------------------
    //----------INIT MODEL PARAMS-------------------
    fr = param_dict['fr'];
    fphi = param_dict['fphi'];
    Rem = param_dict['Rm'];
    L_SHIFT = param_dict['L_SHIFT'];
    lambda = param_dict['lambda'] // questionable
    gamma0 = param_dict['gamma0'] // questionable
    //----------------------------------------------
    //----------INIT PHASE PARAMS-------------------
    mode_ = mode;
    phi_ = phi;
    observer_vec << std::sin(dzeta), 0, std::cos(dzeta);

    // DEFINE THETA_EM AND PHI_EM !!!
    //----------------------------------------------
}


std::vector<double> FixedHeightCalculator::solve_KO_equations(std::vector<double> theta_initial, double l1, double l2)
{
    return std::vector<double>();
}

std::vector<double> FixedHeightCalculator::find_approximate_KO_solution(double l)
{
    return std::vector<double>();
}

double FixedHeightCalculator::find_initial_point()
{
    return 0.0;
}

std::vector<double> FixedHeightCalculator::find_approximate_KO_solution(std::vector<double> theta_initial)
{
    return std::vector<double>();
}

// "on ray" functions ----------------------------------------------------------
Vector3d FixedHeightCalculator::vMoment (double l) {
  /*
  Magnetic momentum unit vector
  */
  Vector3d mvec;
  mvec(0) = std::sin(chi) * std::cos(phi + l / RLC);
  mvec(1) = std::sin(chi) * std::sin(phi + l / RLC);
  mvec(2) = std::cos(chi);
  return mvec;
}

Vector3d FixedHeightCalculator::vR (double l) {
  /*
  Propagation radius vector
  Here strightforward ray propagation is implemented (but more complex cases of refraction can be considere here too) 
  */
  Vector3d n0(3); // unit vector along the ray
  n0(0) = std::sin(theta_em) * std::cos(phi_em);
  n0(1) = std::sin(theta_em) * std::sin(phi_em);
  n0(2) = std::cos(theta_em);
  return Rem * n0 + l * observer_vec;
} 

double psi_m (double l) {
  /*
  Angle between magnetic momentum and point on the ray
  */
  return ANGLE(vR(l), vMoment(l));
}


Vector3d vB (double l) {
  return Bfield(vR(l), vMoment(l));
}

Vector3d vb (double l) {
  return vB(l).normalized();
}

double theta_kb (double l) {
  return ANGLE(vB(l), observer_vec);
}

Vector3d vBetaR (double l) {
  return Rs / constants::c * vOmega.cross(vR(R)); // TODO vOmega
}

Vector3d vUdr (double l) {  
  Vector3d vn;
  Vector3d vm;
  vn = (observer_vec - observer_vec.dot(vb(l)) * vb(l)).normalized();
  vm = (vb(l).cross(vn)).normalized();

  Vector3d temp;
  temp(0) = vBetaR(R).dot(vn);
  temp(1) = vBetaR(R).dot(vm);
  if (temp(0)*temp(0) + temp(1) * temp(1) >= 1.0) {
    temp(0) = 0;
    temp(1) = 0;
    temp(2) = 0.9;
    return temp;
    // throw_error("ERROR: vUdr > 1.");
  }
  temp(2) = std::sqrt(1 - pow(temp(0), 2) - pow(temp(1), 2));
  return temp;
}





