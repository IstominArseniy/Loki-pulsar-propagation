#include "SolverClass.h"
#include "functions.h"
#include "constants.h"
FixedHeightSolver::FixedHeightSolver(double phi, int mode, ObservedRadioPulsar &PSR, FixedHeightModel &model)
{
    model_ = model
    PSR_ = PSR;
    phi_ = phi;
    mode_ = mode;
}

std::vector<double> FixedHeightSolver::solve_KO_equations(std::vector<double> theta_initial, double l1, double l2)
{
    return std::vector<double>();
}

std::vector<double> FixedHeightSolver::find_approximate_KO_solution(std::vector<double> theta_initial)
{
    return std::vector<double>();
}

std::vector<double> FixedHeightSolver::find_approximate_KO_solution(double l)
{
    return std::vector<double>();
}

double FixedHeightSolver::find_initial_point()
{
    return 0.0;
}

// "on ray" functions ----------------------------------------------------------
Vector3d FixedHeightSolver::vMoment (double l) {
  /*
  Magnetic momentum unit vector
  */
  Vector3d mvec;
  mvec(0) = std::sin(PSR_.chi) * std::cos(phi_ + l / PSR_.RLC);
  mvec(1) = std::sin(PSR_.chi) * std::sin(phi_ + l / PSR_.RLC);
  mvec(2) = std::cos(PSR_.chi);
  return mvec;
}

Vector3d FixedHeightSolver::vR (double l) {
  /*
  Propagation radius vector
  Here strightforward ray propagation is implemented (but more complex cases of refraction can be considere here too) 
  */
  Vector3d n0(3); // unit vector along the ray
  n0(0) = std::sin(theta_em) * std::cos(phi_em); //TODO find theta_em and phi_em
  n0(1) = std::sin(theta_em) * std::sin(phi_em);
  n0(2) = std::cos(theta_em);
  return Rem * n0 + l * PSR_.observer_vec;
} 

double FixedHeightSolver::psi_m (double l) {
  /*
  Angle between magnetic momentum and point on the ray
  */
  return ANGLE(vR(l), vMoment(l));
}

Vector3d FixedHeightSolver::vB (double l) {
  return model_.Bfield(vR(l), vMoment(l));
}

Vector3d FixedHeightSolver::vb (double l) {
  return vB(l).normalized();
}

double FixedHeightSolver::theta_kb (double l) {
  return ANGLE(vB(l), PSR_.observer_vec);
}

Vector3d FixedHeightSolver::vBetaR (double l) {
  return Rs / constants::c * PSR_.Omega_vec.cross(vR(R)); 
}

Vector3d FixedHeightSolver::vUdr (double l) {  
  Vector3d vn;
  Vector3d vm;
  vn = (PSR_.observer_vec - PSR_.observer_vec.dot(vb(l)) * vb(l)).normalized();
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
  temp(2) = std::sqrt(1 - std::pow(temp(0), 2) - std::pow(temp(1), 2));
  return temp;
}

/// @param l 
/// @return drift velocity gamma factor
double gammaU (double l) {
  double vx = vUdr(l)(0);
  double vy = vUdr(l)(1);
  return std::pow(1 - vx * vx - vy * vy, -0.5);
}

/// @param l 
/// @return distance from the field line footpoint to the polar cap center
/// @note Only for dipolar magnetic field!!!
double FixedHeightSolver::r_pc(double l){
  return std::pow(std::sin(psi_m(l)), 2) * PSR_.RLC / vR(l).norm();
}

/// @param l 
/// @return polar cap angle cooridnate of the magnetic field line footpoint.
/// Angle is counted from the East-West line on the polar cap
/// @note Applicable only for magnetic fields with zero torsion
double FixedHeightSolver::phi_pc(double l){
  Vector3d m_perp; // Vector, perpendicular to the magnetic axis and e_phi basis vector
  m_perp(0) = -vMoment(l)(2) *  vMoment(l)(0) / std::sqrt(std::pow(vMoment(l)(0), 2) + std::pow(vMoment(l)(1), 2));
  m_perp(1) = -vMoment(l)(2) *  vMoment(l)[1] / std::sqrt(pow(vMoment(l)(0), 2) + std::pow(vMoment(l)(1), 2));
  m_perp(2) = std::sqrt(std::pow(vMoment(l)[0], 2) + std::pow(vMoment(l)[1], 2));
  Vector3d v_perp = vR(l) - vR(l).dot(vMoment(l)) * vMoment(l); // projection of vR, perpendicular to the magntic axis
  if((m_perp.cross(v_perp)).norm() >= 0)
    return constants::PI / 2 + ANGLE(v_perp, m_perp); //REDO ANGLE (?!)
  else
    return constants::PI / 2 - ANGLE(v_perp, m_perp);
}

/// @brief Plasma density transvers profile
/// @param l 
/// @return normalized plasma density
double FixedHeightSolver::gFunc (double l) {
    return model_.density_proifle(x_pc(l), phi_pc(l)); // TODO make 1D/2D differentiation (performace issue)!!!
}

/// @param l 
/// @return Plasma density in physical units (g/cm^3)
double FixedHeightSolver::Ne(double l) {
  double nGJ = PSR_.Omega_vec.dot(vB(l)) * PSR_.B12 * 1e12 / std::pow(vR(l).norm(), 3) / 
  (2 * constants::PI * constants::c * constants::e);
  return model_.lambda * gFunc(l) * nGJ;
}

/// @param l
/// @return Local cyclotron frequency (s^-1)
double FixedHeightSolver::omegaB (double l) {
  return -constants::e * vB(R).norm() * (PSR_.B12*1e12 / std::pow(vR(l).norm(), 3)) / (constants::me * constants::c);
}

/// @param l
/// @return wave frequency in plasma rest frame
double FixedHeightSolver::omegaW (double l) {
  double vx = vUdr(l)(0);
  double vz = vUdr(l)(2);
  double sinth = std::sin(theta_kb(l));
  double costh = std::cos(theta_kb(l));
  return PSR_.omega * (1 - sinth * vx - costh * vz);
}

/// @param l 
/// @return Local plasma frequency (s^-1)
double omegaP (double l) {
  return std::sqrt(4 * constants::PI * constants::e * constants::e * std::abs(Ne(l)) / constants::me);
}

// --------------------------Main Kravtsov-Orlov equation functions------------------------------------

double FixedHeightSolver::delta (double l) {
  double vx = vUdr(l) (0);
  double vy = vUdr(l) (1);
  double sinth = std::sin(theta_kb(l));
  double costh = std::cos(theta_kb(l));
  double sign = sgn (- vy * costh / std::sqrt(std::pow(sinth - vx, 2) + std::pow(costh * vy , 2)));
  return sign * std::acos((sinth - vx) / std::sqrt(std::pow(sinth - vx, 2) + std::pow(costh * vy , 2)));
}

double FixedHeightSolver::BetaB (double l) {
  Vector3d XX;
  Vector3d YY;
  XX = (PSR_.Omega_vec - PSR_.observer_vec.dot(PSR_.Omega_vec) * PSR_.observer_vec).normalized();
  YY = PSR_.observer_vec.cross(XX);
  double bx = XX.dot(vB(l));
  double by = YY.dot(vB(l));
  return std::atan(by / bx);
}

/// @param l
/// @return Q parameter from Kravtsov-Orlov equations
double FixedHeightSolver::Q (double l) {
  double vx = vUdr(l)(0);
  double vy = vUdr(l)(1);
  double vz = vUdr(l)(2);
  double sinth = std::sin(theta_kb(l));
  double costh = std::cos(theta_kb(l));
  return  model_.lambda * omegaB(l) * PSR_.omega * (std::pow(sinth - vx, 2) + std::pow(vy * costh, 2)) * model_.Q_type_avrg()
   / (2 * std::pow(omegaW(l), 2) * (costh * (1 - vx * vx - vy * vy) - vz * (1.0 - sinth * vx))); // ONLY FOR ZERO TEMPERATURE !!! 
}

/// @param R 
/// @return Lambda parameter form Kravtsov-Orlov equations
double FixedHeightSolver::Lambda (double l) {
  double vx = vUdr(l)(0);
  double vy = vUdr(l)(1);
  double sinth = std::sin(theta_kb(l));
  double costh = std::cos(theta_kb(l));
  double avrg = model_.Lamda_type_avrg(std::pow(gammaU(l) * omegaW(l) / omegaB(l), 2));
  return (-1.0 / 2.0) * std::pow(omegaP(R) * gammaU(R) / omegaW(R), 2) * avrg *
   (std::pow(sinth - vx, 2) + std::pow(vy * costh, 2));
}




