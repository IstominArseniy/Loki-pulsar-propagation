#include <boost/numeric/odeint.hpp>

#include "SolverClass.h"
#include "functions.h"
#include "integrator.h"
#include "constants.h"

using Eigen::Vector3d;
using namespace boost::numeric::odeint;

FixedHeightSolver::FixedHeightSolver(double phi, int mode, ObservedRadioPulsar &PSR, FixedHeightModel &model)
{
    model_ = model
    PSR_ = PSR;
    phi_ = phi;
    mode_ = mode;
    std::pair<double, double> emission_coords = find_emission_point();
    theta_em_ = emission_coords[0];
    phi_em_ = emission_coords[1];
}

std::vector<double> FixedHeightSolver::solve_KO_equations(std::vector<double> theta_initial, double l1, double l2)
{
  std::vector<double> dep_vars = theta_initial;
  double eps_abs = 1e-6, eps_rel = 1e-6, h_init = 1.0e-3; // accuracy and initial step
  auto addaptive_stepper = make_controlled(eps_abs, eps_rel, runge_kutta_dopri5 < std::vector<double> >()); // make stepper for ODE integration
  auto ode_range = make_adaptive_time_range(addaptive_stepper, RHS_for_boost, dep_vars, l1, l2, h_init); // make range for integration
  auto it=ode_range.first;
  while(it != ode_range.second){ // integration 
    if(stop_condition(it->second)){ // stop integration if Udr > c(?!)
      break;
    }
    //logs - REDO
    // plot << it->second <<  ", " << it->first[0] << ", " << it->first[1] << ", " << constants::R_star * Globals::omega / (2.0 * constants::c) * Lambda(it->second) << ", " << BetaB (it->second) + delta (it->second) << ", " << gFunc(it->second) << std::endl;
    it++; // make integration step
  }
  return dep_vars;
}

void FixedHeightSolver::RHS_for_boost(const std::vector<double>& f, std::vector<double> &dydx, double l) {
	double coeff = PSR_.R_star * model_.omega / (2.0 * constants::c);

	double LL = Lambda (l);
	double QQ = Q (l);
	double BB = BetaB (l);
	double DD = delta (l);

	dydx[0] = coeff * (-LL / QQ - LL * std::cos(2 * f[0] - 2 * BB - 2 * DD) * std::sinh(2 * f[1]));
	dydx[1] = coeff * LL * std::sin(2 * f[0] - 2 * BB - 2 * DD) * std::cosh(2 * f[1]);
}

std::vector<double> FixedHeightSolver::find_approximate_KO_solution(std::vector<double> theta_initial)
{
  // SEEMS TO BE IMPOSSIBLE
  return std::vector<double>();
}

std::vector<double> FixedHeightSolver::find_approximate_KO_solution(double l)
{
  std::vector<double> theta_final(2);
  double Lambda_integral;
  Lambda_integral = integrate(Lambda, model_.L_SHIFT, l);
  double coefficient = constants::c / constants::R_star / Globals::omega
  double delta_theta = - coefficient / Lambda(model_.L_SHIFT) * (get_derivative_along_the_ray(BetaB, model_.L_SHIFT) + get_derivative_along_the_ray(delta, model_.L_SHIFT)) * std::sin(integral / coefficient);
  double V = -coefficient / Lambda(l) * (get_derivative_along_the_ray(BetaB, l) + get_derivative_along_the_ray(delta, l)) - 1 / Q(l) + 
  coefficient / Lambda(model_.L_SHIFT) * (get_derivative_along_the_ray(BetaB, model_.L_SHIFT) + get_derivative_along_the_ray(delta, model_.L_SHIFT)) * std::cos(integral / coefficient);
  if (mode == 0){
    theta_final[0] = constants::PI / 2 + BetaB(l) + delta(l) + delta_theta;
    theta_final[1] = V;
  }
  else{
    theta_final[0] = BetaB(l) + delta(l) + delta_theta;
    theta_final[1] = -V;
  }
  return std::vector<double>();
}


double FixedHeightSolver::find_initial_point(bool use_binary_search) {
  double freq0 = 0.1;
  if(use_binary_search){ // binary search
    double n_iter=0;

    if(std::abs(get_derivative_along_the_ray(Lambda, model_.L_SHIFT) / std::pow(Lambda(model_.L_SHIFT), 2) * 2 * constants::c / PSR_.Rs / PSR_.omega) > freq0)
      return model_.L_SHIFT; //shift to avoid zero kB angle

    double l_left = model_.L_SHIFT, l_right = PSR_.RLC / 10, l_cur;   
    l_cur = (l_left + l_right) / 2;

    while(std::abs(std::abs(get_derivative_along_the_ray(Lambda, l_cur) / std::pow(Lambda(l_cur), 2) * 2 * constants::c  / PSR_.Rs/ PSR_.omega)  - freq0) > 0.01 && n_iter < 30){
      l_cur = (l_left + l_right) / 2;
      n_iter++;
      if(std::abs(get_derivative_along_the_ray(Lambda, l_cur) / std::pow(Lambda(l_cur), 2) * 2 * constants::c / PSR_.Rs / PSR_.omega) > freq0){
        l_right = l_cur;
      }
      else{
        l_left = l_cur; 
      }
    }
    return l_cur;
  }
  else{ // linear search
    double cr_l = model_.L_SHIFT, step = 10;
    while(std::abs(get_derivative_along_the_ray(Lambda, cr_l+step) / std::pow(Lambda(cr_l+step), 2) * 2 * constants::c / PSR_.Rs / PSR_.omega) < freq0){
      cr_l += step;
    }
    return cr_l;
  }
}

double FixedHeightSolver::find_intensity()
{
    return gFunc(l)
}

// ------------------------------------------------------------"on ray" functions ----------------------------------------------------------
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

std::pair<double, double> FixedHeightSolver::find_emission_point()
{
    return std::pair<double, double>();
}

Vector3d FixedHeightSolver::vR (double l) {
  /*
  Propagation radius vector
  Here strightforward ray propagation is implemented (but more complex cases of refraction can be considere here too) 
  */
  Vector3d n0(3); // unit vector along the ray
  n0(0) = std::sin(theta_em_) * std::cos(phi_em_); 
  n0(1) = std::sin(theta_em_) * std::sin(phi_em_);
  n0(2) = std::cos(theta_em_);
  return Rem * n0 + l * PSR_.observer_vec;
} 

double FixedHeightSolver::psi_m (double l) {
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
  return PSR_.Rs / constants::c * PSR_.Omega_vec.cross(vR(R)); 
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

double FixedHeightSolver::gammaU (double l) {
  double vx = vUdr(l)(0);
  double vy = vUdr(l)(1);
  return std::pow(1 - vx * vx - vy * vy, -0.5);
}

double FixedHeightSolver::r_pc(double l){
  return std::pow(std::sin(psi_m(l)), 2) * PSR_.RLC / vR(l).norm();
}

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

double FixedHeightSolver::gFunc (double l) {
    return model_.density_proifle(x_pc(l), phi_pc(l)); // TODO make 1D/2D differentiation (performace issue)!!!
}

double FixedHeightSolver::Ne(double l) {
  double nGJ = PSR_.Omega_vec.dot(vB(l)) * PSR_.B12 * 1e12 / std::pow(vR(l).norm(), 3) / 
  (2 * constants::PI * constants::c * constants::e);
  return model_.lambda * gFunc(l) * nGJ;
}

double FixedHeightSolver::omegaB (double l) {
  return -constants::e * vB(R).norm() * (PSR_.B12*1e12 / std::pow(vR(l).norm(), 3)) / (constants::me * constants::c);
}

double FixedHeightSolver::omegaW (double l) {
  double vx = vUdr(l)(0);
  double vz = vUdr(l)(2);
  double sinth = std::sin(theta_kb(l));
  double costh = std::cos(theta_kb(l));
  return PSR_.omega * (1 - sinth * vx - costh * vz);
}

double FixedHeightSolver::omegaP (double l) {
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

double FixedHeightSolver::Q (double l) {
  double vx = vUdr(l)(0);
  double vy = vUdr(l)(1);
  double vz = vUdr(l)(2);
  double sinth = std::sin(theta_kb(l));
  double costh = std::cos(theta_kb(l));
  return  model_.lambda * omegaB(l) * PSR_.omega * (std::pow(sinth - vx, 2) + std::pow(vy * costh, 2)) * model_.Q_type_avrg()
   / (2 * std::pow(omegaW(l), 2) * (costh * (1 - vx * vx - vy * vy) - vz * (1.0 - sinth * vx))); // ONLY FOR ZERO TEMPERATURE !!! 
}

double FixedHeightSolver::Lambda (double l) {
  double vx = vUdr(l)(0);
  double vy = vUdr(l)(1);
  double sinth = std::sin(theta_kb(l));
  double costh = std::cos(theta_kb(l));
  double avrg = model_.Lamda_type_avrg(std::pow(gammaU(l) * omegaW(l) / omegaB(l), 2));
  return (-1.0 / 2.0) * std::pow(omegaP(R) * gammaU(R) / omegaW(R), 2) * avrg *
   (std::pow(sinth - vx, 2) + std::pow(vy * costh, 2));
}

