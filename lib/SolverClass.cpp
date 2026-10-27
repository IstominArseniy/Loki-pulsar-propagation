#include <boost/numeric/odeint.hpp>

#include "SolverClass.h"
#include "functions.h"
#include "integrator.h"
#include "constants.h"

#include <fstream>

using Eigen::Vector3d;
using namespace boost::numeric::odeint;

FixedHeightSolver::FixedHeightSolver(double phi_deg, int mode, ObservedRadioPulsar &PSR, FixedHeightModel &model)
{
  model_ = model;
  PSR_ = PSR;
  phi_ = phi_deg * constants::PI / 180;
  mode_ = mode;
  l_ref = 15; // distance, after which the propagation of the O mode can be considered as rectilinear
  r_em = find_emission_point();
  k_em = beta(r_em, vMoment(0));
  r_start, k_start; 
  r_start = translate_rk(r_em, k_em, vMoment(0)).first;
  k_start = translate_rk(r_em, k_em, vMoment(0)).second;
}

// ----------------------- search for emission point ----------------------------------
Vector3d FixedHeightSolver::beta(Vector3d r, Vector3d m){
  Vector3d beta_vec;
  Vector3d Omega_corss_R = PSR_.Rs / constants::c * PSR_.Omega_vec.cross(r);
  Vector3d b_normalized = model_.Bfield(r, m).normalized();
  double kappa = -Omega_corss_R.dot(b_normalized) + std::sqrt(std::pow(Omega_corss_R.dot(b_normalized), 2) - Omega_corss_R.dot(Omega_corss_R) + 1.0 - 1.0/model_.gamma0/model_.gamma0);
  return kappa*b_normalized + Omega_corss_R;
}

Vector3d FixedHeightSolver::find_emission_point()
{ 
  double theta_em = PSR_.chi;
  double phi_em = phi_;
  Vector3d vM = vMoment(0);
  // Vector3d rand_vector = (Vector3d::Random().cross(PSR_.observer_vec)).normalized(); // regularisation of theta_kb
  Vector3d target_vector = PSR_.observer_vec;
  auto func1 = [&](double theta, double phi){
    Vector3d vPoint;
    vPoint << model_.Rem * std::sin(theta) * std::cos(phi), model_.Rem * std::sin(theta) * std::sin(phi), model_.Rem * std::cos(theta);
    Vector3d k = translate_rk(vPoint, beta(vPoint, vM), vM).second;
    return k.cross(target_vector)(0);
  };
  auto func2 = [&](double theta, double phi){
    Vector3d vPoint;
    vPoint << model_.Rem * std::sin(theta) * std::cos(phi), model_.Rem * std::sin(theta) * std::sin(phi), model_.Rem * std::cos(theta);
    Vector3d k = translate_rk(vPoint, beta(vPoint, vM), vM).second;
    return k.cross(target_vector)(1);

  };
  for(int i = 0; i < 15; i ++) {
      double f1x = DX(func1, theta_em, phi_em);
      double f2x = DX(func2, theta_em, phi_em);
      double f1y = DY(func1, theta_em, phi_em);
      double f2y = DY(func2, theta_em, phi_em);
      double f1 = func1(theta_em, phi_em);
      double f2 = func2(theta_em, phi_em);
      double dX = (f1y * f2 - f1 * f2y) / (f1x * f2y - f1y * f2x);
      double dY = (f1x * f2 - f1 * f2x) / (f1y * f2x - f2y * f1x);
      theta_em += dX;
      phi_em += dY;
  }
  Vector3d em_point = {model_.Rem * std::sin(theta_em) * std::cos(phi_em), 
    model_.Rem * std::sin(theta_em) * std::sin(phi_em), 
    model_.Rem * std::cos(theta_em)};
  return em_point;
}
//----------------------------------------------------------------------------------------


// ----------------------- main KO integration -------------------------------------------
std::vector<double> FixedHeightSolver::solve_KO_equations(std::vector<double> theta_initial, double l1, double l2, std::string log_path)
{
  std::vector<double> dep_vars = theta_initial;
  double eps_abs = 1e-6, eps_rel = 1e-6, h_init = 1.0e-8;
  auto addaptive_stepper = make_controlled(eps_abs, eps_rel, runge_kutta_dopri5 < std::vector<double> >()); // make stepper for ODE integration
  auto ode_range = make_adaptive_time_range(addaptive_stepper, [this](const std::vector<double>& f, std::vector<double> &dydx, double l){this->RHS_for_boost(f, dydx, l);}, dep_vars, l1, l2, h_init); // make range for integration
  auto it=ode_range.first;
  if (log_path != ""){ // with logs
    ofstream log_stream(log_path + "/log_" + std::to_string(std::round(100 * phi_*180/constants::PI)/100));
    // log_stream << find_initial_point(true) << std::endl;
    while(it != ode_range.second){ // integration
      if(stop_condition(it->second)){ // TODO: stop integration if Udr > c(?!)
        break;
      }
      // log_stream << it->second <<  ", " << it->first[0] << ", " << it->first[1] << ", " << PSR_.Rs * PSR_.omega_obs / (2.0 * constants::c) * Lambda(it->second) << ", " << BetaBmod (it->second) << ", " <<BetaB (it->second) << ", " << delta (it->second) << ", " << vUdr(it->second)(0) << ", " << vUdr(it->second)(1) << std::endl;
      log_stream << it->second <<  ", " << it->first[0] << ", " << it->first[1] << ", " << BetaBmod (it->second) << ", " << -BetaBmod_derivative(it->second) * constants::c / PSR_.Rs / PSR_.omega_obs / Lambda(it->second) - 1/2/Q(it->second) << std::endl;

      it++; // make integration step
    }
    log_stream.close();

  }
  else{ //no logs
    while(it != ode_range.second){ // integration 
      if(stop_condition(it->second)){ // stop integration if Udr > c(?!)
        break;
      }
      it++; // make integration step
    }
  }
 return dep_vars;
}

bool FixedHeightSolver::stop_condition(double l){
  return false;
}

void FixedHeightSolver::RHS_for_boost(const std::vector<double>& f, std::vector<double> &dydx, double l) {
	double coeff = PSR_.Rs * PSR_.omega_obs / (2.0 * constants::c);
	double LL = Lambda (l);
	double QQ = Q (l);
	double BB = BetaBmod (l);
	dydx[0] = coeff * (-LL / QQ - LL * std::cos(2 * f[0] - 2 * BB) * std::sinh(2 * f[1]));
	dydx[1] = coeff * LL * std::sin(2 * f[0] - 2 * BB) * std::cosh(2 * f[1]);
}

double FixedHeightSolver::find_initial_point(bool use_binary_search) {
  double freq0 = 0.1;
  if(use_binary_search){ // binary search
    double n_iter=0;

    if(std::abs(Lambda_derivative(0) / std::pow(Lambda(0), 2) * 2 * constants::c / PSR_.Rs / PSR_.omega_obs) > freq0)
      return 0; 
    
    double l_left = 0, l_right = PSR_.RLC / 10, l_cur;   
    l_cur = (l_left + l_right) / 2;

    while(std::abs(std::abs(Lambda_derivative(l_cur) / std::pow(Lambda(l_cur), 2) * 2 * constants::c  / PSR_.Rs/ PSR_.omega_obs)  - freq0) > 0.01 && n_iter < 30){
      l_cur = (l_left + l_right) / 2;
      n_iter++;
      if(std::abs(Lambda_derivative(l_cur) / std::pow(Lambda(l_cur), 2) * 2 * constants::c / PSR_.Rs / PSR_.omega_obs) > freq0){
        l_right = l_cur;
      }
      else{
        l_left = l_cur; 
      }
    }
    return l_cur;
  }
  else{ // linear search
    double cr_l = 0, step = 10;
    while(std::abs(Lambda_derivative(cr_l + step) / std::pow(Lambda(cr_l+step), 2) * 2 * constants::c / PSR_.Rs / PSR_.omega_obs) < freq0){
      cr_l += step;
    }
    return cr_l;
  }
}

std::vector<double> FixedHeightSolver::find_approximate_KO_solution(double l)
{
  std::vector<double> theta_final(2);
  double Lambda_integral;
  double starting_point=0;
  Lambda_integral = integrate([this](double l){return this->Lambda(l);}, starting_point, l);
  double coefficient = constants::c / PSR_.Rs / PSR_.omega_obs;
  if (mode_ == 0){
    theta_final[0] = constants::PI / 2 + BetaBmod(l);
    theta_final[1] = (-coefficient / Lambda(l) * BetaBmod_derivative(l) - 1 / 2 / Q(l));
  }
  else{
    theta_final[0] = BetaBmod(l);
    theta_final[1] = -(-coefficient / Lambda(l) * BetaBmod_derivative(l) - 1 / 2 / Q(l));
  }
  return theta_final;
}
//-------------------------------------------------------------------------------------

// ---------------------------------- main functions for KO equations "on ray"--------------------------
Vector3d FixedHeightSolver::vR (double l) {
  /*
  Propagation radius vector
  Here strightforward ray propagation is implemented (but more complex cases of refraction can be considere here too) 
  */
  return r_start + l * k_start;
} 
Vector3d FixedHeightSolver::vK(double l){ 
  return k_start;
}
Vector3d FixedHeightSolver::vMoment (double l) {
  Vector3d mvec;
  mvec(0) = std::sin(PSR_.chi) * std::cos(phi_ + l / PSR_.RLC);
  mvec(1) = std::sin(PSR_.chi) * std::sin(phi_ + l / PSR_.RLC);
  mvec(2) = std::cos(PSR_.chi);
  return mvec;
}
double FixedHeightSolver::delta (double l) {
  Vector3d vr = vR(l);
  Vector3d vm = vMoment(l);
  Vector3d vk = vK(l);
  double vx = vUdr(vr, vm) (0);
  double vy = vUdr(vr, vm) (1);
  double sinth = std::sin(theta_kb(vb(vr, vm), vk));
  double costh = std::cos(theta_kb(vb(vr, vm), vk));
  return std::atan2(-costh * vy, sinth-vx);
}
double FixedHeightSolver::BetaB (double l) {
  Vector3d vr = vR(l);
  Vector3d vm = vMoment(l);
  Vector3d XX;
  Vector3d YY;
  XX = (PSR_.Omega_vec - PSR_.observer_vec.dot(PSR_.Omega_vec) * PSR_.observer_vec).normalized();
  YY = PSR_.observer_vec.cross(XX);
  Vector3d PAvec;
  PAvec = PSR_.observer_vec.cross(vb(vr, vm).cross(PSR_.observer_vec));
  double bx = XX.dot(vB(vr, vm));
  double by = YY.dot(vB(vr, vm));
  return std::atan2(by, bx);  // ??? was atan (by/bx)
}
double FixedHeightSolver::BetaBmod (double l) {
  Vector3d XX;
  Vector3d YY;
  XX = (PSR_.Omega_vec - PSR_.observer_vec.dot(PSR_.Omega_vec) * PSR_.observer_vec).normalized();
  YY = PSR_.observer_vec.cross(XX);
  Vector3d PAvec;
  PAvec = PSR_.observer_vec.cross(vb(vR(l), vMoment(l)).cross(PSR_.observer_vec - beta(vR(l), vMoment(l))));
  double PAx = XX.dot(PAvec);
  double PAy = YY.dot(PAvec);
  return std::atan2(PAy, PAx); 
}
double FixedHeightSolver::Q (double l) {
  Vector3d vr = vR(l);
  Vector3d vm = vMoment(l);
  Vector3d vk = vK(l);
  double vx = vUdr(vr, vm)(0);
  double vy = vUdr(vr, vm)(1);
  double vz = vUdr(vr, vm)(2);
  double costh = std::cos(theta_kb(vb(vr, vm), vk));
  double sinth = std::sin(theta_kb(vb(vr, vm), vk));
  return  model_.lambda * omegaB(vr, vm) * PSR_.omega_obs * (std::pow(sinth - vx, 2) + std::pow(vy * costh, 2)) * model_.Q_type_avrg(std::pow(gammaU(vr, vm) * omegaW(vr, vk, vm) / omegaB(vr, vm), 2))
   / (2 * std::pow(omegaW(vr, vk, vm), 2) * (costh * (1 - vx * vx - vy * vy) - vz * (1.0 - sinth * vx))); // ONLY FOR ZERO TEMPERATURE !!! 
}
double FixedHeightSolver::Lambda (double l) {
  Vector3d vr = vR(l);
  Vector3d vm = vMoment(l);
  Vector3d vk = vK(l);
  double vx = vUdr(vr, vm)(0);
  double vy = vUdr(vr, vm)(1);
  double sinth = std::sin(theta_kb(vb(vr, vm), vk));
  double costh = std::cos(theta_kb(vb(vr, vm), vk));
  double avrg = model_.Lambda_type_avrg(std::pow(gammaU(vr, vm) * omegaW(vr, vk, vm) / omegaB(vr, vm), 2));
  return sgn(avrg) * (-1.0 / 2.0) * std::pow(omegaP(vr, vm) * gammaU(vr, vm) / omegaW(vr, vk, vm), 2) * avrg *
   (std::pow(sinth - vx, 2) + std::pow(vy * costh, 2));
}
double FixedHeightSolver::BetaBmod_derivative(double l)
{
  double dl = 1;
  double sinBeta1 = std::sin(BetaBmod(l));
  double sinBeta2 = std::sin(BetaBmod(l+dl));
  double cosBeta1 = std::cos(BetaBmod(l));
  double cosBeta2 = std::cos(BetaBmod(l+dl));
  double sin_der = (sinBeta2 - sinBeta1) / dl;
  double cos_der = (cosBeta2 - cosBeta1) / dl;
  return std::sqrt(sin_der*sin_der + cos_der*cos_der);
  // return (BetaBmod(l + dl) - BetaBmod(l)) / dl;
}
double FixedHeightSolver::Lambda_derivative(double l)
{
  double dl = 1;
  return (Lambda(l + dl) - Lambda(l)) / dl;
}

//--------------------------------------------------------------------------------------

//------------------------------ refraction --------------------------------------------
std::pair<Vector3d, Vector3d> FixedHeightSolver::translate_rk(Vector3d r0, Vector3d k0, Vector3d m){
  // using 1D approximation
  // double theta_perp0 = ANGLE(m, k0);
  // double r_perp0 = (r0 - r0.dot(m)*m).norm();
  // std::vector<double> f_vec = {theta_perp0, r_perp0};
  // //----------------------integration---------
  // double eps_abs = 1e-6, eps_rel = 1e-6, h_init = 1.0e-8;
  // auto addaptive_stepper = make_controlled(eps_abs, eps_rel, runge_kutta_dopri5 < std::vector<double> >()); // make stepper for ODE integration
  // auto ode_range = make_adaptive_time_range(addaptive_stepper, [this](const vector<double>& f, vector<double> dfdl, double l){this->refraction_RHS_for_boost(f, dfdl, l);}, f_vec, 0, l_ref, h_init); // make range for integration
  // auto it=ode_range.first;
  // while(it != ode_range.second){
  //   it++;
  // }
  // //-------------------------------------------
  // double theta_perp, r_perp;
  // theta_perp = f_vec[0];
  // r_perp = f_vec[1];
  // double delta_theta_perp = theta_perp - theta_perp0;
  // // calculate new r and k
  // Vector3d e = m.cross(k0).normalized();
  // Vector3d r = r_perp/r_perp0 * r0; // TODO: Check if it is a good approximation. Seems not so good
  // Vector3d k = k0 * std::cos(delta_theta_perp) + e.cross(k) * std::sin(delta_theta_perp);
  // // return result
  std::pair<Vector3d, Vector3d> rk_pair = {r0, k0};
  return rk_pair;
}

void FixedHeightSolver::refraction_RHS_for_boost(const std::vector<double>& f, std::vector<double>& dfdl, double l){
  dfdl[0]=1;
  dfdl[1]=1;
}
// double FixedHeightSolver::theta_star(double x, double l_from_star){
//   return 0;
// }

//-----------------------------------------------------------------------------------------

//-------------------------- parameter functions -----------------------------------------
double FixedHeightSolver::psi_m (Vector3d vR, Vector3d vM) {
  return ANGLE(vR, vM);
}
Vector3d FixedHeightSolver::vB (Vector3d vR, Vector3d vM) {
  return model_.Bfield(vR, vM);
}
Vector3d FixedHeightSolver::vb (Vector3d vR, Vector3d vM) {
  return vB(vR, vM).normalized();
}
double FixedHeightSolver::theta_kb (Vector3d vB, Vector3d vK) {
  return ANGLE(vB, vK);
}
Vector3d FixedHeightSolver::vBetaR (Vector3d vR) {
  return PSR_.Rs / constants::c * PSR_.Omega_vec.cross(vR); 
}
Vector3d FixedHeightSolver::vUdr (Vector3d vR, Vector3d vM) {  
  Vector3d vn;
  Vector3d vm;
  vn = (PSR_.observer_vec - PSR_.observer_vec.dot(vb(vR, vM)) * vb(vR, vM)).normalized();
  vm = (vb(vR, vM).cross(vn)).normalized();
  Vector3d temp;
  temp(0) = vBetaR(vR).dot(vn);
  temp(1) = vBetaR(vR).dot(vm);
  if (temp(0)*temp(0) + temp(1) * temp(1) >= 1.0) {
    double coef = 1/std::sqrt(temp(0)*temp(0) + temp(1) * temp(1))*(1e-3);
    temp(0) *= coef;
    temp(1) *= coef;
    temp(2) = 0;
    // throw_error("ERROR: vUdr > 1.");
    std::cout << "Udr > 1" << std::endl;
    return temp;
  }
  temp(2) = std::sqrt(1 - std::pow(temp(0), 2) - std::pow(temp(1), 2));
  return temp;
}
double FixedHeightSolver::gammaU (Vector3d vR, Vector3d vM) {
  double vx = vUdr(vR, vM)(0);
  double vy = vUdr(vR, vM)(1);
  return std::pow(1 - vx * vx - vy * vy, -0.5);
}
double FixedHeightSolver::x_pc(Vector3d vR, Vector3d vM){
  return std::abs(std::sin(psi_m(vR, vM))) * std::sqrt(PSR_.RLC / vR.norm());
}
double FixedHeightSolver::phi_pc(Vector3d vR, Vector3d vM){
  Vector3d m_perp; // Vector, perpendicular to the magnetic axis and e_phi basis vector
  m_perp(0) = -vM(2) *  vM(0) / std::sqrt(std::pow(vM(0), 2) + std::pow(vM(1), 2));
  m_perp(1) = -vM(2) *  vM[1] / std::sqrt(pow(vM(0), 2) + std::pow(vM(1), 2));
  m_perp(2) = std::sqrt(std::pow(vM[0], 2) + std::pow(vM[1], 2));
  Vector3d v_perp = vR - vR.dot(vM) * vM; // projection of vR, perpendicular to the magntic axis
  if((m_perp.cross(v_perp)).norm() >= 0)
    return constants::PI / 2 + ANGLE(v_perp, m_perp); // REDO ANGLE (?!)
  else
    return constants::PI / 2 - ANGLE(v_perp, m_perp);
}
double FixedHeightSolver::gFunc (Vector3d vR, Vector3d vM) {
    return model_.density_profile(x_pc(vR, vM), phi_pc(vR, vM)); // TODO: make 1D/2D differentiation (performace issue)!!!
}

double FixedHeightSolver::Ne(Vector3d vR, Vector3d vM) {
  double nGJ = PSR_.Omega_vec.dot(vB(vR, vM)) * PSR_.B12 * 1e12 / std::pow(vR.norm(), 3) / 
  (2 * constants::PI * constants::c * constants::e);
  return model_.lambda * gFunc(vR, vM) * nGJ;
}

double FixedHeightSolver::omegaB (Vector3d vR, Vector3d vM) {
  return -constants::e * vB(vR, vM).norm() * (PSR_.B12*1e12 / std::pow(vR.norm(), 3)) / (constants::me * constants::c);
}
double FixedHeightSolver::omegaP (Vector3d vR, Vector3d vM) {
  return std::sqrt(4 * constants::PI * constants::e * constants::e * std::abs(Ne(vR, vM)) / constants::me);
}
double FixedHeightSolver::omegaW (Vector3d vR, Vector3d vK, Vector3d vM) {
  double vx = vUdr(vR, vM)(0);
  double vz = vUdr(vR, vM)(2);
  double sinth = std::sin(theta_kb(vB(vR, vM), vK));
  double costh = std::cos(theta_kb(vB(vR, vM), vK));
  return PSR_.omega_obs * (1 - sinth * vx - costh * vz);
}

//-----------------------------------------------------------------------------------------


double FixedHeightSolver::find_intensity()
{
    return gFunc(r_em, vMoment(0)); // TODO: redo; incorrect in the case of refraction
}




// -----------------------Cyclotron Absorption------------------------

double FixedHeightSolver::dtau (double l) {
  double coef = constants::PI * PSR_.Rs / (constants::c * PSR_.omega_obs);
  return coef * std::pow(omegaP(vR(l), vMoment(l)), 2) * model_.fDist(std::abs(omegaB(vR(l), vMoment(l))) / (omegaW(vR(l), vK(l), vMoment(l)) * gammaU(vR(l), vMoment(l))));
}

double FixedHeightSolver::get_tau(){
  return integrate([this](double l){return this->dtau(l);}, 1, PSR_.RLC);
}

// --------------------------- logging ---------------------------------
void FixedHeightSolver::write_params_on_ray(std::string log_path)
{
    if (log_path == "")
      return;
    std::ofstream param_stream(log_path + "/prop_params_" + std::to_string(std::round(100 * phi_*180/constants::PI)/100) + ".csv");
    double l = 0;
    double dl = 2;
    param_stream << "l," << "Lambda," << "delta," << "beta," << "vUdr_x," << "vUdr_y," << "theta_kb," << "x_pc," << "Q," << "beta_eff," << "R/Lambda/2," 
    << "beta_der/Lambda," << "geom_condition," << "p2" << std::endl; 
    while (l <= 1.5 * PSR_.RLC){
        param_stream << l << ",";
        param_stream << Lambda(l) << ",";
        param_stream << delta(l) << ","; 
        param_stream << BetaB(l) << ",";
        param_stream << vUdr(vR(l), vMoment(l))(0) << ",";
        param_stream << vUdr(vR(l), vMoment(l))(1) << ",";
        param_stream << theta_kb(vb(vR(l), vMoment(l)), vK(l)) << ",";
        param_stream << x_pc(vR(l), vMoment(l)) << ",";
        param_stream << Q(l) << ",";
        param_stream << BetaBmod(l) << ", " ;
        param_stream << -BetaBmod_derivative(l) * constants::c / PSR_.Rs / PSR_.omega_obs / Lambda(l) - 1/2/Q(l)<<",";
        param_stream << -BetaBmod_derivative(l) * constants::c / PSR_.Rs / PSR_.omega_obs / Lambda(l)<<",";
        param_stream << std::abs(Lambda_derivative(l) / std::pow(Lambda(l), 2) * 2 * constants::c / PSR_.Rs / PSR_.omega_obs)<<",";
        param_stream << std::pow((PSR_.observer_vec.cross(vb(vR(l), vMoment(l)).cross(PSR_.observer_vec - beta(vR(l), vMoment(l))))).norm(), 2);
        param_stream << std::endl;
        l += dl;
    }
}
// ------------------------------------------------------------------------