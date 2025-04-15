#include <math.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "read_write.h"
#include "constants.h"
#include "functions.h"
#include "process_functions.h"
#include "initialize.h"
#include "integrator.h" 
using namespace std;
using Eigen::Vector3d;

// most of functions should be local

double sgn (double value) {
  if (value >= 0.0) {
    return 1.0;
  } else {
    return -1.0;
  }
}


Vector3d vMoment (double R) {
  /*
  Magnetic momentum unit vector
  */
  Vector3d mvec;
  mvec(0) = sin(Globals::alpha) * cos(Globals::PHI0 + R / Globals::RLC);
  mvec(1) = sin(Globals::alpha) * sin(Globals::PHI0 + R / Globals::RLC);
  mvec(2) = cos(Globals::alpha);
  return mvec;
} 

Vector3d vR (double R) {
  /*
  Propagation radius vector
  */
  Vector3d n0(3); // unit vector along the ray
  n0(0) = sin(Globals::theta_em) * cos(Globals::phi_em);
  n0(1) = sin(Globals::theta_em) * sin(Globals::phi_em);
  n0(2) = cos(Globals::theta_em);
  return Globals::R_em * n0 + R * Globals::o;
} 

double psi_m (double R) {
  /*
  Angle between magnetic momentum and point on the ray
  */
  return ANGLE(vR(R), vMoment(R));
}


/// @param vR - radius vector to the point
/// @param m - magnetic moment vector
/// @return Magnetic field vector in arbitrary point in magnetosphere 
Vector3d Bfield(Vector3d vR, Vector3d m){
  double Rdist = vR.norm();
  Vector3d n = vR.normalized();
  Vector3d Bdipole = 3 * m.dot(n) * n - m; // Dipole component
  Vector3d Bwind; //Wind component
  Bwind(0) = Rdist/Globals::RLC * Globals::fr * n(0) - std::pow(Rdist/Globals::RLC, 2)* 
  Globals::fphi * Globals::fr * (-n(1));
  Bwind(1) = Rdist/Globals::RLC * Globals::fr * n(1) - std::pow(Rdist/Globals::RLC, 2)* 
  Globals::fphi * Globals::fr * n(0);
  Bwind(2) = Rdist/Globals::RLC * Globals::fr * n(2);
  return Bdipole + Bwind;
}


/// @brief Magneic field on the ray
/// @param R 
/// @return B field vector
Vector3d vB (double R) {
  return Bfield(vR(R), vMoment(R));
}

Vector3d vb (double R) {
  return vB(R).normalized();
}

/// @param R 
/// @return angle between magnetic field line and direction to the observer (in radians)
double theta_kb (double R) {
  return ANGLE(vB(R), Globals::o);
}

Vector3d vBetaR (double R) {
  return constants::R_star / constants::c * Globals::vOmega.cross(vR(R));
}

/// @param R 
/// @return E cross B drift speed vector, normalized to the speed of light
Vector3d vUdr (double R) {  
  Vector3d vn;
  Vector3d vm;
  vn = (Globals::o - Globals::o.dot(vb(R)) * vb(R)).normalized();
  vm = (vb(R).cross(vn)).normalized();

  Vector3d temp;
  temp(0) = vBetaR(R).dot(vn);
  temp(1) = vBetaR(R).dot(vm);
  if (temp(0)*temp(0) + temp(1) * temp(1) >= 1.0) {
    throw_error("ERROR: vUdr > 1.");
  }
  temp(2) = std::sqrt(1 - pow(temp(0), 2) - pow(temp(1), 2));
  return temp;
}

/// @param R 
/// @return delta parameter form Kravtsov-Orlov equations
double delta (double R) {
  double vx = vUdr(R) [0];
  double vy = vUdr(R) [1];
  double sinth = sin(theta_kb(R));
  double costh = cos(theta_kb(R));
  double sign = sgn (- vy * costh / sqrt(pow(sinth - vx, 2) + pow(costh * vy , 2)));
  return sign * acos((sinth - vx) / sqrt(pow(sinth - vx, 2) + pow(costh * vy , 2)));
}


/// @param R 
/// @return angular coordinate of the ray vector in the plane, perpendicular to the field line
double BetaB (double R) {
  Vector3d XX;
  Vector3d YY;
  XX = (Globals::vOmega - Globals::o.dot(Globals::vOmega) * Globals::o).normalized();
  YY = Globals::o.cross(XX);
  double bx = XX.dot(vB(R));
  double by = YY.dot(vB(R));
  return atan(by / bx);
}

/// @param R 
/// @return derivative of the delta parameter along the ray
/// @note distnace variable is normalized to the star radius
double delta_derivative(double R){
  double dl = 1;
  return (delta(R + dl) - delta(R))/dl;
}


/// @param R 
/// @return derivative of the BetaB angle along the ray
/// @note distnace variable is normalized to the star radius
double BetaB_derivative(double R){
  double dl = 1;
  return (BetaB(R + dl) - BetaB(R))/dl;
}

/// @param R 
/// @return distance from the field line footpoint to the polar cap center
/// @note Only for dipolar magnetic field!!!
double r_perp(double R){
  return std::pow(std::sin(psi_m(R)), 2) * Globals::RLC / vR(R).norm();
}

/// @param R 
/// @return polar cap angle cooridnate of the magnetic field line footpoint.
/// Angle is counted from the East-West line on the polar cap
/// @note Applicable only for magnetic fields with zero torsion
double phi_pc(double R){
  Vector3d m_perp; // Vector, perpendicular to the magnetic axis and e_phi basis vector
  m_perp(0) = -vMoment(R)(2) *  vMoment(R)(0) / std::sqrt(std::pow(vMoment(R)(0), 2) + std::pow(vMoment(R)(1), 2));
  m_perp(1) = -vMoment(R)(2) *  vMoment(R)[1] / std::sqrt(pow(vMoment(R)(0), 2) + std::pow(vMoment(R)(1), 2));
  m_perp(2) = std::sqrt(std::pow(vMoment(R)[0], 2) + std::pow(vMoment(R)[1], 2));
  Vector3d v_perp = vR(R) - vR(R).dot(vMoment(R)) * vMoment(R); // projection of vR, perpendicular to the magntic axis
  if((m_perp.cross(v_perp)).norm() >= 0)
    return constants::PI / 2 + ANGLE(v_perp, m_perp); //REDO ANGLE (?!)
  else
    return constants::PI / 2 - ANGLE(v_perp, m_perp);
}

/// @brief Plasma density transvers profile
/// @param R 
/// @return normalized plasma density
double gFunc (double R) {
  if(Globals::density_filename == "default"){
    double f = std::pow(std::sin(psi_m(R)), 2) * Globals::RLC / vR(R).norm();
    double theta = ANGLE(vR(R), Globals::vOmega);
    double dtheta = 5.0 * constants::PI / 180.0;
    double gap = 1.0;
    if (Globals::alpha_deg > 80)
      gap = (1. - exp(-pow(constants::PI * 0.5 - theta, 2) / (2.0 * dtheta * dtheta)));
    return (pow(f, 2.5) * exp(-f * f) / (pow(f, 2.5) + pow(Globals::f0, 2.5))) * gap;
  }
  else{
    return Globals::density_interpolation.get_f(r_perp(R), phi_pc(R));
  }
}


/// @param R 
/// @return Plasma density in physical units (g/cm^3)
double Ne(double R) {
  double nGJ = Globals::vOmega.dot(vB(R)) * Globals::B0 / pow(vR(R).norm(), 3) / 
  (2 * constants::PI * constants::c * constants::e);
  return Globals::lambda * gFunc (R) * nGJ;
}

double RM_density(double R){
    return 2.62e-17 * Ne(R) / Globals::lambda * (Globals::B0 / std::pow(vR(R).norm(), 3)) * std::cos(theta_kb(R));
}

/// @param R 
/// @return Local cyclotron frequency (s^-1)
double omegaB (double R) {
  return -constants::e * vB(R).norm() * (Globals::B0 / std::pow(vR(R).norm(), 3)) / (constants::me * constants::c);
}

double omegaW (double R) {
  double vx = vUdr(R)(0);
  double vz = vUdr(R)(2);
  double sinth = sin(theta_kb(R));
  double costh = cos(theta_kb(R));
  return Globals::omega * (1 - sinth * vx - costh * vz);
}

/// @param R 
/// @return Local plasma frequency (s^-1)
double omegaP (double R) {
  return std::sqrt(4 * constants::PI * constants::e * constants::e * std::fabs(Ne(R)) / constants::me);
}


/// @param R 
/// @return Q parameter from Kravtsov-Orlov equations
double Q (double R) {
  double vx = vUdr(R)(0);
  double vy = vUdr(R)(1);
  double vz = vUdr(R)(2);
  double sinth = std::sin(theta_kb(R));
  double costh = std::cos(theta_kb(R));
  return Globals::lambda * omegaB(R) * Globals::omega * (std::pow(sinth - vx, 2) + std::pow(vy * costh, 2))
   / (2 * std::pow(Globals::gamma0, 3) * std::pow(omegaW(R), 2) * (costh * (1 - vx * vx - vy * vy) - vz * (1.0 - sinth * vx)));
}

/// @brief Particle energy distribution function (both for e- and e+)
/// @param gamma 
/// @return f(gamma) | dN = f(gamma) d gamma
double fDist (double gamma) {
  return ((6.0 * Globals::gamma0) / (pow(2.0, 1.0/6.0) * constants::PI)) * (pow(gamma, 4) / (2.0 * pow(gamma, 6) + pow(Globals::gamma0, 6)));
}


/// @param R 
/// @return drift velocity gamma factor
double gammaU (double R) {
  double vx = vUdr(R)(0);
  double vy = vUdr(R)(1);
  return std::pow(1 - vx * vx - vy * vy, -0.5);
}


/// @brief Auxilary function for Lambda parameter calculation (from Kravtsov-Orlov equations).
/// @param gamma 
/// @param R 
/// @return Integral of certain function
double INTEGRAL (double gamma, double R) {
  double cA = pow(gammaU(R) * omegaW(R) / omegaB(R), 2);
  return -(pow(2, 2.0 / 3.0)*(-2*sqrt(3)*atan((2*pow(2, 1.0 / 3.0)*pow(gamma,2) - pow(Globals::gamma0,2))/
          (sqrt(3)*pow(Globals::gamma0,2))) - 2*log(pow(2, 1.0 / 3.0)*pow(gamma,2) + pow(Globals::gamma0,2)) +
          log(pow(2, 2.0 / 3.0)*pow(gamma,4) - pow(2, 1.0 / 3.0)*pow(gamma,2)*pow(Globals::gamma0,2) +
          pow(Globals::gamma0,4))) - pow(2, 1.0 / 3.0)*pow(Globals::gamma0,2)*cA*
          (2*sqrt(3)*atan((2*pow(2, 1.0 / 3.0)*pow(gamma,2) - pow(Globals::gamma0,2))/(sqrt(3)*pow(Globals::gamma0,2))) -
          2*log(pow(2, 1.0 / 3.0)*pow(gamma,2) + pow(Globals::gamma0,2)) +
          log(pow(2, 2.0 / 3.0)*pow(gamma,4) - pow(2, 1.0 / 3.0)*pow(gamma,2)*pow(Globals::gamma0,2) +
          pow(Globals::gamma0,4))) - 2*pow(Globals::gamma0,4)*pow(cA,2)*
          (log(2*pow(gamma,6) + pow(Globals::gamma0,6)) - 3*log(fabs((1 - gamma*sqrt(cA))*(1 + gamma*sqrt(cA))))))/
          (2.*pow(2, 1.0 / 6.0)*constants::PI*pow(Globals::gamma0,3)*(2 + pow(Globals::gamma0,6)*pow(cA,3)));
}

/// @param R 
/// @return Lambda parameter form Kravtsov-Orlov equations
double Lambda (double R) {
  double vx = vUdr(R)(0);
  double vy = vUdr(R)(1);
  double sinth = std::sin(theta_kb(R));
  double costh = std::cos(theta_kb(R));
  double avrg = INTEGRAL(1000000.0, R) - INTEGRAL(0.0, R);
  return (-1.0 / 2.0) * std::pow(omegaP(R) * gammaU(R) / omegaW(R), 2) * avrg *
   (std::pow(sinth - vx, 2) + std::pow(vy * costh, 2));
}

/// @param R 
/// @return derivative of the Lambda parameter along the ray
/// @note distnace variable is normalized to the star radius
double Lambda_derivative(double R){
    double dl = 1;
    return (Lambda(R + dl) - Lambda(R))/dl;
}

/// @param R 
/// @return differential optical dpeth
/// @note should be integrated along the ray to obtain total optical depth
double dtau (double R) {
  return std::pow(omegaP(R), 2) * fDist (std::fabs(omegaB(R)) / (omegaW(R) * gammaU(R)));
}

/// @brief This function find a distance from emission point where oscillations fade out but p.a. is still strictly
/// following beta + delta. 
/// This point is determined from the condition |Lambda_derivative / Lambda^2 * 2 * c / omega| ~ 1 
/// @param use_binary_search - flag, which determines wither to use binary search or not. Binary search can fail when 
/// Lambda is not monotonous
/// @return point, where integration of Kravtsov-Orlov equations can be started
double find_initial_point(bool use_binary_search) {
  double freq0 = 0.1;
  if(use_binary_search){
    double n_iter=0;
    if(std::fabs(Lambda_derivative(Globals::L_SHIFT) / std::pow(Lambda(Globals::L_SHIFT), 2) * 2 * constants::c / constants::R_star / Globals::omega) > freq0)
      return Globals::L_SHIFT; //shift to avoid zero kB angle
    double R_left = Globals::L_SHIFT, R_right = Globals::RLC / 10, R_cur; 
    R_cur = (R_left + R_right) / 2;
    while(std::fabs(std::fabs(Lambda_derivative(R_cur) / std::pow(Lambda(R_cur), 2) * 2 * constants::c  / constants::R_star/ Globals::omega)  - freq0) > 0.01 && n_iter < 30){
      R_cur = (R_left + R_right) / 2;
      n_iter++;
      if(std::fabs(Lambda_derivative(R_cur) / std::pow(Lambda(R_cur), 2) * 2 * constants::c / constants::R_star / Globals::omega) > freq0){
        R_right = R_cur;
      }
      else{
        R_left = R_cur; 
      }
    }
    return R_cur;
  }

  else{
    double cr_R = Globals::L_SHIFT, step = 10;
    while(std::fabs(Lambda_derivative(cr_R+step) / std::pow(Lambda(cr_R+step), 2) * 2 * constants::c / constants::R_star / Globals::omega) < freq0){
      cr_R += step;
    }
    return cr_R;
  }
}

double approximate_solution_theta0(double R, int mode){
  double integral;
  integral = integrate(Lambda, Globals::L_SHIFT, R);
  double add_term = 0;
  if (gFunc(Globals::L_SHIFT) > 1e-8)
    add_term = - constants::c / constants::R_star / Globals::omega / Lambda(Globals::L_SHIFT) * (BetaB_derivative(Globals::L_SHIFT) + delta_derivative(Globals::L_SHIFT)) * sin(Globals::omega / constants::c * constants::R_star * integral);
  if (mode == 1){
    return BetaB(R) + delta(R) + add_term;
  }
  else
    return constants::PI / 2 + BetaB(R) + delta(R) + add_term;
}

double approximate_solution_theta1(double R, int mode){
  double integral;
  integral = integrate(Lambda, Globals::L_SHIFT, R);
  if (gFunc(Globals::L_SHIFT) < 1e-8)
    return 0;
  if (mode == 1)
    return -constants::c / constants::R_star / Globals::omega / Lambda(R) * (BetaB_derivative(R) + delta_derivative(R)) - 1 / Q(R) + 
  constants::c / constants::R_star / Globals::omega / Lambda(Globals::L_SHIFT) * (BetaB_derivative(Globals::L_SHIFT) + delta_derivative(Globals::L_SHIFT)) * cos(Globals::omega / constants::c * constants::R_star * integral);
  else 
    return constants::c / constants::R_star / Globals::omega / Lambda(R) * (BetaB_derivative(R) + delta_derivative(R)) + 1 / Q(R) - 
  constants::c / constants::R_star / Globals::omega / Lambda(Globals::L_SHIFT) * (BetaB_derivative(Globals::L_SHIFT) + delta_derivative(Globals::L_SHIFT)) * cos(Globals::omega / constants::c * constants::R_star * integral);
}
