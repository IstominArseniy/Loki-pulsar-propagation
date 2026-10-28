#include "PSRClass.h"
#include "FixedHeightModelClass.h"
#include <Eigen/Dense>


#pragma once

using Eigen::Vector3d;

class FixedHeightSolver{
    private:
    double phi_;
    double mode_;
    ObservedRadioPulsar PSR_;
    double l_ref;
    Vector3d r_em;
    Vector3d k_em;
    Vector3d r_start;
    Vector3d k_start;
    FixedHeightModel model_;

    public:
    FixedHeightSolver(double phi_deg, int mode, ObservedRadioPulsar &PSR, FixedHeightModel &model);

    /// @brief solve Kravtsov-Orlov equations 
    /// @param theta_inital inital values for comlex polarisation angle (as length 2 vector)
    /// @param l1 - initial point
    /// @param l2 - final point
    /// @param mode emission mode (1 - Omode, 0 - Xmode)
    /// @return theta_final -  final complex polarisation angle (as length 2 vector)
    std::vector<double> solve_KO_equations(std::vector<double> theta_initial, double l1, double l2, std::string log_path="");

    /// @brief RHS for solving Kravtsov-Orlov equations using <boost/numeric/odeint.hpp> integration routines
    /// @param f RHS function
    /// @param dydx derivative by reference
    /// @param l independent variable by reference
    void RHS_for_boost(const std::vector<double>& f, std::vector<double> &dydx, double l);

    void refraction_RHS(const std::vector<double>& f, std::vector<double>& dfdl, double l, Vector3d vM);

    double k_div_n(Vector3d vR, Vector3d vK, Vector3d vM);

    Vector3d dr_derivative(Vector3d vR, Vector3d vK, Vector3d vM);

    Vector3d dk_derivative(Vector3d vR, Vector3d vK, Vector3d vM);

    void refraction1D_RHS(const std::vector<double> &f, std::vector<double> &dfdl, double l);



    /// @param vR point vector
    /// @param vM magnetic moment vector
    /// @brief theta_star^4 parameter for O-mode refraction index
    double theta_star4(Vector3d vR, Vector3d vM);

    /// @param vR point vector
    /// @param vK wave vector
    /// @param vM magnetic moment vector
    /// @brief O mode refrative index
    double nO(Vector3d vR, Vector3d vK, Vector3d vM);

    /// @brief find approximate solution for Kravtsov-Orlov equations in dense plasma regions 
    /// @param l - distance from the star. Theta_inintal is determined by emission mode
    /// @param mode emission mode (1 - Omode, 0 - Xmode)
    /// @return theta_final -  final complex polarisation angle (as length 2 vector)
    std::vector<double> find_approximate_KO_solution(double l);

    /// @brief This function finds a distance from emission point where oscillations fade out but p.a. is still strictly
    /// following beta + delta. 
    /// This point is determined from the condition |Lambda_derivative / Lambda^2 * 2 * c / omega| ~ 1 
    /// @param use_binary_search - flag, which determines wither to use binary search or not. Binary search can fail when 
    /// Lambda is not monotonous
    /// @return point, where integration of Kravtsov-Orlov equations can be started
    double find_initial_point(bool use_binary_search);

    /// @brief find emission intensity
    /// @return emission intensity in arbitrary uits
    double find_intensity();

    /// @return tau - optical thickness (I = I0 * e^(-tau))
    double get_tau();

    void write_params_on_ray(std::string log_path);



    // --------------------------Main Kravtsov-Orlov equation functions------------------------------------

    /// @brief additional phase parameter related to ExB drift 
    /// @param l distance alnog the ray 
    /// @return delta phase 
    double delta (double l);

    /// @brief polar angle of magnetic filed projected to the plane perpendicular to the ray
    /// @param l distnace along the ray
    /// @return BetaB angle in radians
    double BetaB (double l);

    /// @brief equvivalent of BetaB + delta 
    /// @param l distnace along the ray
    /// @return BetaB + delta angle in radians
    double BetaBmod (double l);

    /// @param l distance alnog the ray
    /// @return Q parameter from Kravtsov-Orlov equations
    double Q (double l);

    /// @param l distance along the ray
    /// @return Lambda parameter form Kravtsov-Orlov equations
    double Lambda (double l);
    // -----------------------------------------------------------------"On ray" functions-------------------------------------------------
    private:

    /// @param l distance along the ray
    /// @return magnetic moment unit vector (Egien 3d Vector)
    Vector3d vMoment (double l);
    
    /// @param r - radius vector
    /// @param m - magnetic moment unit vector
    /// @return particles velocity vector
    Vector3d beta(Vector3d r, Vector3d m);

    /// @brief Evolve r and k vectors, solving geometrial optics equations
    /// @param r0 initial wave r
    /// @param k0 initial wave k
    /// @param m magnetic moment unit vector
    /// @return pair of final (r, k) vector
    std::pair<Vector3d, Vector3d>  translate_rk(Vector3d r0, Vector3d k0, Vector3d vM);

    std::pair<Vector3d, Vector3d> translate_rk_1D(Vector3d r0, Vector3d k0, Vector3d m);


    /// @brief finds spherical angular coordinates of the emission point (on the distance Rem) using ...
    /// @return 3d Vector to the emission point
    Vector3d find_emission_point();

    /// @param l distance along the ray
    /// @return vector to the point on the ray as Eigen 3d Vector
    Vector3d vR (double l);

    Vector3d vK(double l);

    /// @param vR point vector
    /// @param vM magnetic moment vector
    /// @return angle between magnetic axis and r vector
    double psi_m(Vector3d vR, Vector3d vM);

    /// @param vR point vector
    /// @param vM magnetic moment vector
    /// @return Magnetic field vector as Eigen 3d Vector
    Vector3d vB(Vector3d vR, Vector3d vM);

    /// @param vR point vector
    /// @param vM magnetic moment vector
    /// @return normalized magnetic field vector
    Vector3d vb(Vector3d vR, Vector3d vM);


    /// @param vB magnetic field vector
    /// @param vK k vector
    /// @return angle between wave vecotr and magnetic field (in radians)
    double theta_kb (Vector3d vB, Vector3d vK);

    /// @param vR point vector
    /// @return Omega_vec corss r_vec 
    Vector3d vBetaR (Vector3d vR);

    /// @brief ExB drift particle velocity component (V = V_|| * b_vec + Udr)
    /// @param vR point vector
    /// @param vM magnetic moment vector
    /// @return drift velocity Eigen 3d Vector  
    Vector3d vUdr (Vector3d vR, Vector3d vM);

    /// @param vR point vector
    /// @param vM magnetic moment vector
    /// @return drift velocity gamma factor
    double gammaU (Vector3d vR, Vector3d vM);

    /// @param vR point vector
    /// @param vM magnetic moment vector
    /// @return distance from the field line footpoint to the polar cap center
    /// @note Only for dipolar magnetic field!!!
    double x_pc(Vector3d vR, Vector3d vM);

    /// @param vR point vector
    /// @param vM magnetic moment vector
    /// @return polar cap angle cooridnate of the magnetic field line footpoint.
    /// Angle is counted from the East-West line on the polar cap
    /// @note Applicable only for magnetic fields with zero torsion
    double phi_pc(Vector3d vR, Vector3d vM);

    /// @param vR point vector
    /// @param vM magnetic moment vector
    /// @param l distance along the ray
    /// @return normalized plasma density
    double gFunc (Vector3d vR, Vector3d vM);

    /// @param l  distance along the ray
    /// @return Plasma density in physical units (g/cm^3)
    double Ne(Vector3d vR, Vector3d vM);

    /// @param vR point vector
    /// @param vM magnetic moment vector
    /// @return Local cyclotron frequency (s^-1)
    double omegaB (Vector3d vR, Vector3d vM);

    /// @param vR point vector
    /// @param vM magnetic moment vector
    /// @return Local plasma frequency (s^-1)
    double omegaP (Vector3d vR, Vector3d vM);


    /// @param vR point vector
    /// @param vK wave vector
    /// @param vM magnetic moment vector
    /// @return wave frequency in plasma rest frame
    double omegaW (Vector3d vR, Vector3d vK, Vector3d vM);



    /// @brief ...
    /// @param l 
    /// @return ...
    bool stop_condition(double l);




    /// @brief polar angle of magnetic filed projected to the plane perpendicular to the ray
    /// @param l distnace along the ray
    /// @return Derivative of BetaB angle in radians along the ray
    double BetaB_derivative (double l);

    /// @param l distance alnog the ray
    /// @return Derivative of Q parameter from Kravtsov-Orlov equations along the ray
    double Q_derivative (double l);

    /// @param l distance along the ray
    /// @return Derivative of Lambda parameter form Kravtsov-Orlov equations along the ray
    double Lambda_derivative (double l);

    /// @param l distance along the ray
    /// @return Derivative of delta parameter form Kravtsov-Orlov equations along the ray
    double delta_derivative(double l);

    /// @param l distance along the ray
    /// @return Derivative of beta+delta parameter form Kravtsov-Orlov equations along the ray
    double BetaBmod_derivative(double l);

    // -----------------------Cyclotron Absorption------------------------

    /// @param l - distance along the ray
    /// @return dtau - differential optical thickness (I = I0 * e^(-tau), tau = integral dtau)
    double dtau (double l);

};

