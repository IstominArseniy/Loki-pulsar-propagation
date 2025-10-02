#include <map>
#include <Eigen/Dense>
#include "PSRClass.h"

using Eigen::Vector3d;

class FixedHeightModel{
    private:
    double Rm;
    double fr;
    double fphi;
    double L_SHIFT;
    double lambda;
    double gamma0;
    ObservedRadioPulsar PSR_;

    public:
    FixedHeightModel(std::map<std::string, double> model_dict, ObservedRadioPulsar PSR);
    double density_profile(double x_pc, double phi_pc);
    Vector3d Bfield(Vector3d vR, Vector3d m);

    /// @brief Particle energy distribution function (both for e- and e+)
    /// @param gamma 
    /// @return f(gamma) | dN = f(gamma) d gamma
    double fDist (double gamma);
    
    /**
     * @brief Function to calculate <1/gamma^3 * 1/(1 - gamma^2 A^2) >(A) average, which is present in Lambda epxression for Kravtosov-Orlov equations
     * 
     * fixed fDist Assumed  - TODO make direct conversion from fDist to this function
     * @param A - variable
     * @return double average
     */
    double Lambda_type_avrg(double A);

    /**
     * @brief Function to calculate <1/(1 - gamma^2 A^2) >(A) average, which is present in ImEps epxression for Kravtosov-Orlov equations
     * 
     * fixed fDist Assumed  - TODO make direct conversion from fDist to this function
     * @param A - variable
     * @return double average
     */
    double ImEps_type_avrg(double A);

    double Q_type_avrg(double A);

};