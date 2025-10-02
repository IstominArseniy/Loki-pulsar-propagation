#include <map>
#include <vector>
#include <Eigen/Dense>
using Eigen::Vector3d;

class CppFixedHeightModel{
    private:
    std::map<std::string, double> psr_dict;
    std::map<std::string, double> param_dict;
    double fr;
    double fphi;

    public:
    CppFixedHeightModel (std::map<std::string, double> psr_dict, std::map<std::string, double> param_dict); // Consturctor - will be exposed to python
    std::map<std::string, double> find_ILVPA(double phi, int mode); // will be exposed to Python
    Vector3d Bfield(Vector3d vR, Vector3d m); // Magentic field model
};


class FixedHeightCalculator{
    private:
    //-------PSR PARMAS-------------
        double B12;
        double freqGhz;
        double omega_obs;
        double Period;
        double Omega;
        double chi_deg;
        double beta_deg;
        double chi;
        double beta;
        double dzeta;
        double Rs;
        double RLC;
    //-------------------------------
    //-------MODEL PARAMS------------
        double fr;
        double fphi;
        double Rem;
        double lambda;
        double gamma0;
        double L_SHIFT;
    //-------------------------------
    //-------PHASE PARAMS--------------
        double phi;
        int mode;
        double theta_em;
        double phi_em;
        Vector3d observer_vec;
    //-------------------------------

    public:
        FixedHeightCalculator(double phi, int mode, std::map<std::string, double> psr_dict, std::map<std::string, double> param_dict);
        // add Bfiled model and g(x, phi) 

        /// @brief solve Kravtsov-Orlov equations
        /// @param phi pulsar phase
        /// @param theta_inital inital values for comlex polarisation angle (as length 2 vector)
        /// @param mode emission mode (1 - Omode, 0 - Xmode)
        /// @return theta_final -  final complex polarisation angle (as length 2 vector)
        std::vector<double> solve_KO_equations(std::vector<double> theta_initial, double l1, double l2);

        /// @brief find approximate solution for Kravtsov-Orlov equations in dense plasma regions 
        /// @param theta_inital inital values for comlex polarisation angle (as length 2 vector)
        /// @param l1 - initial point
        /// @param l2 - final point
        /// @param mode emission mode (1 - Omode, 0 - Xmode)
        /// @return theta_final -  final complex polarisation angle (as length 2 vector)
        std::vector<double> find_approximate_KO_solution(std::vector<double> theta_initial);

        /// @brief find approximate solution for Kravtsov-Orlov equations in dense plasma regions 
        /// @param l - distance from the star. Theta_inintal is determined by emission mode
        /// @param mode emission mode (1 - Omode, 0 - Xmode)
        /// @return theta_final -  final complex polarisation angle (as length 2 vector)
        std::vector<double> find_approximate_KO_solution(double l);

        double find_initial_point();

    private:
    
        /**
         * @brief Magnetic moment vector unit vector
         * @param l distance alnog the ray
         * @return Vector3d magnetic moment
         */
        Vector3d vMoment (double l);
};


