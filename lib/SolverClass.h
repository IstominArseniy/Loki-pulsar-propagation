#include <PSRClass.h>
#include <FixedHeightModelClass.h>

class FixedHeightSolver{
    private:
    double phi_;
    double mode_;
    double PSR_;
    FixedHeightModel model_;

    public:
    FixedHeightSolver(double phi, int mode, ObservedRadioPulsar &PSR, FixedHeightModel &model);
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
}