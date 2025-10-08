#include <map>
#include <vector>
#include <cmath>
#include "constants.h"
#include "functions.h"
#include "CppInterfaceClass.h"

using Eigen::Vector3d;

CppInterface::CppInterface(std::map<std::string, double> psr_dict, std::map<std::string, double> param_dict){ // Constructor form dict
    psr_dict_ = psr_dict;
    param_dict_ = param_dict;
    PSR_ = PSR(psr_dict);
    model_ = FixedHeightModel(psr_dict);
}

std::map<std::string, double> CppInterface::find_ILVPA(double phi, int mode)
{
    solver = FixedHeightSolver(phi, mode, psr_dict_, param_dict_);
    double l1 = solver.find_initial_point();
    double l2 = PSR_.RLC; // Probably redo as solver.find_final_point()
    std::vector<double> theta_init = solver.find_approximate_KO_solution(l1);
    std::vector<double> theta_final = solver.solve_KO_equations(theta_init, l1, l2);
    double I = solver.find_intensity();
    double V = I * std::tanh(2.0 * theta_final[1]);
    double PA = theta_final[0] * 180.0 / constants::PI;
    std::map<std::string, double> > result;
    result['I'] = I;
    result['L'] = std::sqrt(I**2 - V**2);
    result['V'] = V;
    result['PA'] = PA;
    return result;
}


