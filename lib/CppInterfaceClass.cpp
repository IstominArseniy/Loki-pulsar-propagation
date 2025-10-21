#include <map>
#include <vector>
#include <cmath>
#include <iostream>
#include "constants.h"
#include "functions.h"
#include "CppInterfaceClass.h"
#include "SolverClass.h"


using Eigen::Vector3d;

CppInterface::CppInterface(std::map<std::string, double> psr_dict, std::map<std::string, double> param_dict){ // Constructor form dict
    psr_dict_ = psr_dict;
    param_dict_ = param_dict;
    PSR_ = ObservedRadioPulsar(psr_dict_);
    model_ = FixedHeightModel(param_dict_, PSR_);
}

std::map<std::string, double> CppInterface::find_ILVPA(double phi, int mode)
{
    FixedHeightSolver solver(phi, mode, PSR_, model_);
    double l1 = solver.find_initial_point(false);
    // double l2 = PSR_.RLC; // Probably redo as solver.find_final_point()
    double l2 = 3 * get_R_escape();
    std::vector<double> theta_init = solver.find_approximate_KO_solution(l1);
    std::vector<double> theta_final = solver.solve_KO_equations(theta_init, l1, l2);
    double I = solver.find_intensity();
    double tau = solver.get_tau();
    I *= std::exp(-tau);
    double V = I * std::tanh(2.0 * theta_final[1]);
    double PA = theta_final[0] * 180.0 / constants::PI;
    std::map<std::string, double> result;
    result["I"] = I;
    result["L"] = std::sqrt(I*I - V*V);
    result["V"] = V;
    result["PA"] = PA;
    return result;
}

std::vector<std::map<std::string, double> > CppInterface::calculate_profile(double phi1, double phi2, double phi_step, int mode, bool use_multiprocessing, double n_threads){
    std::vector<std::map<std::string, double> > profile;
    if (use_multiprocessing == false){
        double phi_tmp = phi1;
        while(phi_tmp <= phi2){
            auto result = find_ILVPA(phi_tmp, mode);
            profile.push_back(result);
            phi_tmp += phi_step;
        }
    }
    else{
        std::cout << "NOT IMPLEMENTED YET!!!";
    }
    return profile;
}

double CppInterface::get_R_escape()
{
    return 1.0e3 * std::pow(model_.lambda / 1.0e4, 1.0/3.0) * std::pow(model_.gamma0 / 100.0, -6.0/5.0) * std::pow(PSR_.B12, 2.0/5.0) * std::pow(PSR_.freqGHz, -2.0/5.0) * std::pow(PSR_.Period, -1.0/5.0);

}
