#include <Eigen/Dense>
#include <map>

#pragma once

using Eigen::Vector3d;

class ObservedRadioPulsar{
    public:
    ObservedRadioPulsar();
    ObservedRadioPulsar(std::map<std::string, double> psr_dict); // constructor from dictionary
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
    Vector3d observer_vec;
    Vector3d Omega_vec;
};