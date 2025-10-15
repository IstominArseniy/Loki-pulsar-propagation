#include "PSRClass.h"
#include "constants.h"

using Eigen::Vector3d;

ObservedRadioPulsar::ObservedRadioPulsar(){}

ObservedRadioPulsar::ObservedRadioPulsar(std::map<std::string, double> psr_dict)
{
    B12 = psr_dict["B12"];
    freqGHz = psr_dict["freqGHz"];
    omega_obs = 2 * constants::PI * freqGHz * 1e9;
    Period = psr_dict["Period"];
    Omega = 2 * constants::PI / Period;
    chi_deg = psr_dict["chi_deg"];
    beta_deg = psr_dict["beta_deg"];
    chi = chi_deg * constants::PI / 180;
    beta = beta_deg * constants::PI / 180;
    dzeta = chi - beta;
    Rs = psr_dict["Rs"];
    RLC = (constants::c / Omega) / Rs;

    Omega_vec << 0, 0, Omega;
    observer_vec << std::sin(dzeta), 0, std::cos(dzeta);
}
