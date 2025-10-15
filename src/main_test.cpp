#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <map>
// #include "../lib/CppInterfaceCalss.h"
#include "../lib/PSRClass.h"
#include "../lib/FixedHeightModelClass.h"
#include "../lib/SolverClass.h"
#include "../lib/CppInterfaceClass.h"

int main(int argc, char* argv[]) {
    std::map<std::string, double> psr_dict{
        {"B12", 1}, {"Period", 1}, {"freqGHz", 1}, {"chi_deg", 45}, {"beta_deg", 2}, {"Rs", 1.2e6}
    };
    ObservedRadioPulsar PSR(psr_dict);
    std::map<std::string, double> model_dict{
        {"fr", 1}, {"fphi", 1}, {"Rem", 50}, {"L_SHIFT", 1}, {"lambda", 1000}, {"gamma0", 100}
    };
    FixedHeightModel model(model_dict, PSR);
    CppInterface interface (psr_dict, model_dict);
    std::map<std::string, double> answer = interface.find_ILVPA(-5, 0);
    std::cout << answer["V"] << std::endl;
    return 0;
}
