#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <map>
// #include "../lib/CppInterfaceCalss.h"
#include "../lib/PSRClass.h"

int main(int argc, char* argv[]) {
    std::map<std::string, double> psr_dict{
        {"B12", 1}, {"Period", 1}, {"freqGhz", 1}, {"chi_deg", 45}, {"beta_dega", 2}, {"Rs", 1.2e6}
    };
    ObservedRadioPulsar PSR(psr_dict);
    return 0;
}
