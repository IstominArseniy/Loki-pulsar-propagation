#include <map>
#include <vector>
#include <Eigen/Dense>

#include "PSRClass.h"
#include "FixedHeightModelClass.h"

using Eigen::Vector3d;

class CppInterface{
    private:
    std::map<std::string, double> psr_dict;
    std::map<std::string, double> param_dict;
    ObservedRadioPulsar PSR_;
    FixedHeightModel model_;
    
    public:
    CppInterface (std::map<std::string, double> psr_dict, std::map<std::string, double> param_dict); // Consturctor - will be exposed to python
    CppInterface (ObservedRadioPulsar Pulsar, ...); // Consturctor from PSR class
    std::map<std::string, double> find_ILVPA(double phi, int mode); // will be exposed to Python
};



