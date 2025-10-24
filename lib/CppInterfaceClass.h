#include <map>
#include <vector>
#include <Eigen/Dense>

#include "PSRClass.h"
#include "FixedHeightModelClass.h"


using Eigen::Vector3d;

class CppInterface{
    private:
    ObservedRadioPulsar PSR_;
    FixedHeightModel model_;
    
    public:
    CppInterface ();
    CppInterface (std::map<std::string, double> psr_dict, std::map<std::string, double> param_dict); // Consturctor - will be exposed to python
    CppInterface (ObservedRadioPulsar PSR, std::map<std::string, double> param_dict); 
    void init_from_file(std::string filename);

    // CppInterface (ObservedRadioPulsar Pulsar, ...); // Consturctor from PSR class
    std::map<std::string, double> find_ILVPA(double phi, int mode); // will be exposed to Python
    double get_R_escape();
    std::vector<std::map<std::string, double> > calculate_profile(double phi1, double phi2, double phi_step, int mode, bool use_multiprocessing=true, double n_threads=12);
};



