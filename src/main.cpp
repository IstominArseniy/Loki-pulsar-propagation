#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <map>
#include <fstream>

#include "../lib/PSRClass.h"
#include "../lib/FixedHeightModelClass.h"
#include "../lib/SolverClass.h"
#include "../lib/CppInterfaceClass.h"
#include "../lib/read_write.h"

int main(int argc, char* argv[]) {
    //create interface for calculation
    CppInterface calculation_interface;
    //reading input output paths
    std::string input_path, output_path;
    read_in_out(input_path, output_path, argc, argv); // place input and output paths from cmd -i -o flags into input_path and output_path variables
    calculation_interface.init_from_file(input_path);
    make_dir(output_path);
    std::ofstream output(output_path+"/" +"output.dat");
    //read phase start, end ans step
    double phi_start = read_from_file(input_path, "phi_start");
    double phi_end = read_from_file(input_path, "phi_end");
    double phi_step = read_from_file(input_path, "phi_step");
    //read mode to calculate
    int emission_mode = read_from_file(input_path, "mode");
    // computation of the profile
    auto result = calculation_interface.calculate_profile(phi_start, phi_end, phi_step, emission_mode, false);
    //writing data into output
    for(auto& elem : result){
        output << elem["phase"] << " " << elem["I"] << " " << elem["V"] << " " << elem["PA"] << std::endl;
    }
    return 0;
}
