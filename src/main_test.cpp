#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <ctime>

#include "../lib/PSRClass.h"
#include "../lib/FixedHeightModelClass.h"
#include "../lib/SolverClass.h"
#include "../lib/CppInterfaceClass.h"
#include "../lib/read_write.h"
#include "../lib/constants.h"

int main(int argc, char* argv[]) {
    //create interface for calculation
    CppInterface calculation_interface;
    //reading input output paths
    std::string input_path, output_path;
    read_in_out(input_path, output_path, argc, argv); // place input and output paths from cmd -i -o flags into input_path and output_path variables
    calculation_interface.init_from_file(input_path);
    //read phase start, end ans step
    double phi_start = read_from_file(input_path, "phi_start");
    double phi_end = read_from_file(input_path, "phi_end");
    double phi_step = read_from_file(input_path, "phi_step");
    //read mode to calculate
    int emission_mode = read_from_file(input_path, "mode");
    //important info
    std::time_t t_start = std::time(nullptr);
    std::cout << "Calculation started!" << std::endl;
    std::cout << "R_esc = " << calculation_interface.get_R_escape() << std::endl;
    std::cout << "R_LC = " << calculation_interface.get_RLC() << std::endl;
    FixedHeightSolver solver(0, 0, calculation_interface.PSR_, calculation_interface.model_);
    std::vector<double> theta_init = {solver.BetaB(0) + constants::PI/3, 0};
    std::vector<double> theta_final;
    theta_final = solver.solve_KO_equations(theta_init, 0, 1000, "./");
    std::cout << "theta1 = " << theta_final[0] << " " << " theta2 = " << theta_final[1] << std::endl;
    std::time_t t_end = std::time(nullptr);
    std::cout << "Ellapsed time = " << t_end - t_start << std::endl;
    return 0;
}