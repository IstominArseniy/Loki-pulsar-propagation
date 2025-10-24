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
    CppInterface calculation_interface;
    std::string input_path, output_path;
    read_in_out(input_path, output_path, argc, argv); // place input and output paths from cmd -i -o flags into input_path and output_path variables
    calculation_interface.init_from_file(input_path);
    make_dir(output_path);
    std::ofstream output(output_path+"/" +"output.dat");
    auto result = calculation_interface.calculate_profile(-15, 15, 1, 0, false);
    for(auto& elem : result){
        output << elem["phase"] << " " << elem["I"] << " " << elem["V"] << " " << elem["PA"] << std::endl;
    }
    return 0;
}
