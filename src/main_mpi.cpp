/*
command to run this programm 
mpirun --use-hwthread-cpus -np 10 ./bin/loki_mpi -i inputs_outputs/loki.input -o inputs_outputs/my_output
*/
 


#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <map>
#include <fstream>

#include "mpi.h"

#include "../lib/PSRClass.h"
#include "../lib/FixedHeightModelClass.h"
#include "../lib/SolverClass.h"
#include "../lib/CppInterfaceClass.h"
#include "../lib/read_write.h"
#include "../lib/functions.h"

int main(int argc, char* argv[]) {
    //reading input output paths
    std::string input_path, output_path;
    read_in_out(input_path, output_path, argc, argv); // place input and output paths from cmd -i -o flags into input_path and output_path variables
    //create interface for calculation
    make_dir(output_path);
    make_dir(output_path+"/logs");
    CppInterface calculation_interface(output_path+"/logs");
    calculation_interface.init_from_file(input_path);
    //read phase start, end ans step
    double phi_start_global = read_from_file(input_path, "phi_start");
    double phi_end_global = read_from_file(input_path, "phi_end");
    double phi_step = read_from_file(input_path, "phi_step");
    //read mode to calculate
    int emission_mode = read_from_file(input_path, "mode");
    //init MPI
    int rank, size;
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (MPI_Comm_size(MPI_COMM_WORLD, &size) != MPI_SUCCESS) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (rank == 0){
        std::cout << "Calculation started with " << size << " processes." << std::endl;
        std::cout << "R_esc = " << calculation_interface.get_R_escape() << std::endl;
        std::cout << "R_LC = " << calculation_interface.get_RLC() << std::endl;
        std::cout << "Half-openning angle = " << calculation_interface.get_rho() << std::endl;
        if (emission_mode == 0){
            std::cout << "Emission mode = X " << std::endl;
        }
        else{
            std::cout << "Emission mode = O" << std::endl;
        }
    }
    //split computation domain
    auto phases = split_phases(phi_start_global, phi_end_global, phi_step, size, rank);
    // computation of the profile
    
    auto result = calculation_interface.calculate_profile(phases.first, phases.second, phi_step, emission_mode, false);
    auto init_data = calculation_interface.calculate_initial_profile(phases.first, phases.second, phi_step, emission_mode);
    //wrting output data to file
    MPI_Barrier(MPI_COMM_WORLD); // wait untill all processes have finished
    if(rank == 0){ // recieve data from all processes
        std::ofstream output(output_path+"/" +"output.dat"); // output file
        std::ofstream output_init(output_path+"/" +"output_init.dat"); // output file for initial data
        for(auto& elem : result){ // output from rank 0 process
            output << elem["phase"] << " " << elem["I"] << " " << elem["V"] << " " << elem["PA"] << std::endl;
        }
        for(auto& elem : init_data){ // output from rank 0 process (for initial profile)
            output_init << elem["phase"] << " " << elem["I"] << " " << elem["V"] << " " << elem["PA"] << std::endl;
        }
        for(int rank_id=1; rank_id<size; rank_id++){
            int N_counts;
            MPI_Recv(&N_counts, 1, MPI_INT, rank_id, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // recieve number of phase counts
            std::vector<double> vectorized_data(N_counts*5);
            MPI_Recv(vectorized_data.data(), N_counts*5, MPI_DOUBLE, rank_id, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // recieve data
            std::vector<double> vectorized_init_data(N_counts*5);
            MPI_Recv(vectorized_init_data.data(), N_counts*5, MPI_DOUBLE, rank_id, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // recieve init data
            auto recived_result = CppInterface::restore_data(vectorized_data);
            auto recived_init_data = CppInterface::restore_data(vectorized_init_data);
            for(auto& elem : recived_result){
                output << elem["phase"] << " " << elem["I"] << " " << elem["V"] << " " << elem["PA"] << std::endl;
            }
            for(auto& elem : recived_init_data){ 
                output_init << elem["phase"] << " " << elem["I"] << " " << elem["V"] << " " << elem["PA"] << std::endl;
            }
        }
        output.close();
    }
    else{ // send data to proecss with rank 0
        int N_counts = result.size();
        MPI_Ssend(&N_counts, 1, MPI_INT, 0, 0, MPI_COMM_WORLD); // send number of phase counts
        auto vectorized_data = CppInterface::vectorize_data(result);
        MPI_Ssend(vectorized_data.data(), N_counts*5, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); // send data
        auto vectorized_init_data = CppInterface::vectorize_data(init_data);
        MPI_Ssend(vectorized_init_data.data(), N_counts*5, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); // send initial data
    }
    MPI_Finalize();
    return 0;
}
