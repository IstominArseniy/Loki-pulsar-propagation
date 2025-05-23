/*
command to run this programm 
mpirun --use-hwthread-cpus -np 12 ./bin/loki -i bin/loki.input -o bin/my_output
*/
 

#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <boost/numeric/odeint.hpp>
#include "mpi.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <filesystem>

#include "../lib/constants.h"
#include "../lib/read_write.h"

#include "../lib/functions.h"
#include "../lib/process_functions.h"
#include "../lib/initialize.h"

#include "../lib/integrator.h"
#include "../lib/RHS.h"

#include "../lib/diffeqsolver.h"

namespace fs = std::filesystem;
using namespace std;


int main(int argc, char* argv[]) {
  using namespace boost::numeric::odeint;
  initialize(argc, argv);
  //CREATE FOLDER FOR CALCULATION IF NOT EXIST-----------------------
  std::string global_data_path = Globals::out_path+"/"+Globals::RUN_ID+"_global_data";
  if (!fs::is_directory(global_data_path) || !fs::exists(global_data_path)) { // Check if folder exists
    fs::create_directory(global_data_path); // create folder
  } 
  //-----------------------------------------------------------------

  double phi_start_global = read_from_file(Globals::input_name, "phi_start");
  double phi_end_global = read_from_file(Globals::input_name, "phi_end");
  double phi_t_step = read_from_file(Globals::input_name, "phi_step");
  // MPI INICIALISATION-----------------------------------------------
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
  //-------------------------------------------------------------------
  //PHASE BOARDERS COMPUTATION-----------------------------------------
  double phi_t_start, phi_t_end;
  auto phases = find_thread_phases(phi_start_global, phi_end_global, phi_t_step, size, rank);
  phi_t_start = phases.first;
  phi_t_end = phases.second;
  //--------------------------------------------------------------------
  //INITIAL INFORMATION-------------------------------------------------
  if(rank==0){
    cout << "RUN_ID: " << Globals::RUN_ID << "\n\n";
    cout << "\nR_esc: " << Globals::RESCAPE << endl;
    cout << "R_A: " << Globals::ROMODE << endl;
    cout << "R_lc: " << Globals::RLC << endl << endl;
  }
  //--------------------------------------------------------------------

  // SIMULATION STARTS HERE />
  for (double phi_t = phi_t_start; phi_t < phi_t_end; phi_t += phi_t_step) { // Phase switch
    //INITIALISATION---------------------------------------------------------------------------------
    Globals::PHI0 = phi_t * constants::PI / 180.0;
    ofstream output(global_data_path + "/" + Globals::RUN_ID + "_" + to_string(Globals::PHI0 * 180 / constants::PI) + ".dat");
    setInitPoints();
    double x1, x2;
    std::vector<double> dep_vars(2);
    // Avoiding initial osc. region
    // x1 = 10;
    x1 = find_initial_point(false);
    x2 = std::min(1.5 * Globals::RESCAPE, Globals::RLC - 1.1 * Globals::R_em);
    x2 = 1.5 * Globals::RESCAPE;
    // Initial values />
    dep_vars[0] = approximate_solution_theta0(x1, Globals::mode);
    dep_vars[1] = approximate_solution_theta1(x1, Globals::mode);
    // std::cout << phi_t << " " << x1 << " " << dep_vars[0] << " " << dep_vars[1] << std::endl; 
    // </ Initial values
    double PA = dep_vars[0] * 180 / constants::PI;
    double tau = constants::PI * constants::R_star * integrate(dtau, x1, Globals::RLC - 1.1 * Globals::R_em) / (constants::c * Globals::omega);
    double II0 = std::pow(gFunc(0), 1);
    double II = II0 * exp (-tau);
    // double II = II0;
    double PA0_rad = dep_vars[0];
    double VV = II * tanh(2.0 * dep_vars[1]);
    output << phi_t << " " << II0 << " " << VV << " " << PA << std::endl;
    //---------------------------------------------------------------------------------------------------
    //INTEGRATION----------------------------------------------------------------------------------------
    // std::cout << phi_t << " " << r_perp(0) << " " <<  phi_pc(0) << " " << gFunc(0) <<  std::endl;
    double eps_abs = 1e-6, eps_rel = 1e-6, h_init = 1.0e-3;
    auto addaptive_stepper = make_controlled(eps_abs, eps_rel, runge_kutta_dopri5 < std::vector<double> >());
    std::string path = Globals::out_path + "/" + Globals::RUN_ID + "_theta_data/" + Globals::RUN_ID + "_theta_"+ to_string(Globals::PHI0 * 180 / constants::PI) + ".dat";
    ofstream plot(path);
    auto ode_range = make_adaptive_time_range(addaptive_stepper, RHS_for_boost, dep_vars, x1, x2, h_init);
    auto it=ode_range.first;
    double last_step = x1;
    std::vector<double> dep_vars_prev_step(2);
    double cr_step_len = h_init;
    while(it != ode_range.second){
      if(stop_condition(it->second)){
        break;
      }
      plot << it->second <<  ", " << it->first[0] << ", " << it->first[1] << ", " << constants::R_star * Globals::omega / (2.0 * constants::c) * Lambda(it->second) << ", " << BetaB (it->second) + delta (it->second) << std::endl;
      // std::cout << it->second << " " << it->second + 5.1 * cr_step_len << " " << stop_condition(it->second + 5.1 * cr_step_len) << std::endl;
      last_step = it->second;
      dep_vars_prev_step[0] = it->first[0];
      dep_vars_prev_step[1] = it->first[1];
      it++;
      cr_step_len = it->second - last_step;
    }
    // integrate_adaptive(addaptive_stepper, RHS_for_boost, dep_vars, x1, x2, h_init, Observer(plot));
    //double RM = integrate(RM_dencity, Globals::R_em, Globals::RLC); // RM calculation
    //---------------------------------------------------------------------------------------------------
    //OUTPUTS--------------------------------------------------------------------------------------------
    VV = II * tanh(2.0 * dep_vars_prev_step[1]);
    PA = dep_vars_prev_step[0] * 180.0 / constants::PI;
    // cout << "\tI: " << II << "\n\tV: " << VV << "\n\tPA: " << -PA << endl << endl;
    output << phi_t << " " << II << " " << VV << " " << PA << std::endl ;
    output << (dep_vars[0] - PA0_rad) / std::pow((constants::c / Globals::freqGHz / 1e9), 2) << std::endl;
    output << dep_vars[0] << std::endl;
    output.close();
    //--------------------------------------------------------------------------------------------------
  }
  //DATA REARANGEMENT---------------------------------------------
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank == 0){ // should be done only by one thread
    ofstream output0(Globals::out_path + "/" + Globals::RUN_ID + "_0.dat");
    ofstream output1(Globals::out_path + "/" + Globals::RUN_ID + "_1.dat");
    /*
      POSSIBLE ADDITIONAL OUTPUTS
      ofstream output3(Globals::out_path + "/" + Globals::RUN_ID + "_RMs.dat");
      ofstream output4(Globals::out_path + "/" + Globals::RUN_ID + "PAs.dat");
    */
    for (const auto & entry : fs::directory_iterator(global_data_path)){
        ifstream in(entry.path());
        double phase, I0, V0, PA0, I1, V1, PA1;
        in >> phase >> I0 >> V0 >> PA0;
        in >> phase >> I1 >> V1 >> PA1;
        in.close();
        output0 << phase << " " << I0 << " " << V0 << " " << PA0 << std::endl;
        output1 << phase << " " << I1 << " " << V1 << " " << PA1 << std::endl;
    }
    output0.close();
    output1.close();
  }
  //-------------------------------------------------------------
  MPI_Finalize();
	return 0;
}
