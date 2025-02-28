#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <iterator>
using namespace std;

#include "constants.h"
#include "read_write.h"

void throw_error(const string msg) {
  cout << msg << endl;
  exit (EXIT_FAILURE);
}
bool is_number(const string& s) {
  if (s.empty())
    return false;
  bool trigger = false;
  if (s[0] == '.') {
    return false;
  }
  for (int i = 0; i < s.size(); i ++) {
    char c = s[i];
    if (c == '.' && !trigger) {
      trigger = !trigger;
      continue;
    } else if (c == '.' && trigger) {
      return false;
    }
    if (i == 0 && c == '-') {
      continue;
    }
    if (!std::isdigit(c))
    return false;
  }
  return true;
}
double read_from_file (string input_name, const string param) {
  ifstream infile(input_name);
  string str;
  stringstream msg;
  if (infile.is_open()) {
    while (getline(infile, str)) {
      // cout << iter << ": " << str << "\n";
      istringstream iss(str);
      vector<string> words{istream_iterator<string>{iss}, istream_iterator<string>{}};
      if(words.size() < 1) continue;
      // cout << words[0] << "\n";
      if (words[0] == param) {
        if(is_number(words[1])) {
          return atof(words[1].c_str());
        } else {
          msg << "ERROR: Converting string to number for parameter '" << param << "' in the input file.";
          throw_error(msg.str());
        }
        break;
      }
    }
  } else {
    throw_error("ERROR: Cannot open the input file.");
  }
  msg << "ERROR: Cannot find the given parameter '" << param << "' in the input file.";
  throw_error(msg.str());
  return 0;
}
string read_from_file_str (string input_name, const string param, const string def_val) {
  ifstream infile(input_name);
  string str;
  stringstream msg;
  if (infile.is_open()) {
    while (getline(infile, str)) {
      istringstream iss(str);
      vector<string> words{istream_iterator<string>{iss}, istream_iterator<string>{}};
      if(words.size() < 1) continue;
      if (words[0] == param) {
        if(!is_number(words[1])) {
          return words[1].c_str();
        } else {
          msg << "ERROR: param '" << param << "' is not a string in the input file.";
          throw_error(msg.str());
        }
        break;
      }
    }
  } else {
    throw_error("ERROR: Cannot open the input file.");
  }
  return def_val;
}

void read_in_out(string &in, string &out, int argc, char* argv[]) {
  bool found_input = false, found_output = false;
  for (int i=0; i<argc; ++i) {
		if (typeid(argv[i]) == typeid(char*)) {
			if (strncmp(argv[i], "-i", 3) == 0) {
				if (i+1 < argc) {
					in = argv[i+1];
          found_input = true;
        }
				else {
          throw_error("ERROR: No input file given.");
				}
			}
      if (strncmp(argv[i], "-o", 3) == 0) {
				if (i+1 < argc) {
					out = argv[i+1];
          found_output = true;
        }
				else {
          throw_error("ERROR: No output path given.");
				}
			}
		}
	}
  if (!found_input) {
    throw_error("ERROR: No input file given.");
  } 
  // else {
  //   cout << "INPUT: " << in << "\n";
  // }
  if (!found_output) {
    out = "output";
  }
}

pair<double, double> find_thread_phases(double phi_start_global, double phi_end_global, double phi_step, int size, int rank){
  /*
  function to find whitch phases should be processed by this thread 
  returns phi_initial and phi_final for thread
  */
  pair <double, double> phases;
  int Nsteps_global = int((phi_end_global - phi_start_global) / phi_step);
  int Nsteps = Nsteps_global / size;
  int residual = Nsteps_global % size;
  if(rank < residual){
    phases.first = phi_start_global + rank * (Nsteps + 1) * phi_step;
    phases.second = phi_start_global + (rank + 1) * (Nsteps + 1) * phi_step;
  }
  else{
    phases.first = phi_start_global + ((Nsteps + 1) * residual + Nsteps * (rank - residual)) * phi_step;
    phases.second = phi_start_global + ((Nsteps + 1) * residual + Nsteps * (rank - residual + 1)) * phi_step;
  }   
  return phases;
}