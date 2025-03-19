#pragma FUNCTIONS
#include <Eigen/Dense>
using namespace std;
using Eigen::Vector3d;


double SCALAR (vector <double> vec1, vector <double> vec2);
vector <double> CROSS (vector <double> vec1, vector <double> vec2);
double NORM (vector <double> vec);
vector <double> SUM (vector <double> vec1, vector <double> vec2);
vector <double> TIMES (double a, vector <double> vec);
vector <double> NORMALIZE (vector <double> vec);
double ANGLE (Vector3d vec1, Vector3d vec2);
double Arcsinh(double x);