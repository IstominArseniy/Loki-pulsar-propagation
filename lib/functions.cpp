#include <vector>
#include <math.h>
#include <Eigen/Dense>
#include "functions.h"
using namespace std;

using Eigen::Vector3d;

double ANGLE (Vector3d vec1, Vector3d vec2) {
  return 2 * atan2((vec2.normalized() - vec1.normalized()).norm(), (vec2.normalized() + vec1.normalized()).norm());
} // Angle

double Arcsinh(double x) {
  return log(x + sqrt(pow(x, 2.0) + 1.0));
}
