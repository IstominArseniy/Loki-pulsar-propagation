#include <vector>
#include <math.h>
#include <Eigen/Dense>
#include <functional>
#include "functions.h"

using Eigen::Vector3d;

double ANGLE (Vector3d vec1, Vector3d vec2) {
  return 2 * std::atan2((vec2.normalized() - vec1.normalized()).norm(), (vec2.normalized() + vec1.normalized()).norm());
} // Angle

double Arcsinh(double x) {
  return std::log(x + std::sqrt(std::pow(x, 2.0) + 1.0));
}

double sgn (double value) {
  if (value >= 0.0) {
    return 1.0;
  } else {
    return -1.0;
  }
}

double get_derivative_along_the_ray(std::function<double(double)> func, double l, double dl)
{
  return (func(l + dl) - func(l)) / dl; // Very approximate, but high accuracy is not needed (probably better still to make at least 2nd order)
}

double DX(std::function<double(double, double)> func, double x, double y) {
    double h = 0.00001;
    double fm2 = func(x - 2 * h, y);
    double fp2 = func(x + 2 * h, y);
    double fm1 = func(x - h, y);
    double fp1 = func(x + h, y);
    return (fm2 - 8 * fm1 + 8 * fp1 - fp2) / (12 * h);
}
double DY (std::function<double(double, double)> func, double x, double y) {
    double h = 0.00001;
    double fm2 = func(x, y - 2 * h);
    double fp2 = func(x, y + 2 * h);
    double fm1 = func(x, y - h);
    double fp1 = func(x, y + h);
    return (fm2 - 8 * fm1 + 8 * fp1 - fp2) / (12 * h);
}