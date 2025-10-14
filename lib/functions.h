#pragma once

#include <functional>
#include <Eigen/Dense>
using namespace std;
using Eigen::Vector3d;

double ANGLE (Vector3d vec1, Vector3d vec2);
double Arcsinh(double x);
double sgn (double value);
double get_derivative_along_the_ray(std::function<double(double)> func, double l, double dl=1);
double DX(std::function<double(double, double)> func, double x, double y);
double DY(std::function<double(double, double)> func, double x, double y);