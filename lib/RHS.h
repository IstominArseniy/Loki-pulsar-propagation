#include<vector>
#pragma RHS

void RHS(double R, double *f, double *dydx);
void RHS_for_boost(const std::vector<double>& f, std::vector<double> &dydx, double R);