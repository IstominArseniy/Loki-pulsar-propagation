#include <math.h>
#include <vector>
#include <iostream>
#include "RHS.h"
#include "constants.h"
#include "read_write.h"
#include "process_functions.h"
#include "initialize.h"

void RHS(double R, double *f, double *dydx) {
	double coeff = constants::R_star * Globals::omega / (2.0 * constants::c);

	double LL = Lambda (R);
	double QQ = Q (R);
	double BB = BetaB (R);
	double DD = delta (R);

	dydx[0] = coeff * (-LL / QQ - LL * cos(2 * f[0] - 2 * BB - 2 * DD) * sinh(2 * f[1]));
	dydx[1] = coeff * LL * sin(2 * f[0] - 2 * BB - 2 * DD) * cosh(2 * f[1]);
}

void RHS_for_boost(const std::vector<double>& f, std::vector<double> &dydx, double R) {
	double coeff = constants::R_star * Globals::omega / (2.0 * constants::c);

	double LL = Lambda (R);
	double QQ = Q (R);
	double BB = BetaB (R);
	double DD = delta (R);

	dydx[0] = coeff * (-LL / QQ - LL * cos(2 * f[0] - 2 * BB - 2 * DD) * sinh(2 * f[1]));
	dydx[1] = coeff * LL * sin(2 * f[0] - 2 * BB - 2 * DD) * cosh(2 * f[1]);
}


// double* myRHS(double R, double *f, double) {
// 	double coeff = constants::R_star * Globals::omega / (2.0 * constants::c);
// 	double LL = Lambda (R);
// 	double QQ = Q (R);
// 	double BB = BetaB (R);
// 	double DD = delta (R);

// 	dydx[0] = coeff * (-LL / QQ - LL * cos(2 * f[0] - 2 * BB - 2 * DD) * sinh(2 * f[1]));
// 	dydx[1] = coeff * LL * sin(2 * f[0] - 2 * BB - 2 * DD) * cosh(2 * f[1]);
// }
