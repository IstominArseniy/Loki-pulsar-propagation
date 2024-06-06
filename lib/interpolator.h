#include <functional>
#include <vector>
#include <string>

double bilinear_interp(std::vector<double>& xs, std::vector<double>& ys, std::vector<std::vector<double> >& data, double x, double y);


class interpolator2D{
private:
        int Nrs, Nphis;
        std::vector<double> rs;
        std::vector<double> phis;
        std::vector<std::vector<double> > fs;
public:
        interpolator2D ();
        void init(std::string fname);
        double get_f(double r, double phi);
};

void rk4_adaptive_odeint(
	double ystart[],
	int nvar,
	double x1,
	double x2,
	double eps,
	double h1,
	double hmin,
	bool save_trajectory, std::vector<std::vector<double> >& trajectory,
	std::function<void(double, double *, double *)> derivs,
	std::function<bool(double, double *)> early_stop_condition
	);