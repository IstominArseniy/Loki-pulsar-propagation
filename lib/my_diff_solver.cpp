/*#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include "diffeqsolver.h"
#include "NRutil.h"
#include "process_functions.h"
#include "initialize.h"
#include "constants.h"
using namespace std;

vector<double> rk4(vector<double> y, vector<double> dydx, double x, double h, vector<double> (*derivs)(double, vector<double> )) {
    int i;
    double xh,hh,h6;
    int n = y.size();

    vector<double> dym(n);
    vector<double> dyt(n);
    vector<double> yt(n);
    vector<double> yout(n);

    hh=h*0.5;
    h6=h/6.0;
    xh=x+hh;
    for (i=0;i<n;i++)
        yt[i]=y[i]+hh*dydx[i];
    dyt = (*derivs)(xh, yt);
    for (i=0;i<n;i++)
        yt[i]=y[i]+hh*dyt[i];
    dym = (*derivs)(xh, yt);
    for (i=0;i<n;i++) {
        yt[i]=y[i]+h*dym[i];
        dym[i] += dyt[i];
    }
    dyt = (*derivs)(x+h, yt);
    for (i=0;i<n;i++)
        yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
    return yout;
}

#define PGROW -0.20
#define PSHRNK -0.25
#define FCOR 0.06666666
#define SAFETY 0.9
#define ERRCON 1.89e-4
#define HMAX 10.0

#define MAXSTP 10000 // normally 10000
#define TINY 1.0e-30

void rkqc(
          vector<double> y,
          vector<double> dydx,
          double x,
          double htry,
          double eps,
          vector<double> yscal,
          double hdid,
          double hnext,
          vector<double> (*derivs)(double, vector<double> )
          ) {
    int i;
    int n = y.size();
    double x_new,hh,h,temp,errmax;

    vector<double> y1(n);
    vector<double> y2(n);
    vector<double> dydx1(n);

    for (i=0;i<n;i++) {
        ysav[i]=y[i];
        dysav[i]=dydx[i];
    }
    h=htry;
    for (;;) {
        // calculation with a two halfsteps
        hh=0.5*h;
        y1 = rk4(y, dydx, n, x, hh, derivs);
        dydx1 = derivs(n, y1);
        y_temp = rk4(н1, dydx1, n, x+hh, hh, derivs);
        x_new = x+h;
        if (x_new == x) {
            cout << "Step size too small in routine RKQC\n";
            break;
        }
        // calculation with a one whole step
        y2 = rk4(y, dydx, x, h, derivs);

        errmax=0.0;
        for (i=0;i<n;i++) {
            errmax = max(errmax, fabs((y1[i] - y2[i]))/yscal[i]);
        }
    //cout << "x:" << *x << "\t\tp.a.:" << y[0]*180.0/M_PI << "\t\tV:" << tanh(2*y[1]) << "\t\tbeta_b + delta:" << ( BetaB(*x) + delta(*x) ) * 180.0 / M_PI + 90.0 << "\n";
    //cout << "x:" << *x << "\t\tratio:" << ( Lambda(*x) / Q(*x) ) / ( Lambda(*x) * cos(2.0 * (y[0] - BetaB(*x) - delta(*x))) * sinh(2 * y[1]) ) << "\n";
        errmax /= eps;
        if (errmax <= 1.0){
            *hdid=h;
            *hnext=(errmax > ERRCON ?
                    SAFETY*h*exp(PGROW*log(errmax)) : 4.0*h);

            // Making sure that stepsize is not too long

            if (*hnext > HMAX) *hnext = HMAX;
            break;
        }
        h=SAFETY*h*exp(PSHRNK*log(errmax));
    }
    for (i=0;i<n;i++)
        y[i] += ytemp[i]*FCOR;
}

#undef PGROW
#undef PSHRNK
#undef FCOR
#undef SAFETY
#undef ERRCON

void odeint(
        double ystart[],
        int nvar,
        double x1,
        double x2,
        double eps,
        double h1,
        double hmin,
        int nok,
        int nbad,
        void (*derivs)(double, double *, double *)
        ) {
    int nstp,i;
    double xsav, x, hnext, hdid, h;
    double BIG_EPS = 0.0;  
    int kmax = 2, kount = 0;
    double *xp, **yp, dxsav;
    xp = Nvector(0, kmax);
    yp = Nmatrix(0, nvar, 0, kmax);

    double * yscal;
    double * y;
    double * dydx;
    yscal = new double [nvar];
    y = new double [nvar];
    dydx = new double [nvar];

    string path = Globals::out_path + "/" + Globals::RUN_ID;
    ofstream plot(path + "_theta.dat");

    x=x1;
    h=(x2 > x1) ? fabs(h1) : -fabs(h1);
    nok = nbad = kount = 0;
    for (i=0;i<nvar;i++) y[i]=ystart[i];
    if (kmax > 0) xsav=x-dxsav*2.0;
    for (nstp=1;nstp<=MAXSTP;nstp++) {
        (*derivs)(x, y, dydx);
        for (i=0;i<nvar;i++)
            yscal[i]=fabs(y[i])+TINY + BIG_EPS; //+fabs(dydx[i]*h)
        if (kmax > 0) {
            if (fabs(x-xsav) > fabs(dxsav)) {
                if (kount < kmax-1) {
                    xp[++kount]=x;
                    for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
                    xsav=x;
                }
            }
        }
        if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;

        (rkqc)(y, dydx, nvar, &x, h, eps, yscal, &hdid, &hnext, derivs);

        plot << x <<  ", " << y[0] << ", " << y[1] << " " << constants::R_star * Globals::omega / (2.0 * constants::c) * Lambda(x) << endl;

        if (hdid == h) nok++; else nbad++;
        if ((x-x2)*(x2-x1) >= 0.0) {
            for (i=0; i<nvar; i++) {
                ystart[i]=y[i];
            }
            if (kmax) {
                xp[++kount]=x;
                for (i=1; i<=nvar; i++) yp[i][kount]=y[i];
            }
            plot.close();
            return;
        }
        if (fabs(hnext) <= hmin) {
            cout << "Step size too small in ODEINT\n";
            break;
        }
        h=hnext;
    }
    cout << "Too many steps in routine ODEINT\n";
    plot.close();
    delete [] yscal;
    delete [] y;
    delete [] dydx;
    }

#undef MAXSTP
#undef TINY
*/