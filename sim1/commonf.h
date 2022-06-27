#ifndef COMMONF_H
#define COMMONF_H

#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <Rmath.h>
#include <iostream>

using namespace Rcpp;

double hpwc_double(NumericVector x, NumericVector cuts, NumericVector levels, int logInd);
double Hpwc_double(NumericVector x, NumericVector cuts, NumericVector levels, int logInd);
double ppwc_double(double q, NumericVector cuts, NumericVector levels, int lower, int logInd);
double pi11_jkjpkp(double marg_jk, double marg_jpkp, double gamma);
#endif




