#ifndef __MYDWEIBULL_H
#define __MYDWEIBULL_H

#include <RcppArmadillo.h>

using namespace arma;
vec my_dweibull(double val, double shape, vec vec_scale, int log);

#endif


