#include <RcppArmadillo.h>
#include <cmath>

using namespace arma;

// This function evaluates the density of Weibull distribution for a vector of 
// scale paramaters
vec my_dweibull(double val, double shape, vec vec_scale, int log) {
  int n = vec_scale.n_elem;
  vec res(n);
  for (int i = 0; i < n; i++) {
    res[i] = R::dweibull(val, shape, vec_scale[i], log);
  }
  return res;
}