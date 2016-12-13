#include <RcppArmadillo.h>

using namespace arma;

// This function generates a sample from some given discrete probability distribution
int sample_singleton(vec p) {
  int elem;
  double x, u;

  u = R::runif(0, 1);

  elem = 0;
  x = p(0);
  while(x < u) {
    elem++;
    x += p(elem);
  }
  // Return 0, 1, ... , p.n_elem - 1
  return(elem);
}
