#include <RcppArmadillo.h>

using namespace arma;

//' Compute stationary distributions of environment process and state process
//' 
//' @param p an array of transition matrices of state process
//' @param q a matrix of transition matrix of environment process
//' @param nstate total number of states
//' 
//' @return A list with the following components:
//'   \item{e_stat}{stationary distribution of environment process}
//'   \item{s_stat}{stationary distribution of state process}
//' @export
// [[Rcpp::export]]
Rcpp::List stat_dist(arma::cube p, arma::mat q, int nstate) {
  int d = nstate * nstate, e_prev, e_cur, s_prev, s_cur, i, j;
  mat trans(d, d), joint_stat(nstate, nstate);
  vec b(d, fill::zeros);
  b(d - 1) = 1;

  for (i = 0; i < d; i++) {
    e_prev = i / nstate;
    s_prev = i % nstate;
    for (j = 0; j < d; j++) {
      e_cur = j  / nstate;
      s_cur = j % nstate;
      trans(i, j) = p.slice(e_cur)(s_prev, s_cur) * q(e_prev, e_cur);
    }
  }
  trans.diag() -= 1;
  trans.col(d - 1).ones();
  joint_stat = reshape(solve(trans.t(), b), nstate, nstate);

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("e_stat") = sum(joint_stat), 
    Rcpp::Named("s_stat") = sum(joint_stat, 1));
  return out;
}