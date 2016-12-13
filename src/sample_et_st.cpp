// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

#include "sample_singleton.h"
#include "my_dweibull.h"

using namespace arma;

//' Sample environment process and state process of a subject
//'
//' Produce a Gibbs sample of environment process and state process of a subject
//' 
//' @param rt a vector of response time sequence of a subject
//' @param wb_shape current update of Weibull shape parameter
//' @param wb_scale a vector of current update of Weibull scale parameters in all states
//' @param p an array of current update of transition matrices of state process
//' @param q a matrix of current update of transition matrix of environment process
//' @param e_stat a vector of current update of stationary distribution of environment process
//' @param s_stat a vector of current update of stationary distribution of state process
//' @param ntrial total number of trials
//' @param nstate total number of states
//' 
//' @return A list with the following components:
//'   \item{et}{a vector of Gibbs sample of environment process}
//'   \item{st}{a vector of Gibbs sample of state process}
//' @export
// [[Rcpp::export]]
Rcpp::List sample_et_st(arma::vec rt, double wb_shape, arma::vec wb_scale, 
                        arma::cube p, arma::mat q, arma::vec e_stat, 
                        arma::vec s_stat, int ntrial, int nstate) {
  int k, i, j, t;
  mat p_use(nstate, nstate), p_e(nstate, nstate), p_s(nstate, nstate), 
  cond_temp(nstate, nstate), temp_prop(nstate, nstate), temp(nstate, nstate);
  vec p_e1(nstate, fill::zeros), cond_temp_0(nstate), temp_prop_0(nstate), 
  temp_0(nstate);
  rowvec lik(nstate);

  // Placeholder for prob dists used in conditional simulation pieces
  mat cond_1(nstate, nstate * nstate), cond_last(nstate * nstate, nstate);
  cube cond_middle(nstate * nstate, nstate * nstate, ntrial - 2);
  // Placeholder for current updated filter
  mat filter_update(nstate, nstate);

  // Start FFBS dynamic programming procedure
  // Forward filtering
  // Compute filter_0
  temp_prop_0 = my_dweibull(rt(0), wb_shape, wb_scale, 0) % s_stat;
  temp_0 = temp_prop_0 / accu(temp_prop_0);

  // From filter_0 to filter_1, also compute cond_1
  k = 0;
  for (i = 0; i < nstate; i++) {
    p_use = p.slice(i);
    p_e1 += e_stat(i);
    for (j = 0; j < nstate; j ++) {
      cond_temp_0 = p_use.col(j) % p_e1 % temp_0;
      filter_update(i, j) = sum(cond_temp_0);
      cond_1.col(k) = cond_temp_0 / filter_update(i, j);
      k++;
    }
  }

  lik = my_dweibull(rt(1), wb_shape, wb_scale, 0).t();
  temp_prop = repmat(lik, nstate, 1) % filter_update;
  temp = temp_prop / accu(temp_prop);

  // From filter_1 to filter_(T-1), also compute conditional probs in cond_middle
  for(t = 1; t < ntrial - 1; t++) {
    k = 0;
    for(i = 0; i < nstate; i++) {
      p_use = p.slice(i);
      p_e = repmat(q.col(i), 1, nstate);
      for(j = 0; j < nstate; j++) {
        p_s = repmat(p_use.col(j), 1, nstate).t();
        cond_temp = p_s % p_e % temp;
        filter_update(i, j) = accu(cond_temp);
        cond_middle.slice(t - 1).col(k) = vectorise((cond_temp / filter_update(i, j)).t());
        k++;
      }
    }
    lik = my_dweibull(rt(t + 1), wb_shape, wb_scale, 0).t();
    temp_prop = repmat(lik, nstate, 1) % filter_update;
    temp = temp_prop / accu(temp_prop);
  }

  // From filter_(T-1) to filter_T, also compute cond_last
  // Need to store filter_T as the starting point of backward sampling
  vec filter_T(nstate);
  k = 0;
  for(i = 0; i < nstate; i++) {
    cond_temp = repmat(q.col(i), 1, nstate) % temp;
    filter_T(i) = accu(cond_temp);
    cond_last.col(k) = vectorise((cond_temp / filter_T(i)).t());
    k++;
  }

  // Backward sampling
  int etst_combo;
  Rcpp::IntegerVector et(ntrial), st(ntrial);
  et(ntrial - 1) = sample_singleton(filter_T);
  etst_combo = sample_singleton(cond_last.col(et(ntrial - 1)));
  et(ntrial - 2) = etst_combo / nstate;
  st(ntrial - 1) = etst_combo % nstate;
  for (i = ntrial - 3; i > -1; i--) {
    etst_combo = sample_singleton(cond_middle.slice(i).col(etst_combo));
    et(i) = etst_combo / nstate;
    st(i + 1) = etst_combo % nstate;
  }
  st(0) = sample_singleton(cond_1.col(etst_combo));

  Rcpp::List out = Rcpp::List::create(Rcpp::Named("et") = et, 
                                      Rcpp::Named("st") = st);
  return out;
}
