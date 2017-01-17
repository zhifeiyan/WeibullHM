//
// sample_et_st.cpp
//
// Created by Zhifei Yan
// Last update 2017-1-14
//

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
//' @param nenv   total number of environments
//' 
//' @return A list with the following components:
//'   \item{et}{a vector of Gibbs sample of environment process}
//'   \item{st}{a vector of Gibbs sample of state process}
//' @export
// [[Rcpp::export]]
Rcpp::List sample_et_st(arma::vec rt, double wb_shape, arma::vec wb_scale, 
                        arma::cube p, arma::mat q, arma::vec e_stat, 
                        arma::vec s_stat, int ntrial, int nstate, 
                        int nenv) {
  int k, i, j, t;
  mat p_use(nstate, nstate),
  p_e(nenv, nstate),
  p_s(nenv, nstate),
  cond_temp(nenv, nstate),
  temp_prop(nenv, nstate),
  temp(nenv, nstate);

  vec p_e1(nstate, fill::zeros),
  cond_temp_0(nstate),
  temp_prop_0(nstate),
  temp_0(nstate);
  rowvec lik(nstate);

  // Placeholder for prob dists used in conditional simulation pieces
  mat cond_1(nstate, nenv * nstate), cond_last(nenv * nstate, nenv);
  cube cond_middle(nenv * nstate, nenv * nstate, ntrial - 2);
  // Placeholder for current updated filter after the first updating step
  mat filter_update(nenv, nstate);

  // Start FFBS dynamic programming procedure
  // Forward filtering
  // Compute filter_0
  temp_prop_0 = my_dweibull(rt(0), wb_shape, wb_scale, 0) % s_stat;
  temp_0 = temp_prop_0 / accu(temp_prop_0);

  // From filter_0 to filter_1, also compute cond_1
  k = 0;
  for (i = 0; i < nenv; i++) {
    p_use = p.slice(i);
    p_e1 += e_stat(i);
    for (j = 0; j < nstate; j ++) {
      cond_temp_0 = p_use.col(j) % p_e1 % temp_0;
      filter_update(i, j) = sum(cond_temp_0);
      cond_1.col(k) = cond_temp_0 / filter_update(i, j);
      k++;
    }
    p_e1.zeros();
  }

  lik = my_dweibull(rt(1), wb_shape, wb_scale, 0).t();
  temp_prop = repmat(lik, nenv, 1) % filter_update;
  temp = temp_prop / accu(temp_prop);

  // From filter_1 to filter_(T-1), also compute T - 2 conditional probs in cond_middle
  for(t = 1; t < ntrial - 1; t++) {
    k = 0;
    for(i = 0; i < nenv; i++) {
      p_use = p.slice(i);
      p_e = repmat(q.col(i), 1, nstate);
      for(j = 0; j < nstate; j++) {
        p_s = repmat(p_use.col(j).t(), nenv, 1);
        cond_temp = p_s % p_e % temp;
        filter_update(i, j) = accu(cond_temp);
        cond_middle.slice(t - 1).col(k) = vectorise((cond_temp / filter_update(i, j)).t());
        k++;
      }
    }
    lik = my_dweibull(rt(t + 1), wb_shape, wb_scale, 0).t();
    temp_prop = repmat(lik, nenv, 1) % filter_update;
    temp = temp_prop / accu(temp_prop);
  }

  // From filter_(T-1) to filter_T, also compute cond_last
  // Need to store filter_T as the starting point of backward sampling
  vec filter_T(nenv);
  for(i = 0; i < nenv; i++) {
    cond_temp = repmat(q.col(i), 1, nstate) % temp;
    filter_T(i) = accu(cond_temp);
    cond_last.col(i) = vectorise((cond_temp / filter_T(i)).t());
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
