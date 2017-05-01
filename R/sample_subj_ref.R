#
# sample_subj_ref.R
#
# Created by Zhifei Yan
# Last update 2017-4-22
#

#' Sample subject random effects
#' 
#' Produce a Gibbs sample of subject random effect of a subject in a state
#' 
#' @param y a matrix of current update of Weibull log scale parameters
#' @param mu a vector of current update of intercept effects in all states
#' @param var_logscale a vector of current update of variances of Weibull log scale parameters in all states
#' @param ref_std a vector of current update of standard deviations of subject random effects
#' @param i a specific subject
#' @param s a specific state
#' @return A Gibbs sample of subject random effect of a subject in a state
#' @export
sample_subj_ref <- function(y, mu, var_logscale, ref_std, i, s) {
  y <- y[i, s]
  mu <- mu[s]
  ref_var <- ref_std[s] ^ 2
  var_logscale <- var_logscale[s]

  rnorm(1, (y - mu) * ref_var / (var_logscale + ref_var), 
        sqrt(ref_var * var_logscale / (var_logscale + ref_var)))
}