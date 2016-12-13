#
# sample_ref_std.R
#
# Created by Zhifei Yan
# Last update 2016-12-13
#

#' Sample the standard deviation of subject random effects
#' 
#' Produce a Gibbs sample of the standard deviation of subject random effects in a state
#' 
#' @param alpha a matrix of current update of subject random effects
#' @param ref_std current value of the standard deviations of subject random effects
#' @param s a specific state
#' @param nsubj total number of subjects
#' @param prior_hc_scale scale parameter of half-Cauchy prior
#' @param prop_std random walk proposal standard deviation
#' 
#' @return A vector whose first element is a Gibbs sample and second element 
#'         is a binary indicator value of acceptance or rejection of the proposed 
#'         sample in random walk Metropolis–Hastings algorithm
#' @export
sample_ref_std <- function(alpha, ref_std, s, nsubj, prior_hc_scale, 
                           prop_std) {
  ref_std_old <- ref_std[s]
  ref_std_new <- rnorm(1, ref_std_old, prop_std[s])
  if (ref_std_new <= 0) return(c(ref_std_old, 0))
  
  log_ratio <- - nsubj * log(ref_std_new) - 
               sum(alpha[, s] ^ 2) / 2 / ref_std_new ^ 2 - 
               log(ref_std_new ^ 2 + prior_hc_scale ^ 2) + 
               nsubj * log(ref_std_old) + 
               sum(alpha[, s] ^ 2) / 2/ ref_std_old ^ 2 + 
               log(ref_std_old ^ 2 + prior_hc_scale ^ 2)
  if (log(runif(1, 0, 1)) < log_ratio) {
    return(c(ref_std_new, 1))
  } else {
    return(c(ref_std_old, 0))
  }
}