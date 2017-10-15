#
# sample_std_logscale.R
#
# Created by Zhifei Yan
# Last update 2017-10-13
#

#' Sample standard deviation of Weibull log scale parameters in a given state
#' 
#' Produce a Gibbs sample of standard deviation of Weibull log scale parameters in a given state
#' 
#' @param log_scale a vector of Weibull log scale parameters for all subjects 
#'                  in the given state
#' @param mean_log_scale a vector of means of Weibull log scale parameters
#'                       for all subjects in the given state
#' @param std current value of the standard deviation
#' @param nsubj total number of subjects
#' @param prior_hc_scale scale parameter of half-Cauchy prior
#' @param prop_std a vector of random walk proposal standard deviations for all states
#' @param s a specific state
#' 
#' @return A sample from the full conditional of standard deviation of Weibull log 
#'         scale parameter in a given state
#' @export
sample_std_logscale <- function(log_scale, mean_log_scale, std, nsubj, 
                                prior_hc_scale, prop_std, s) {
  std_old <- std
  std_new <- rnorm(1, std_old, prop_std[s])

  if (std_new <= 0) return(c(std_old, 0))

  log_ratio <- sum(dnorm(log_scale, mean_log_scale, std_new, log = TRUE)) -
               log(prior_hc_scale ^ 2 + std_new ^ 2) -
               sum(dnorm(log_scale, mean_log_scale, std_old, log = TRUE)) + 
               log(prior_hc_scale ^ 2 + std_old ^ 2)

  if (log(runif(1, 0, 1)) < log_ratio) {
    return(c(std_new, 1))
  } else {
    return(c(std_old, 0))
  }
}

# sample_var_logscale <- function(log_scale, mean_log_scale, nsubj, 
#                                 prior_invgamma_shape, prior_invgamma_scale) {
#   rinvgamma(1, shape = prior_invgamma_shape + nsubj / 2, 
#             scale = 0.5 * sum((log_scale - mean_log_scale)^2) + 
#             prior_invgamma_scale)
# }
