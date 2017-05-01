#
# sample_var_logscale.R
#
# Created by Zhifei Yan
# Last update 2017-4-22
#

#' Sample the variance of Weibull log scale parameter in a state
#' 
#' Produce a Gibbs sample of the variance of Weibull log scale parameter in a state
#' 
#' @param log_scale a vector of log Weibull scale parameters across all 
#'                  subjects in a state
#' @param mean_log_scale a vector of means of log Weibull scale parameters
#'                       across all subjects in a state
#' @param nsubj total number of subjects
#' @param prior_invgamma_shape shape parameter of the inverse-gamma prior
#' @param prior_invgamma_scale scale parameter of the inverse-gamma prior
#' 
#' @return A sample from the full conditional of the variance of Weibull log 
#'         scale parameter in a state
#' @export
sample_var_logscale <- function(log_scale, mean_log_scale, nsubj, 
                                prior_invgamma_shape, prior_invgamma_scale) {
  rinvgamma(1, shape = prior_invgamma_shape + nsubj / 2, 
            scale = 0.5 * sum((log_scale - mean_log_scale)^2) + 
            prior_invgamma_scale)
}