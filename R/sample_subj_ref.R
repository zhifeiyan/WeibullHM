#' Sample subject random effects
#' 
#' Produce a Gibbs sample of subject random effect of a subject in a state
#' 
#' @param y a matrix of current update of Weibull log scale parameters
#' @param mu a vector of current update of intercept effects in all states
#' @param sigma a vector of standard deviations of Weibull log scale parameters
#' @param ref_std a vector of current update of standard deviations of subject random effects
#' @param i a specific subject
#' @param s a specific state
#' @return A Gibbs sample of subject random effect of a subject in a state
#' @export
sample_subj_ref <- function(y, mu, sigma, ref_std, i, s) {
  y <- y[i, s]
  mu <- mu[s]
  ref_var <- ref_std[s] ^ 2
  sigma_sq <- sigma[s] ^ 2

  rnorm(1, (y - mu) * ref_var / (sigma_sq + ref_var), 
        sqrt(ref_var * sigma_sq / (sigma_sq + ref_var)))
}