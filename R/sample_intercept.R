#' Sample intercept effects
#' 
#' Produce a Gibbs sample of intercept effects of all states
#' 
#' @param y a matrix of current update of Weibull log scale parameters
#' @param alpha a matrix of current update of subject random effects
#' @param sigma_lambda a vector of standard deviations of Weibull log scale parameters
#' @param m_mu a vector of means of multivariate normal prior of intercept effects
#' @param sigma_mu_inv inverse covariance matrix of multivariate normal prior 
#'        of intercept effects
#' @param nsubj total number of subjects
#' 
#' @return A Gibbs sample of intercept effects of all states
#' @export
sample_intercept <- function(y, alpha, sigma_lambda, m_mu, 
                             sigma_mu_inv, nsubj) {
  xty <- apply(y - alpha, 2, sum) / (sigma_lambda ^ 2)
  xtx <- nsubj * diag(1 / sigma_lambda ^ 2)

  cov_post <- solve(sigma_mu_inv + xtx)
  m_post <- cov_post %*% (sigma_mu_inv %*% m_mu + xty)

  draw <- mvrnorm(1, m_post, cov_post)
  while (! all(draw == cummax(draw))) {
    draw <- mvrnorm(1, m_post, cov_post)
  }

  draw
}