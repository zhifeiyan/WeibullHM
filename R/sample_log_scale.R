#' Sample Weibull log scale parameter
#' 
#' Produce a Gibbs sample of Weibull log scale parameter of a subject in a state
#' 
#' @param rt a vector of response time sequence of a subject
#' @param st a vector of current update of state process of a subject
#' @param y current value of Weibull log scale parameter
#' @param wb_shape current value of Weibull shape parameter
#' @param mu current value of mean of Weibull log scale parameter
#' @param sigma a vector of standard deviations of Weibull log scale parameters
#' @param proposal_std a vector of random walk proposal standard deviations 
#' @param s a specific state
#' 
#' @return A vector whose first element is a Gibbs sample and second element 
#'         is a binary indicator value of acceptance or rejection of the proposed 
#'         sample in random walk Metropolis–Hastings algorithm
#' @export
sample_log_scale <- function(rt, st, y, wb_shape, mu, sigma, proposal_std, s){
  sigma <- sigma[s]
  proposal_std <- proposal_std[s]

  y_old <- y
  y_new <- rnorm(1, y_old, proposal_std)

  log_ratio <- (sum(dweibull(rt[st == s - 1], wb_shape, 
                             exp(y_new), log = TRUE)) + 
                dnorm(y_new, mu, sigma, log = TRUE)) - 
               (sum(dweibull(rt[st == s - 1], wb_shape, 
                             exp(y_old), log = TRUE)) + 
                dnorm(y_old, mu, sigma, log = TRUE))

  if(log(runif(1, 0, 1)) < log_ratio) {
    return(c(y_new, 1))
  } else {
    return(c(y_old, 0))
  }
}
