#' Sample Weibull log shape parameter
#' 
#' Produce a Gibbs sample of Weibull log shape parameter of a subject
#' 
#' @param rt a vector of response time sequence of a subject
#' @param st a vector of current update of state process of a subject
#' @param y a vector of current values of Weibull log scale parameters in all states
#' @param eta current value of Weibull log shape parameter of a subject
#' @param prior_mu mean of normal prior distribution of Weibull log shape parameter
#' @param prior_sigma standard deviation of normal prior distribution of Weibull
#'        log shape parameter
#' @param proposal_std random walk proposal standard deviation
#' @param nstate total number of states
#' 
#' @return A vector whose first element is a Gibbs sample and second element 
#'         is a binary indicator value of acceptance or rejection of the proposed 
#'         sample in random walk Metropolis–Hastings algorithm
#' @export
sample_log_shape <- function(rt, st, y, eta, prior_mu, 
                             prior_sigma, proposal_std, nstate){
  eta_old <- eta
  eta_new <- rnorm(1, eta_old, proposal_std)

  log_ratio <- dnorm(eta_new, prior_mu, prior_sigma, log = TRUE) - 
               dnorm(eta_old, prior_mu, prior_sigma, log = TRUE)

  for (s in 1:nstate) {
    log_ratio <- log_ratio + 
                 sum(dweibull(rt[st == s - 1], exp(eta_new), 
                              exp(y[s]), log = TRUE)) - 
                 sum(dweibull(rt[st == s - 1], exp(eta_old), 
                              exp(y[s]), log = TRUE))
  }

  if(log(runif(1,0,1)) < log_ratio){
    return(c(eta_new, 1))
  } else {
    return(c(eta_old, 0))
  }

}