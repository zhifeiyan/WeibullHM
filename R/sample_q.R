#
# sample_q.R
#
# Created by Zhifei Yan
# Last update 2016-1-15
#

# This function calculates the stationary distribution of an nstate-Markov
# chain with a transition matrix trans
stat_dist_onechain <- function(trans, nstate) {
  diag(trans) <- diag(trans) - 1
  trans[, nstate] <- 1
  b <- c(rep(0, nstate - 1), 1)
  return(solve(t(trans), b))
}

#' Sample transition matrix of environment process
#' 
#' Produce a Gibbs sample of transition matrix of environment process of a subject
#' 
#' @param q a matrix of current update of transition matrix of environment process
#' @param et a vector of current update of environment process of a subject
#' @param prior a matrix of concentration parameters of Dirichlet priors,
#'        each row corresponds to a row of transition matrix
#' @param ntrial total number of trials
#' @param nenv total number of environments
#' 
#' @return A list with the following components:
#'  \item{q}{a Gibbs sample of transition matrix of environment process}
#'  \item{accept}{a binary indicator vector of acceptance or rejection 
#'                of proposed sample of each row of transition matrix 
#'                in the independent Metropolis–Hastings algorithm}
#' @export
sample_q <- function(q, et, prior, ntrial, nenv) {
  e_stat_old <- stat_dist_onechain(q, nenv)
  q_temp <- q
  accept_vec <- rep(0, nenv)
  et_droplast <- et[-ntrial]

  # Update each of the rows of q matrix sequentially
  for (k in 0:(nenv - 1)) {
    alpha <- prior[k + 1, ]
    id <- which(et_droplast == k)
    if(length(id) == 0) {
      n_vec <- rep(0, nenv)
    } else {
      et_subset <- et[id + 1]
      n_vec <- rep(NA, nenv)
      for (i in 0:(nenv - 1)) {
        n_vec[i + 1] <- sum(et_subset == i)
      }
    }
    q_temp[k + 1, ] <- as.vector(rdirichlet(1, n_vec + alpha))
    e_stat_new <- stat_dist_onechain(q_temp, nenv)
    log_ratio <- log(e_stat_new[et[1] + 1]) - log(e_stat_old[et[1] + 1])
    if(log(runif(1, 0, 1)) < log_ratio) {
      accept_vec[k + 1] <- 1
      q <- q_temp
      e_stat_old <- e_stat_new
    } else {
      q_temp <- q
    }
  }
  return(list(q = q, accept = accept_vec))
}