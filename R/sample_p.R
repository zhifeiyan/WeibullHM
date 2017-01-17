#
# sample_p.R
#
# Created by Zhifei Yan
# Last update 2016-1-15
#

#' Sample transition matrices of state process
#' 
#' Produce a Gibbs sample of transition matrices of state process under all environments of a subject
#' 
#' @param p_all a matrix of current update of transition matrices of 
#'              state process, transition matrices under different 
#'              environments are stacked vertically
#' @param q     a matrix of current update of transition matrix of 
#'              environment process
#' @param prior_mat a matrix of concentration parameters of Dirichlet priors,
#'                  each row corresponds to a row of transition matrix
#' @param et a vector of current update of environment process of a subject
#' @param st a vector of current update of state process of a subject
#' @param ntrial total number of trials
#' @param nstate total number of states
#' @param nenv   total number of environments
#' 
#' @return A list with the following components:
#'  \item{p}{a Gibbs sample of transition matrices of state process}
#'  \item{accept}{a binary indicator vector of acceptance or rejection 
#'                of proposed sample of each row of transition matrices 
#'                in the independent Metropolis–Hastings algorithm}
#' @export
sample_p <- function(p_all, q, prior_mat, et, st, ntrial, nstate, nenv){
  # p matrices transformation relation between two different data structure
  p_array <- aperm(array(t(p_all), c(nstate, nstate, nenv)), c(2, 1, 3))
  s_stat_old <- stat_dist(p_array, q, nstate, nenv)$s_stat
  p_all_temp <- p_all
  accept_vec <- rep(0, nenv * nstate)
  et_st_id <- (nstate * et + st)[-ntrial]
  
  # Update each of the rows in all p matrices sequentially
  for (k in 0:(nenv * nstate - 1)){
    alpha <- prior_mat[k + 1, ]
    id <- which(et_st_id == k)
    if(length(id) == 0){
      n_vec <- rep(0, nstate)
    } else {
      st_subset <- st[id + 1]
      n_vec <- rep(NA, nstate)
      for (i in 0:(nstate - 1)) {
        n_vec[i + 1] <- sum(st_subset == i)
      }
    }
    p_new <- as.vector(rdirichlet(1, n_vec + alpha))
    p_all_temp[k + 1, ] <- p_new
    p_array_temp <- aperm(array(t(p_all_temp), c(nstate, nstate, nenv)), 
                          c(2, 1, 3))
    s_stat_new <- stat_dist(p_array_temp, q, nstate, nenv)$s_stat

    log_ratio <- log(s_stat_new[st[1] + 1]) - log(s_stat_old[st[1] + 1])
    if(log(runif(1,0,1)) < log_ratio) {
      accept_vec[k + 1] <- 1
      p_all <- p_all_temp
      s_stat_old <- s_stat_new
    } else {
      p_all_temp <- p_all
    }
  }
  return(list(p = p_all, accept = accept_vec))
}