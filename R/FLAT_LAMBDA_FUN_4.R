#' Log Posterior of Poisson Gap Parameter lambda_{kj} — Internal
#'
#' Computes the log posterior probability of the Poisson gap parameter \eqn{\lambda_{kj}}
#' between consecutive motif positions \eqn{j} and \eqn{j+1}, using a gamma prior
#' and a truncated Poisson likelihood.
#'
#' @param R A character vector of sequences.
#' @param W_samp A numeric vector of motif assignments.
#' @param lambda_kj The current value of the Poisson rate parameter \eqn{\lambda_{kj}}.
#' @param A_samp A list of motif position matrices.
#' @param pri_beta Shape parameter of the gamma prior.
#' @param pri_nu Rate parameter of the gamma prior.
#' @param rep_mt Integer. Index of the motif being evaluated.
#' @param rep_j Integer. Index of the current motif position.
#'
#' @return A numeric value: the log posterior probability of \code{lambda_kj}.
#' @keywords internal
log_prob_lambda_kj = function(R, W_samp, lambda_kj, A_samp,pri_beta, pri_nu, rep_mt, rep_j, motif_len_w){
  res = log(dgamma(lambda_kj, shape = pri_beta, rate = pri_nu))
  for (rep_i in 1:total_n) {
    wi = W_samp[rep_i]
    if(wi == rep_mt){Ri = R[rep_i]
    len_seq = str_length(Ri)
    aij = A_samp[[rep_mt]][rep_i,rep_j]
    aij1 = A_samp[[rep_mt]][rep_i,rep_j+1]


    max_y = len_seq - motif_len_w + (rep_j+1) - aij - 1

    y_vals = 0:max_y
    poisson_probs = dpois(y_vals, lambda_kj)  # 泊松分布概率
    normalization_constant = sum(poisson_probs)  # 分母

    res = res + log(poisson_probs[aij1 - aij]) - log(normalization_constant)}
  }
  return(res)
}

#' Metropolis-Hastings Update for Lambda Vector — Internal
#'
#' Performs a Metropolis-Hastings update step for the Poisson gap parameters
#' \eqn{\lambda_{kj}} in a motif. A log-normal proposal is used with acceptance
#' probability based on the log posterior.
#'
#' @param R A character vector of sequences.
#' @param W_samp A numeric vector of motif assignments.
#' @param lambda_samp A list of current lambda vectors for each motif.
#' @param A_samp A list of motif position matrices.
#' @param pri_beta Shape parameter of the gamma prior.
#' @param pri_nu Rate parameter of the gamma prior.
#' @param rep_mt Integer. Index of the motif to update.
#'
#' @return A numeric vector of updated lambda values for the selected motif.
#' @keywords internal
lambda_samp_fun = function(R,W_samp, lambda_samp, A_samp,pri_beta, pri_nu, rep_mt, motif_len_w){
  lambda_samp_unobs = rep(NA,motif_len_w-1)

  for (rep_j in 1:(motif_len_w-1)) {
    lambda_kj = lambda_samp[[rep_mt]][rep_j]
    log_lambda_kj = log(lambda_kj)

    log_lambda_kj_star = log_lambda_kj + rnorm(1,0,.1)
    lambda_kj_star = exp(log_lambda_kj_star)

    ## mh accecpt algorithm
    log_pi = (log_prob_lambda_kj(R, W_samp, lambda_kj_star, A_samp,pri_beta, pri_nu, rep_mt, rep_j, motif_len_w)-
                log_prob_lambda_kj(R,W_samp, lambda_kj, A_samp,pri_beta, pri_nu, rep_mt, rep_j, motif_len_w))
    alpha = min(1, exp(log_pi));s = runif(1)
    if(s>alpha){lambda_samp_unobs[rep_j] = lambda_kj}else{
      lambda_samp_unobs[rep_j] = lambda_kj_star
    }
  }

  return(lambda_samp_unobs)
}
