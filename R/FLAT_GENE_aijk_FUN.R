#' Sample a_{ijk} from Shifted Truncated Poisson — Internal
#'
#' Samples the motif position \eqn{a_{ijk}} based on a Poisson-distributed gap 
#' from the previous position \eqn{a_{ij(k-1)}}. This is used in motif models 
#' where parts of the motif are separated by variable-length gaps.
#'
#' @param a_prev Integer. The previous motif position \eqn{a_{ij(k-1)}}.
#' @param lambda Numeric. The Poisson rate parameter for the gap.
#' @param L_i Integer. Length of the current sequence.
#' @param J_k Integer. Total number of motif positions.
#' @param j Integer. Index of the current motif position being sampled.
#'
#' @return An integer representing the sampled motif position \eqn{a_{ijk}}.
#' @keywords internal
generate_aijk <- function(a_prev, lambda, L_i, J_k, j) {

  max_x <- L_i - J_k + j - a_prev - 1
  if (max_x < 0) stop("Invalid range for a_{ijk}: max_x is negative")
  
  x_vals <- 0:max_x
  poisson_probs <- dpois(x_vals, lambda)  
  normalization_constant <- sum(poisson_probs)  
  
  possible_aijk <- a_prev + 1 + x_vals
  shifted_probs <- poisson_probs / normalization_constant 
  
  if(shifted_probs[1] == 1) {sampled_aijk = possible_aijk} else{
    sampled_aijk <- sample(possible_aijk, size = 1, prob = shifted_probs)
  }
  return(sampled_aijk)
}