#' Gibbs Update for Latent Motif Positions (A matrix) — Internal
#'
#' Performs a Gibbs sampling update of the latent motif positions (A matrix)
#' for a specific motif. This function uses a Poisson gap model and position-specific
#' likelihoods based on the motif parameters.
#'
#' @param R A character vector of sequences.
#' @param W_samp A numeric vector indicating motif label assignments for each sequence.
#' @param lambda_samp A list of Poisson gap parameters for each motif.
#' @param A_samp A list of A matrices (motif positions) for each motif.
#' @param dict A character vector representing the sequence alphabet (e.g., \code{c("A", "C", "G", "T")}).
#' @param theta_0_samp A vector of background letter probabilities.
#' @param theta_samp A list of position-specific motif probability matrices (PWMs).
#' @param rep_mt Integer. Index of the motif being updated.
#' @param col_range A vector of motif position indices to update (e.g., 1:motif_len_w).
#' @param row_range A vector of sequence indices to update.
#'
#' @return A matrix of updated motif positions for the selected motif.
#' @keywords internal
A_samp_fun = function(R, W_samp, lambda_samp, A_samp,dict,theta_0_samp,theta_samp, rep_mt, col_range, row_range, motif_len_w){
  A_samp_unobs = A_samp[[rep_mt]]
  for (rep_i in row_range) {
    Ri = R[rep_i]
    len_seq = str_length(Ri)
    tmp_index = 1:str_length(Ri)
    Ri = str_sub(Ri,tmp_index,tmp_index)
    wi = W_samp[rep_i]

    if(wi == rep_mt){
      for (rep_j in col_range) {
        if(rep_j == 1){
          Range = 1:(A_samp_unobs[rep_i,rep_j+1]-1)
        }else if((rep_j > 1)&(rep_j < motif_len_w)){
          Range = (A_samp_unobs[rep_i,rep_j-1]+1):(A_samp_unobs[rep_i,rep_j+1]-1)
        }else{
          Range = (A_samp_unobs[rep_i,rep_j-1]+1):len_seq
        }

        if(length(Range)>1){
          res = rep(1,length(Range))
          for (rep_Ai in Range) {
            A_potential = rep(NA,motif_len_w)
            A_potential[rep_j] = rep_Ai
            A_potential[-rep_j] = A_samp_unobs[rep_i,-rep_j]

            Res = logllk_fun_kj(Ri,rep_mt,rep_j,A_potential,theta_0_samp,theta_samp,dict)


            if((rep_j<motif_len_w)){
              max_y = len_seq - motif_len_w + (rep_j+1) - (rep_Ai) - 1
              y_vals = 0:max_y
              poisson_probs = dpois(y_vals, lambda_samp[[rep_mt]][rep_j])
              normalization_constant = sum(poisson_probs)
              shifted_probs = (poisson_probs / normalization_constant)
              Res = (Res + log(poisson_probs[A_samp_unobs[rep_i,(rep_j+1)]-rep_Ai])-
                       log(normalization_constant))
            }
            if((rep_j>1)){
              poisson_prob = dpois((rep_Ai -  A_samp_unobs[rep_i,rep_j-1] - 1), lambda_samp[[rep_mt]][rep_j-1])
              Res = Res + log(poisson_prob)
            }

            res[which(Range==rep_Ai)] = exp(Res)
          }
          post_prob = res/sum(res)

          a = sample(Range,size = 1, prob = post_prob)
          A_samp_unobs[rep_i,rep_j] = a
        }else{
          A_samp_unobs[rep_i,rep_j] = Range
        }
      }
    }
  }
  return(A_samp_unobs)
}
