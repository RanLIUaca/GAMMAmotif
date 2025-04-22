#' Sample Motif Position Weight Matrix Theta — Internal
#'
#' Samples the position weight matrix (PWM) \code{theta} for a specific motif,
#' using a Dirichlet posterior based on observed positions in assigned sequences.
#'
#' @param R A character vector of sequences.
#' @param W_samp A numeric vector of motif assignments.
#' @param A_samp A list of motif position matrices.
#' @param dict A character vector of the sequence alphabet.
#' @param pri_alw A numeric vector of prior Dirichlet concentration parameters (motif-wise).
#' @param rep_mt Integer. Index of the motif to sample theta for.
#'
#' @return A numeric matrix of sampled motif PWM for the selected motif.
#' @keywords internal
theta_samp_fun = function(R,W_samp,A_samp,dict,pri_alw,rep_mt, motif_len_w){
  A_mt = A_samp[[rep_mt]]
  total_n=length(R)

  # sample parameter
  theta_samp = matrix(nrow = length(dict), ncol = motif_len_w)
  for (rep_j in 1:motif_len_w) {
    #rep_j=1
    # post_parameter
    H_a=rep(0,length(dict))
    names(H_a)=dict
    for (rep_i in 1:total_n) {
      #rep_i=1
      Ri = R[rep_i]
      tmp_index = 1:str_length(Ri)
      Ri = str_sub(Ri,tmp_index,tmp_index)
      wi = W_samp[rep_i]

      if(wi==rep_mt){
        aij = A_mt[rep_i,rep_j]
        Ri_aij = Ri[aij]
        tmp_index = 1:str_length(Ri_aij)
        Ri_aij = str_sub(Ri_aij,tmp_index,tmp_index)
        Ri_aij = factor(Ri_aij, levels = dict)
        H_a=H_a+table(Ri_aij)
      }
    }

    post_alpha_w = pri_alw[rep_mt]+H_a

    temp = rdirichlet(1,post_alpha_w)
    theta_samp[,rep_j] = temp
  }

  return(theta_samp)
}
