#' Sample Latent Motif Assignments W — Internal
#'
#' Samples the unobserved motif assignment variable \code{W} for a given set of sequences,
#' based on posterior probabilities using current estimates of background and motif parameters.
#'
#' @param R A character vector of sequences.
#' @param UW_loc A numeric vector of sequence indices for which \code{W} needs to be sampled.
#' @param A_samp A list of motif position matrices.
#' @param dict A character vector of the alphabet (e.g., \code{c("A", "C", "G", "T")}).
#' @param theta_0_samp A numeric vector of background probabilities.
#' @param theta_samp A list of position weight matrices (PWMs) for motifs.
#' @param pri_wp A numeric vector of prior probabilities for each motif (including background).
#'
#' @return A numeric vector of sampled motif assignments for the unobserved entries.
#' @keywords internal
w_samp_fun = function(R,UW_loc,A_samp,dict,theta_0_samp,theta_samp,pri_wp){
  W_samp_unobs = rep(NA, length(UW_loc))
  for (rep_i in UW_loc) {
    # rep_i=4
    Rui = R[rep_i]
    tmp_index = 1:str_length(Rui)
    Rui = str_sub(Rui,tmp_index,tmp_index)
    
    llk_one_seq = rep(NA,(num_motif+1))
    for (rep_potent_wui in 1:(num_motif+1)) {
      if(rep_potent_wui == (num_motif+1)){
        llk_one_seq[rep_potent_wui] = exp(logllk_fun_0(Rui,theta_0_samp,dict))
      }else{
        aui = A_samp[[rep_potent_wui]][rep_i,]
        llk_one_seq[rep_potent_wui] = exp(logllk_fun_k(Rui,rep_potent_wui,aui,theta_0_samp,theta_samp,dict))
      }
    }
     
    post_wp = (pri_wp*llk_one_seq)/sum(pri_wp*llk_one_seq)
    W_samp_unobs[which(UW_loc==rep_i)] =  sample((1:(num_motif+1)),size = 1,replace=T,prob = post_wp)
  }
  return(W_samp_unobs)
}