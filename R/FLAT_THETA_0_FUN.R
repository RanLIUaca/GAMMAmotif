#' Count Letter Frequencies for Background — Internal
#'
#' Counts the frequencies of each letter in a sequence, assuming the whole sequence
#' is generated from the background model (no motif).
#'
#' @param Ri A character vector representing a sequence.
#' @param dict A character vector representing the alphabet (e.g., \code{c("A", "C", "G", "T")}).
#'
#' @return A named table with counts of each letter in \code{dict}.
#' @keywords internal
h_0 = function(Ri,dict){
  Ri_factor = factor(Ri, levels = dict)
  return(table(Ri_factor))
}

#' Count Background Letter Frequencies Excluding Motif — Internal
#'
#' Counts the frequencies of letters in a sequence excluding the motif positions.
#' This is used to update the background model when a motif is present.
#'
#' @param Ri A character vector representing a sequence.
#' @param ai A numeric vector of motif position indices.
#' @param dict A character vector representing the alphabet.
#'
#' @return A named table with counts of background letters.
#' @keywords internal
h_k = function(Ri,ai,dict){
  Ri_els = Ri[-ai]
  Ri_els_factor = factor(Ri_els, levels = dict)
  return(table(Ri_els_factor))
}


#' Sample Background Letter Distribution Theta_0 — Internal
#'
#' Computes the posterior Dirichlet parameters for the background distribution
#' based on sequences either fully background or partially background (excluding motifs),
#' and samples \code{theta_0} from the posterior.
#'
#' @param R A character vector of sequences.
#' @param W_samp A numeric vector of motif assignments for each sequence.
#' @param A_samp A list of motif position matrices.
#' @param dict A character vector of the alphabet.
#' @param pri_al0 A numeric vector of Dirichlet prior parameters.
#'
#' @return A numeric vector sampled from the posterior Dirichlet distribution for \code{theta_0}.
#' @keywords internal
theta_0_samp_fun = function(R,W_samp,A_samp,dict,pri_al0){
  total_n=length(R)
  # post_parameter
  H_a=rep(0,length(dict))
  names(H_a)=dict
  for (rep_i in 1:total_n) {
    Ri = R[rep_i]
    tmp_index = 1:str_length(Ri)
    Ri = str_sub(Ri,tmp_index,tmp_index)
    wi = W_samp[rep_i]
    
    if((wi==(num_motif+1))){
      H_a=H_a+h_0(Ri,dict)
    }else if((wi!=(num_motif+1))){
      ai = A_samp[[wi]][rep_i,]
      H_a=H_a+h_k(Ri,ai,dict)
    }
  }
  post_alpha_0 = pri_al0+H_a
  
  # sample parameter
  theta_0_samp = c(rdirichlet(1,post_alpha_0))
  
  return(theta_0_samp)
}