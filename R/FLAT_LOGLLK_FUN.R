#' Background Model Log-Likelihood — Internal
#'
#' Computes the log-likelihood of a sequence under the background model 
#' (i.e., no motif), using the background letter distribution.
#'
#' @param Ri A character vector representing a single sequence.
#' @param theta_0_samp A vector of background letter probabilities.
#' @param dict A character vector of the sequence alphabet (e.g., \code{c("A", "C", "G", "T")}).
#'
#' @return A numeric value representing the log-likelihood.
#' @keywords internal
logllk_fun_0 = function(Ri,theta_0_samp,dict){
  Ri_factor = factor(Ri, levels = dict)
  res = sum((table(Ri_factor))*log(theta_0_samp))
  return(res)
}

#' Log-Likelihood under Motif Model — Internal
#'
#' Computes the log-likelihood of a sequence under a motif model, given the motif 
#' assignment and motif position indices.
#'
#' @param Ri A character vector representing a single sequence.
#' @param wi Integer. The motif index assigned to the sequence.
#' @param ai A numeric vector of motif start positions.
#' @param theta_0_samp A vector of background letter probabilities.
#' @param theta_samp A list of motif position weight matrices.
#' @param dict A character vector of the sequence alphabet.
#'
#' @return A numeric value representing the log-likelihood.
#' @keywords internal
logllk_fun_k = function(Ri,wi,ai,theta_0_samp,theta_samp,dict){
  Ri_w = Ri[ai]; Ri_els = Ri[-ai]
  Ri_els_factor = factor(Ri_els, levels = dict)
  res = sum((table(Ri_els_factor))*log(theta_0_samp))
  for (rep_j in 1:length(ai)) {
    Ri_aij =  Ri_w[rep_j]
    Ri_aij_factor = factor(Ri_aij, levels = dict)
    res = res+sum((table(Ri_aij_factor))*log(theta_samp[[wi]][,rep_j]))
  }
  return(res)
}

#' Log-Likelihood for Motif Position j — Internal
#'
#' Computes the partial log-likelihood of a sequence under the motif model, 
#' focusing on a specific motif position \code{j}, while treating the rest as background.
#'
#' @param Ri A character vector representing a single sequence.
#' @param wi Integer. The motif index assigned to the sequence.
#' @param col_j Integer. The column index (motif position) to evaluate.
#' @param ai A numeric vector of motif start positions.
#' @param theta_0_samp A vector of background letter probabilities.
#' @param theta_samp A list of motif position weight matrices.
#' @param dict A character vector of the sequence alphabet.
#'
#' @return A numeric value representing the partial log-likelihood.
#' @keywords internal
logllk_fun_kj = function(Ri,wi,col_j,ai,theta_0_samp,theta_samp,dict){
  Ri_w = Ri[ai]; Ri_els = Ri[-ai]
  Ri_els_factor = factor(Ri_els, levels = dict)
  res = sum((table(Ri_els_factor))*log(theta_0_samp))
  
  Ri_aij =  Ri_w[col_j]
  Ri_aij_factor = factor(Ri_aij, levels = dict)
  res = res+sum((table(Ri_aij_factor))*log(theta_samp[[wi]][,col_j]))

  return(res)
}



#' Log-Likelihood for a Single Sequence — Internal
#'
#' Computes the log-likelihood of a single sequence under either the 
#' background model or a specific motif model, depending on its motif assignment.
#'
#' @param Ri A character vector representing a single sequence.
#' @param wi Integer. Motif assignment (if \code{wi == num_motif + 1}, use background model).
#' @param ai A numeric vector of motif positions for the sequence.
#' @param dict A character vector of the sequence alphabet.
#' @param theta_0_samp A vector of background letter probabilities.
#' @param theta_samp A list of motif position weight matrices.
#'
#' @return A numeric value representing the log-likelihood.
#' @keywords internal
logllk_fun_one_seq = function(Ri,wi,ai,dict,theta_0_samp,theta_samp){
  if((wi==(num_motif+1))){
    res=logllk_fun_0(Ri,theta_0_samp,dict)
  }else if((wi!=(num_motif+1))){
    res=logllk_fun_k(Ri,wi,ai,theta_0_samp,theta_samp,dict)
  }
  return(res)
}

#' Total Log-Likelihood over All Sequences — Internal
#'
#' Computes the total log-likelihood across all sequences, under either 
#' the background or motif model depending on each sequence's assignment.
#'
#' @param R A character vector of sequences.
#' @param W_samp A numeric vector of motif assignments.
#' @param A_samp A list of position matrices for each motif.
#' @param dict A character vector of the sequence alphabet.
#' @param theta_0_samp A vector of background letter probabilities.
#' @param theta_samp A list of motif position weight matrices.
#'
#' @return A numeric value representing the total log-likelihood.
#' @keywords internal
logllk_fun = function(R,W_samp,A_samp,dict,theta_0_samp,theta_samp){
  res=0
  total_n=length(R)
  for (rep_i in 1:total_n) {
    #rep_i=1
    Ri = R[rep_i]
    tmp_index = 1:str_length(Ri)
    Ri = str_sub(Ri,tmp_index,tmp_index)
    wi = W_samp[rep_i]

    if((wi==(num_motif+1))){
      res=res+logllk_fun_0(Ri,theta_0_samp,dict)
    }else if((wi!=(num_motif+1))){
      ai = A_samp[[wi]][rep_i,]
      res=res+logllk_fun_k(Ri,wi,ai,theta_0_samp,theta_samp,dict)
      # ai = A_samp[[1]][rep_i,]; logllk_fun_k(Ri,1,ai,theta_0_samp,theta_samp,dict)
    }
    # print(c(rep_i,res))
  }
  return(res)
}