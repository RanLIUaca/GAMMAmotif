#' Bayesian Inference for Motif Discovery using MCMC
#'
#' This function performs Bayesian inference for motif discovery in biological sequences
#' using Gibbs sampling and Metropolis-Hastings updates. It estimates motif-specific parameters,
#' latent motif positions, and gap distributions for multiple motifs in observed sequences.
#'
#' @param Data A data frame with two columns: \code{seq} (character vector of sequences) and \code{W_obs} (numeric vector of observed motif labels; use \code{NA} for unobserved).
#' @param dict A character vector representing the dictionary/alphabet (e.g., \code{c("A", "C", "G", "T")}).
#' @param num_motif An integer specifying the number of motifs to infer.
#' @param Motif_length An integer vector of length \code{num_motif}, specifying the length of each motif.
#' @param sum_N An integer specifying the total number of MCMC iterations to perform.
#' @param burnin An integer specifying the number of burn-in iterations (before shift sampling is activated).
#' @param mh_step_w An integer specifying the interval at which Metropolis-Hastings shift proposals are made.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{Theta_0_SAMP}}{MAP estimate of the background distribution.}
#'   \item{\code{Theta_SAMP}}{List of MAP-estimated motif PWMs.}
#'   \item{\code{W}}{MAP estimate of motif assignments.}
#'   \item{\code{A}}{List of MAP-estimated motif position matrices.}
#'   \item{\code{plot}}{List of ggplot2 objects: log-likelihood trace and sequence logos.}
#' }
#'
#' @details
#' including:
#' \itemize{
#'   \item \code{generate_aijk()}
#'   \item \code{w_samp_fun()}
#'   \item \code{theta_0_samp_fun()}
#'   \item \code{theta_samp_fun()}
#'   \item \code{A_samp_fun()}
#'   \item \code{lambda_samp_fun()}
#'   \item \code{A_theta_lambda_samp_fun()}
#'   \item \code{logllk_fun()}
#'   \item \code{GAMMA_res()}
#' }
#' These must be sourced or included in the package's \code{R/} directory.
#'
#' @import MCMCpack
#' @import stringr
#' @import MASS
#' @import ggseqlogo
#' @import ggplot2
#' @export

GAMMA_FUN = function(Data, dict,num_motif, Motif_length, sum_N, burnin, mh_step_w){
  stop_jump_step = burnin
  total_n = length(Data[,1])
  colnames(Data) = c('seq', 'W_obs')
  R = Data$seq; W_obs = Data$W_obs
  Len_seq = rep(NA, total_n)
  for (rep_n in 1:total_n) {
    Ri = R[rep_n]
    Len_seq[rep_n] = str_length(Ri)
  }


  ################################################
  ### Generate parameters and latent varibels ###
  ###############################################

  #### initial value ###
  ##prior parameters
  pri_al0=1
  pri_alw=rep(1,num_motif)
  pri_wp=c(rep((1/(num_motif+1)),(num_motif+1)))
  pri_beta = 1; pri_nu = 1

  alpha_0 = rep(pri_al0,length(dict))
  alpha = list()
  for (rep_mt in 1:num_motif) {
    alpha[[rep_mt]] = rep(pri_alw[rep_mt],length(dict))
  }


  # initial values
  theta_0_ini = c(rdirichlet(1,alpha_0))

  theta_ini = list()
  for (rep_motif in 1:num_motif) {
    motif_len_w  = Motif_length[rep_motif]
    theta_ini[[rep_motif]] = matrix(nrow = length(dict), ncol = motif_len_w)
    for (i in 1:motif_len_w) {
      temp = rdirichlet(1,alpha[[rep_motif]])
      theta_ini[[rep_motif]][,i] = temp
    }
  }


  W_ini = Data$W_obs; UW_loc = which(is.na(W_ini))
  W_ini[UW_loc] = sample((1:(num_motif+1)),size = length(UW_loc),replace=T,prob = pri_wp)

  lambda_ini = list()
  for (rep_motif in 1:num_motif) {
    motif_len_w  = Motif_length[rep_motif]
    if(motif_len_w>1) {
      lambda_ini[[rep_motif]] = rep((round(mean(Len_seq-1)/(motif_len_w-1))+.5), (motif_len_w-1))}else{
        lambda_ini[[rep_motif]] = NA
      }
    }

  A_ini = list()
  # source('FUN/FLAT_GENE_aijk_FUN.R')
  for (rep_motif in 1:num_motif) {
    # rep_motif = 1
    motif_len_w  = Motif_length[rep_motif]
    A_ini[[rep_motif]] = matrix(nrow = total_n, ncol = motif_len_w)
    for (i in 1:total_n) {
      # i = 1
      a = sample(1:(Len_seq[i]-motif_len_w+1),size = 1,replace=T)
      A_ini[[rep_motif]][i,1] = a
      if(motif_len_w>1){
        for (rep_loc_mt in 2:motif_len_w) {
          # rep_loc_mt = 3
          A_ini[[rep_motif]][i,rep_loc_mt] = generate_aijk(A_ini[[rep_motif]][i,rep_loc_mt-1],lambda_ini[[rep_motif]][motif_len_w-1],
                                                           Len_seq[i],motif_len_w,rep_loc_mt)
        }
      }
    }
  }


  Logllk_SAMP = rep(NA,sum_N)
  Theta_0_SAMP = matrix(NA, nrow=len_dict, ncol = sum_N)
  W_SAMP = matrix(NA, nrow=total_n, ncol = sum_N)

  Theta_SAMP = list(); A_SAMP = list(); Lambda_SAMP = list()
  for(rep_mt in 1:num_motif){
    motif_len_w  = Motif_length[rep_motif]
    Theta_SAMP[[rep_mt]] = rep(list(matrix(nrow = len_dict,ncol = motif_len_w)),sum_N)
    A_SAMP[[rep_mt]] = rep(list(matrix(nrow = total_n,ncol = motif_len_w)),sum_N)
    Lambda_SAMP[[rep_mt]] = rep(list(rep(NA,(motif_len_w-1))),sum_N)
  }


  for (rep_N in 1:sum_N) {
    # rep_N=1
    t1 = proc.time()
    print(rep_N)
    if(rep_N==1){
      W_samp = W_ini
      theta_0_samp = theta_0_ini

      theta_samp = theta_ini
      lambda_samp = lambda_ini
      A_samp = A_ini
    }else{
      W_samp[UW_loc] = w_samp_fun(R,UW_loc,A_samp,dict,theta_0_samp,theta_samp,pri_wp)
      theta_0_samp = theta_0_samp_fun(R,W_samp,A_samp,dict,pri_al0)

      temp_theta_samp = list()
      for (rep_mt in 1:num_motif) {
        #rep_mt = 1
        motif_len_w  = Motif_length[rep_motif]

        theta_samp[[rep_mt]] = theta_samp_fun(R,W_samp,A_samp,dict,pri_alw,rep_mt,motif_len_w)
        temp_theta_samp[[rep_mt]] = as.data.frame(theta_samp[[rep_mt]])
        row.names(temp_theta_samp[[rep_mt]]) = dict

        A_samp[[rep_mt]] = A_samp_fun(R, W_samp, lambda_samp, A_samp,dict,theta_0_samp,theta_samp, rep_mt,1:motif_len_w, 1:total_n, motif_len_w)
        if(motif_len_w>1){lambda_samp[[rep_mt]] = lambda_samp_fun(R,W_samp, lambda_samp, A_samp,pri_beta, pri_nu, rep_mt, motif_len_w)}
      }

      ### shift for theta A lambda #####
      if(((rep_N%%mh_step_w)==0)&(rep_N<stop_jump_step)){
        # print('start shift')
        A_theta_lambda_samp = list()
        for (rep_mt in 1:num_motif) {
          #rep_mt = 1
          # print(rep_mt)
          # t3 = proc.time()
          motif_len_w  = Motif_length[rep_motif]
          if((motif_len_w>1)){
          A_theta_lambda_samp[[rep_mt]] = A_theta_lambda_samp_fun(R,W_samp,lambda_samp, A_samp, dict,theta_0_samp,theta_samp, pri_alw, pri_beta, pri_nu, rep_mt, motif_len_w)
          A_samp[[rep_mt]] = A_theta_lambda_samp[[rep_mt]][[1]]
          theta_samp[[rep_mt]] = A_theta_lambda_samp[[rep_mt]][[2]]
          lambda_samp[[rep_mt]] = A_theta_lambda_samp[[rep_mt]][[3]]

          }
        }

      }

    }



    logllk = logllk_fun(R,W_samp,A_samp,dict,theta_0_samp,theta_samp)


    Logllk_SAMP[rep_N] = logllk
    Theta_0_SAMP[,rep_N] = theta_0_samp
    W_SAMP[,rep_N] = W_samp
    for(rep_mt in 1:num_motif){
      motif_len_w  = Motif_length[rep_motif]

      Theta_SAMP[[rep_mt]][[rep_N]] = theta_samp[[rep_mt]]
      A_SAMP[[rep_mt]][[rep_N]] = A_samp[[rep_mt]]
      if(motif_len_w>1) {Lambda_SAMP[[rep_mt]][[rep_N]] = lambda_samp[[rep_mt]]}
    }

  }


  Res = GAMMA_res(Logllk_SAMP, Theta_0_SAMP, W_SAMP, Theta_SAMP, A_SAMP, Lambda_SAMP, sum_N, burnin)

  return(Res)
}





