#' Posterior Prediction Analysis for Motif Discovery
#'
#' This function performs posterior predictive analysis using sampled parameters from
#' a Bayesian motif discovery model. It estimates the most probable motif assignments
#' and positions given posterior samples of background and motif-specific distributions.
#'
#' @param Data A data frame with two columns:
#'   \describe{
#'     \item{\code{seq}}{A character vector of sequences.}
#'     \item{\code{W_obs}}{A numeric vector of observed motif labels (use \code{NA} for unobserved).}
#'   }
#' @param theta_0_samp A numeric matrix of sampled background distributions across MCMC iterations.
#' @param theta_samp A list of sampled motif PWMs (position weight matrices) across iterations.
#' @param Motif_length An integer vector specifying the length of each motif.
#' @param burnin An integer specifying the number of burn-in iterations to discard.
#' @param sum_N An integer specifying the total number of MCMC iterations used.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{logllk_plot}}{A \code{ggplot2} object showing the log-likelihood trace across MCMC iterations.}
#'   \item{\code{est_W}}{MAP estimate of motif assignments for each sequence.}
#'   \item{\code{est_A}}{A list of MAP-estimated motif position matrices for each motif.}
#' }
#'
#' @details
#' This function assumes that MCMC sampling has already been performed using a function like \code{GAMMA_FUN()},
#' and that posterior samples of model parameters (\code{theta_0_samp} and \code{theta_samp}) are available.
#'
#' Several helper functions must be available in the environment or sourced before using this function, including:
#' \itemize{
#'   \item \code{generate_aijk()}
#'   \item \code{w_samp_fun()}
#'   \item \code{A_samp_fun()}
#'   \item \code{lambda_samp_fun()}
#'   \item \code{A_theta_lambda_samp_fun()}
#'   \item \code{logllk_fun()}
#' }
#'
#' @import ggplot2
#' @import stringr
#' @export
pred_ana = function(Data,theta_0_samp,theta_samp,Motif_length, burnin, sum_N){
  ###############################################
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
  W_SAMP = matrix(NA, nrow=total_n, ncol = sum_N)

  A_SAMP = list(); Lambda_SAMP = list()
  for(rep_mt in 1:num_motif){
    motif_len_w  = Motif_length[rep_motif]
    A_SAMP[[rep_mt]] = rep(list(matrix(nrow = total_n,ncol = motif_len_w)),sum_N)
    Lambda_SAMP[[rep_mt]] = rep(list(rep(NA,(motif_len_w-1))),sum_N)
  }


  for (rep_N in 1:sum_N) {
    # rep_N=1
    t1 = proc.time()
    print(rep_N)
    if(rep_N==1){
      W_samp = W_ini

      lambda_samp = lambda_ini
      A_samp = A_ini
    }else{
      W_samp[UW_loc] = w_samp_fun(R,UW_loc,A_samp,dict,theta_0_samp,theta_samp,pri_wp)

      for (rep_mt in 1:num_motif) {
        #rep_mt = 1
        motif_len_w  = Motif_length[rep_motif]
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
    W_SAMP[,rep_N] = W_samp
    for(rep_mt in 1:num_motif){
      motif_len_w  = Motif_length[rep_motif]

      A_SAMP[[rep_mt]][[rep_N]] = A_samp[[rep_mt]]
      if(motif_len_w>1) {Lambda_SAMP[[rep_mt]][[rep_N]] = lambda_samp[[rep_mt]]}
    }

  }


  ############## Define the Estimation for Theta_SAMP ################
  ### define the estimation as MAP ###
  loc_maxmizer = burnin + which.max(Logllk_SAMP[((burnin+1):sum_N)])
  map_W = W_SAMP[,loc_maxmizer]
  map_A = list()

  for (rep_mt in 1:num_motif) {
    map_A[[rep_mt]] =  A_SAMP[[rep_mt]][[loc_maxmizer]]
  }

  #####################
  ### plot #############
  #######################
  plot_list = list()

  ##################################
  ############# logllk curve ###########
  ###############################
  data_logllk = data.frame(x = 1:length(Logllk_SAMP),value = Logllk_SAMP)

  # ????????ͼ
  logllk_plot = ggplot(data_logllk, aes(x = x, y = value)) +
    geom_line(color = "#3498DB") +
    labs(title = "Loglikelihood", x = " ", y = " ") +
    theme_bw()+theme_classic()+
    theme(panel.grid=element_blank(),plot.title = element_text(size = 11),
          axis.text=element_text(size=11),
          legend.key = element_blank(),
          legend.title = element_blank(),
          legend.position = "none"
    )
  plot_list [[1]] = logllk_plot


  return(list(logllk_plot = logllk_plot, est_W = map_W, est_A = map_A))
}





