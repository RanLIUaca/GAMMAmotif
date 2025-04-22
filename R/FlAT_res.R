#' Extract MAP Estimates and Plot Results from MCMC Sampling — Internal
#'
#' Extracts Maximum A Posteriori (MAP) estimates for parameters from MCMC samples,
#' including background distribution (\code{Theta_0_SAMP}), motif assignments (\code{W_SAMP}),
#' motif position weight matrices (\code{Theta_SAMP}), and motif positions (\code{A_SAMP}).
#' Also generates a list of plots including:
#' - Log-likelihood trace plot
#' - Sequence logo plots for each motif
#'
#' @param Logllk_SAMP A numeric vector of log-likelihood values from MCMC iterations.
#' @param Theta_0_SAMP A matrix of sampled background distributions (columns = iterations).
#' @param W_SAMP A matrix of motif assignment samples (columns = iterations).
#' @param Theta_SAMP A nested list of motif PWM samples.
#' @param A_SAMP A nested list of motif position matrices.
#' @param Lambda_SAMP A list of lambda samples (not used directly in this function).
#' @param sum_N Total number of MCMC iterations.
#' @param burnin Number of burn-in iterations to discard before MAP extraction.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{Theta_0}}{MAP estimate of the background distribution.}
#'   \item{\code{Theta}}{List of MAP-estimated motif PWMs.}
#'   \item{\code{W}}{MAP estimate of motif assignments.}
#'   \item{\code{A}}{List of MAP-estimated motif position matrices.}
#'   \item{\code{plot}}{List of ggplot2 objects: log-likelihood trace and sequence logos.}
#' }
#' @keywords internal
GAMMA_res = function(Logllk_SAMP, Theta_0_SAMP, W_SAMP, Theta_SAMP, A_SAMP, Lambda_SAMP, sum_N, burnin)  {
  ############## Define the Estimation for Theta_SAMP ################
  ### define the estimation as MAP ###
  loc_maxmizer = burnin + which.max(Logllk_SAMP[((burnin+1):sum_N)])
  map_Theta_0_SAMP = Theta_0_SAMP[,loc_maxmizer]
  map_W = W_SAMP[,loc_maxmizer]
  map_Theta_SAMP = list()
  map_A = list()

  for (rep_mt in 1:num_motif) {
    map_Theta_SAMP[[rep_mt]] = as.matrix(Theta_SAMP[[rep_mt]][[loc_maxmizer]])
    row.names(map_Theta_SAMP[[rep_mt]]) = dict

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

  #########################
  ### parameters logo ####
  ##########################
  start_plot = 1
  for(rep_mt in 1:num_motif){
    # rep_mt = 2
    motif_len_w  = Motif_length[rep_mt]
    if(motif_len_w>1) {start_plot = start_plot + 1
    plot_list[[start_plot]] = ggseqlogo(map_Theta_SAMP[[rep_mt]]) +
      theme(legend.position = "none", axis.text = element_text(size = 11),plot.title = element_text(size = 11))+
      labs(title = paste0('No.', rep_mt, " Motif Estimation"))
    }
  }
  return(list(Theta_0 = map_Theta_0_SAMP, Theta = map_Theta_SAMP, W = map_W, A = map_A, plot = plot_list))
}













