#' Joint Log Posterior of A, Theta, and Lambda — Internal
#'
#' Computes the joint log posterior for a given motif, combining:
#' - sequence log-likelihood,
#' - Dirichlet priors on theta,
#' - Gamma priors on lambda,
#' - Poisson gap likelihood for motif positions.
#'
#' @param R A character vector of sequences.
#' @param W_samp A numeric vector of motif assignments.
#' @param lambda A numeric vector of Poisson gap parameters for the motif.
#' @param A A list of motif position matrices (motif-wise).
#' @param dict A character vector of the alphabet (e.g., \code{c("A", "C", "G", "T")}).
#' @param theta_0_samp Background distribution.
#' @param theta A list of motif position weight matrices.
#' @param pri_alw A numeric vector of Dirichlet prior concentration parameters.
#' @param pri_beta Gamma prior shape parameter for lambda.
#' @param pri_nu Gamma prior rate parameter for lambda.
#' @param rep_mt Integer. Index of the motif being evaluated.
#'
#' @return A numeric value: the joint log posterior.
#' @keywords internal
log_pi_A_theta_lambda_fun = function(R,W_samp,lambda, A, dict,theta_0_samp,theta,pri_alw, pri_beta, pri_nu, rep_mt, motif_len_w){
  # A = A_star; theta = theta_star; lambda = lambda_star
  res = logllk_fun(R,W_samp,A,dict,theta_0_samp,theta)
  for(rep_j in 1:motif_len_w){
    res =(res +log(ddirichlet(theta[[rep_mt]][,rep_j], rep(pri_alw[rep_mt],len_dict))))
  }

  for(rep_j in 1:(motif_len_w-1)){
    res =(res + log(dgamma(lambda[rep_j], shape = pri_beta, rate = pri_nu)))
  }


  for (rep_i in 1:total_n) {
    # rep_i=2
    Ri = R[rep_i]
    len_seq = str_length(Ri)
    tmp_index = 1:str_length(Ri)
    Ri = str_sub(Ri,tmp_index,tmp_index)
    wi = W_samp[rep_i]

    if(wi == rep_mt){
      for (rep_j in 2:motif_len_w) {
        #rep_j=1
        Range = (A[[rep_mt]][rep_i,rep_j-1]+1):(len_seq-motif_len_w + rep_j)

        aij = A[[rep_mt]][rep_i,rep_j]
        aij_pre = A[[rep_mt]][rep_i,rep_j-1]

        if(length(Range)>1){
          max_y = len_seq - motif_len_w + (rep_j) - (aij_pre) - 1
          # 鍒嗘瘝閮ㄥ垎 (褰掍竴鍖栧父鏁?)
          y_vals = 0:max_y
          poisson_probs = dpois(y_vals, lambda[rep_j-1])  # 娉婃澗鍒嗗竷姒傜巼
          normalization_constant = sum(poisson_probs)  # 鍒嗘瘝
          shifted_probs = (poisson_probs / normalization_constant)  # 褰掍竴鍖栫殑姒傜巼
          res = (res + log(poisson_probs[aij-aij_pre])-
                   log(normalization_constant))
        }
      }
    }
  }

  return(res)
}

#' Proposal Log Probability of A — Internal
#'
#' Computes the log probability of a proposed A (motif positions) under a uniform
#' proposal distribution for the selected column.
#'
#' @param R A character vector of sequences.
#' @param W_samp Motif assignments.
#' @param lambda Lambda vector for the motif.
#' @param A A list of motif position matrices.
#' @param dict Alphabet vector.
#' @param theta_0_samp Background model.
#' @param theta List of motif PWMs.
#' @param rep_mt Integer. Index of motif.
#' @param col_name Integer. Column index (i.e., motif position) being proposed.
#'
#' @return Log proposal probability.
#' @keywords internal
log_prob_A_lambda_theta_fun = function(R, W_samp, lambda, A, dict,theta_0_samp,theta, rep_mt,col_name, motif_len_w){
  # A= A_star; theta = theta_star
  res = 0
  for (rep_i in 1:total_n) {
    # rep_i=2;  rep_mt = 1
    Ri = R[rep_i]
    len_seq = str_length(Ri)
    tmp_index = 1:str_length(Ri)
    Ri = str_sub(Ri,tmp_index,tmp_index)
    wi = W_samp[rep_i]

    if(wi == rep_mt){
      #col_name=1
      if(col_name == 1){
        Range = 1:(A[[rep_mt]][rep_i,col_name+1]-1)
      }else if((col_name > 1)&(col_name < motif_len_w)){
        Range = (A[[rep_mt]][rep_i,col_name-1]+1):(A[[rep_mt]][rep_i,col_name+1]-1)
      }else{
        Range = (A[[rep_mt]][rep_i,col_name-1]+1):len_seq
      }

      if(length(Range)>1){
        res = res + log(1/length(Range))
      }
    }
  }
  return(res)
}

#' Proposal Log Probability of Theta from A — Internal
#'
#' Computes the log probability of the proposed theta columns using posterior Dirichlet
#' distribution based on current motif positions A.
#'
#' @param R A character vector of sequences.
#' @param W_samp A numeric vector of motif assignments.
#' @param A A matrix of motif positions for a given motif.
#' @param dict Alphabet vector.
#' @param theta PWM matrix for the motif.
#' @param pri_alw Dirichlet prior concentration vector.
#' @param rep_mt Integer. Index of motif.
#' @param col_name Integer vector of motif positions (columns) being updated.
#'
#' @return Log proposal probability of theta.
#' @keywords internal
log_prob_theta_A_fun = function(R,W_samp,A,dict,theta,pri_alw,rep_mt,col_name, motif_len_w){
  # sample parameter
  # A = A_star[[rep_mt]]; theta = theta_star[[rep_mt]]
  Post_alpha_w = matrix(nrow = length(dict), ncol = motif_len_w)
  for (rep_j in col_name) {
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
        aij = A[rep_i,rep_j]
        Ri_aij = Ri[aij]
        tmp_index = 1:str_length(Ri_aij)
        Ri_aij = str_sub(Ri_aij,tmp_index,tmp_index)
        Ri_aij = factor(Ri_aij, levels = dict)
        H_a=H_a+table(Ri_aij)
      }
    }

    post_alpha_w = pri_alw[rep_mt]+H_a
    Post_alpha_w[,rep_j] = post_alpha_w
  }

  res = 0
  for (rep_j in col_name) {
    res = res + log(ddirichlet(theta[,rep_j], Post_alpha_w[,rep_j]))
  }
  return(res)
}



#' Metropolis-Hastings Sampling of A, Theta, Lambda — Internal
#'
#' Proposes a joint update (block shift) of A (motif positions), Theta (motif PWM),
#' and Lambda (Poisson gap parameters) for a selected motif using a Metropolis-Hastings step.
#'
#' The method:
#' - selects a block of consecutive motif positions,
#' - shifts the block left or right,
#' - proposes new parameters,
#' - evaluates acceptance based on joint posterior and proposal probabilities.
#'
#' @param R A character vector of sequences.
#' @param W_samp Motif assignments.
#' @param lambda_samp A list of current lambda vectors.
#' @param A_samp A list of current motif position matrices.
#' @param dict Alphabet vector.
#' @param theta_0_samp Background distribution.
#' @param theta_samp A list of PWM matrices.
#' @param pri_alw Dirichlet prior concentration vector.
#' @param pri_beta Gamma prior shape for lambda.
#' @param pri_nu Gamma prior rate for lambda.
#' @param rep_mt Integer. Index of motif being updated.
#'
#' @return A list with updated \code{A}, \code{theta}, \code{lambda}, and a flag \code{no_jump}
#' indicating whether the proposal was accepted.
#' @keywords internal
A_theta_lambda_samp_fun = function(R,W_samp,lambda_samp, A_samp, dict,theta_0_samp,theta_samp, pri_alw, pri_beta, pri_nu, rep_mt, motif_len_w){
  # rep_mt = 1
  candi_sf_len = 1:(motif_len_w-1)
  jumped = FALSE  # Flag to check if we need to exit the outer loop

  for (rep_sf_len in candi_sf_len) {
    if (jumped) break  # If a jump has occurred, exit the outer loop
    ## Define the shift blog length (all candicates of length need to be run)
    sf_len = rep_sf_len

    ## Choose the shift blog start column index ramdomly
    candi_sf_start = 1:(motif_len_w-sf_len + 1)
    if(length(candi_sf_start) == 1){
      sf_start = candi_sf_start
    }else{
      sf_start = sample(1:(motif_len_w-sf_len + 1),size = 1,replace=T,prob = rep(1/(motif_len_w-sf_len+1),(motif_len_w-sf_len+1)))
    }

    ## Define the shift blog direction for moving (all candicates of length need to be run)
    for (rep_delta in c(-1,1)) {
      if (jumped) break  # If a jump has occurred, exit the inner loop
      delta = rep_delta
      if(delta == 1){new_col = sf_start; old_col =  sf_start  + sf_len }else if(delta == -1){
        new_col = sf_start  + sf_len -1; old_col =  sf_start - 1 }

      # delta = -1
      A_star = A_samp; A_old = A_samp; A_temp = A_samp
      theta_star = theta_samp; theta_old = theta_samp
      lambda_star = lambda_samp[[rep_mt]]; lambda_old = lambda_samp[[rep_mt]]

      nume = 1; deno = 1
      if(((delta == -1)&(sf_start > 1))||((delta == 1)&(sf_start < (motif_len_w-sf_len+1)))){
        for (rep_i in 1:total_n) {
          # print(rep_i)
          # rep_i = 32
          Ri = R[rep_i]
          len_seq = str_length(Ri)
          wi = W_samp[rep_i]

          if((  wi == rep_mt )) {
            if(delta == 1){
              now_pt = A_temp[[rep_mt]][rep_i,sf_start]
              if(sf_start == 1){
                pre_pt = 0
              }else if(sf_start > 1){
                pre_pt = A_temp[[rep_mt]][rep_i,sf_start-1]
              }
              if(now_pt >(pre_pt+1)){

                A_star[[rep_mt]][rep_i,((sf_start+1):(sf_start+sf_len))] = A_star[[rep_mt]][rep_i,(sf_start:(sf_start+sf_len-1))];
                if(new_col == 1) {
                  loc_range = 1:(A_star[[rep_mt]][rep_i,(new_col+1)]-1)
                }else{
                  loc_range = (A_star[[rep_mt]][rep_i,(new_col-1)]+1):(A_star[[rep_mt]][rep_i,(new_col+1)]-1)
                }

                if(length(loc_range)==1){A_star[[rep_mt]][rep_i,new_col] = loc_range}else{
                  A_star[[rep_mt]][rep_i,new_col] = sample(loc_range,size = 1,replace=T,prob = rep(1/length(loc_range),length(loc_range)))
                }

              }
            }else if(delta == -1){
              now_pt = A_temp[[rep_mt]][rep_i,(sf_start+sf_len-1)]
              if(sf_start == (motif_len_w-sf_len+1)){
                post_pt = len_seq + 1
              }else if(sf_start < (motif_len_w-sf_len+1)){
                post_pt = A_temp[[rep_mt]][rep_i,(sf_start+sf_len)]
              }
              if(now_pt <(post_pt-1)){

                A_star[[rep_mt]][rep_i,((sf_start-1):(sf_start+sf_len-2))] = A_star[[rep_mt]][rep_i,(sf_start:(sf_start+sf_len-1))];
                if(new_col == motif_len_w) {
                  loc_range = (A_star[[rep_mt]][rep_i,(new_col-1)]+1):len_seq
                }else{
                  loc_range = (A_star[[rep_mt]][rep_i,(new_col-1)]+1):(A_star[[rep_mt]][rep_i,(new_col+1)]-1)
                }
                if(length(loc_range)==1){A_star[[rep_mt]][rep_i,new_col] = loc_range}else{
                  A_star[[rep_mt]][rep_i,new_col] = sample(loc_range,size = 1,replace=T,prob = rep(1/length(loc_range),length(loc_range)))
                }
              }
            }
          }
        }


        theta_star[[rep_mt]][,((sf_start+delta):(sf_start+sf_len-1+delta))] = theta_star[[rep_mt]][,(sf_start:(sf_start+sf_len-1))]
        theta_star[[rep_mt]][,new_col] = theta_samp_fun(R,W_samp,A_star,dict,pri_alw,rep_mt, motif_len_w)[,new_col]

        for(rep_j in 1:(motif_len_w-1)){
         log_lambda_old = log(lambda_old[rep_j])
         log_lambda_star = log_lambda_old + rnorm(1,0,.1)
         lambda_star[rep_j] = exp(log_lambda_star)
        }
        # browser()
        nume = (log_pi_A_theta_lambda_fun(R,W_samp,lambda_star, A_star, dict,theta_0_samp,theta_star,pri_alw, pri_beta, pri_nu, rep_mt,motif_len_w)
                + log_prob_A_lambda_theta_fun(R, W_samp, lambda_star, A_old, dict,theta_0_samp,theta_star, rep_mt,old_col, motif_len_w)
                + log_prob_theta_A_fun(R,W_samp,A_old[[rep_mt]],dict,theta_old[[rep_mt]],pri_alw,rep_mt,old_col, motif_len_w)
        )

        deno = (log_pi_A_theta_lambda_fun(R,W_samp,lambda_old, A_old, dict,theta_0_samp,theta_old,pri_alw, pri_beta, pri_nu, rep_mt, motif_len_w)
                + log_prob_A_lambda_theta_fun(R, W_samp, lambda_old, A_star, dict,theta_0_samp,theta_old, rep_mt,new_col, motif_len_w)
                + log_prob_theta_A_fun(R,W_samp,A_star[[rep_mt]],dict,theta_star[[rep_mt]],pri_alw,rep_mt,new_col, motif_len_w)
        )
      }

      # calculate acceptance rate
      acc_rate = min(1, exp(nume-deno))
      # generate a random number
      acc_rn = runif(1)

      no_jump=0;
      if((acc_rn<acc_rate)&(nume!=deno)){
        # print(c(paste0('Jump for no.motif=', rep_mt)));
        # print(paste0('theta_A jump from ', paste((sf_start:(sf_start  + sf_len -1)),collapse = ','), ' to ', paste(((sf_start+delta):(sf_start  + sf_len -1+delta)),collapse = ',')))
        # print(c(round(acc_rn,3),round(acc_rate,3)))
        no_jump=1;
        jumped = TRUE  # Set flag to exit outer loop
        break  # Exit inner loop
      }
    }
  }


  if((acc_rn<acc_rate)&(nume!=deno)){
    # browser()
    return(list(A_star[[rep_mt]],theta_star[[rep_mt]],lambda_star,no_jump))
  }else{return(list(A_old[[rep_mt]],theta_old[[rep_mt]],lambda_old,no_jump))}
}
