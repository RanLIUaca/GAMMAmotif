% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FLAT_MH_A_THETA_LAMBDA_FUN_10.R
\name{log_prob_theta_A_fun}
\alias{log_prob_theta_A_fun}
\title{Proposal Log Probability of Theta from A — Internal}
\usage{
log_prob_theta_A_fun(
  R,
  W_samp,
  A,
  dict,
  theta,
  pri_alw,
  rep_mt,
  col_name,
  motif_len_w
)
}
\arguments{
\item{R}{A character vector of sequences.}

\item{W_samp}{A numeric vector of motif assignments.}

\item{A}{A matrix of motif positions for a given motif.}

\item{dict}{Alphabet vector.}

\item{theta}{PWM matrix for the motif.}

\item{pri_alw}{Dirichlet prior concentration vector.}

\item{rep_mt}{Integer. Index of motif.}

\item{col_name}{Integer vector of motif positions (columns) being updated.}
}
\value{
Log proposal probability of theta.
}
\description{
Computes the log probability of the proposed theta columns using posterior Dirichlet
distribution based on current motif positions A.
}
\keyword{internal}
