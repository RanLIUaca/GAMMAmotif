% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FLAT_LOGLLK_FUN.R
\name{logllk_fun_one_seq}
\alias{logllk_fun_one_seq}
\title{Log-Likelihood for a Single Sequence — Internal}
\usage{
logllk_fun_one_seq(Ri, wi, ai, dict, theta_0_samp, theta_samp)
}
\arguments{
\item{Ri}{A character vector representing a single sequence.}

\item{wi}{Integer. Motif assignment (if \code{wi == num_motif + 1}, use background model).}

\item{ai}{A numeric vector of motif positions for the sequence.}

\item{dict}{A character vector of the sequence alphabet.}

\item{theta_0_samp}{A vector of background letter probabilities.}

\item{theta_samp}{A list of motif position weight matrices.}
}
\value{
A numeric value representing the log-likelihood.
}
\description{
Computes the log-likelihood of a single sequence under either the
background model or a specific motif model, depending on its motif assignment.
}
\keyword{internal}
