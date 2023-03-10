% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_helpers.R
\name{get_cutoffdist}
\alias{get_cutoffdist}
\title{Covariance cutoff distance}
\usage{
get_cutoffdist(fit, L, eps)
}
\arguments{
\item{fit}{LDA fit providing covariance functions of training data.}

\item{L}{Number of latent variables.}

\item{eps}{Tolerance up to which values are regarded as 0.}
}
\description{
Retrieves distance up to which covariances are > 0
}
