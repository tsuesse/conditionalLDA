#***************************************************************************
# loglik of multivariate normal----
#***************************************************************************
#' Loglikelihood
#' @description Computes loglikelihood for given parameters of a multivariate normal distribution.
#' @param y dimension reduced variable values
#' @param mu class means
#' @param C Covariance Matrix
#' @param logdet logarithm of matrix determinant

loglik <- function(y, mu, C, logdet = NULL) {
  res <- as.matrix(y - mu)

  if (is.null(logdet)) {
    cholC <- chol(C)
    logdet <- 2 * sum(log(diag(cholC)))
  }

  # logarithmised form of multivariate normal density function
  ll <- -0.5 * (t(res) %*% Matrix::solve(C, res))[1] # is it Matrix::solve?? It is definitely faster
  ll <- ll - 0.5 * (logdet + (dim(C)[1]) * log(2 * pi))[1]
}

#***************************************************************************
# loglik_train_sp----
#***************************************************************************
#' Loglikelihood training
#' @description Computes loglikelihood of training observations
#' @param fit LDA fit with covariance functions and coords of training data
#' @param L Number of latent variables
#' @param adjust_sigmas Should variances and/or covariances be adjusted? Default is 1 = adjust both, 2 = adjust only variances, NULL = no adjustments.
loglik_train_sp <- function(fit, L) {

  coords_old <- as.matrix(fit$coords)
  n_old <- dim(coords_old)[1]

  SIGMA_old_list <- list()

  for (l in 1:L) {
    SIGMA_old_list[[l]] <- geoR::varcov.spatial(
      coords = coords_old,
      cov.model = fit$variograms[[l]]$cov.model,
      kappa = fit$variograms[[l]]$kappa,
      nugget = fit$variograms[[l]]$nugget,
      cov.pars = fit$variograms[[l]]$cov.pars
    )$varcov
    SIGMA_old_list[[l]] <- Matrix::drop0(SIGMA_old_list[[l]], tol = 1e-10) # make it sparser


    # otherwise do nothing leave SIGMA.list unchanged
  } # end for (l in 1:L)

  #************
  # LL computed for joint matrix

  # combine l covariance matrizes into one block diagonal matrix
  #SIGMA_old <- Matrix::bdiag(SIGMA_old_list)

  #ind <- rep(1:n_old, each = L) + rep(n_old * (0:(L - 1)), times = n_old)
  #SIGMA_old <- SIGMA_old[ind, ind]

  # retrieve class means and dimension reduced data
  #mu_o <- c(t(fit$means[fit$classes, ]))
  #y_o <- c(t(as.matrix(fit$Y_old)))
  # compute Loglikelihood
  #LLmodel <- loglik(y_o, mu_o, SIGMA_old)

  #************
  # alternative way with LL computed separately for each l
  LLmodel <- rep(NA, times = L)
  for (l in 1:L) {
    SIGMA_old <- as.matrix(SIGMA_old_list[[l]])
    mu_o <- c(t(fit$means[l][fit$classes, ]))
    y_o <- c(t(as.matrix(fit$Y_old[ , l])))
    LLmodel[l] <- loglik(y_o, mu_o, SIGMA_old)
  }
  LLmodel <- sum(LLmodel)
  #************
}
