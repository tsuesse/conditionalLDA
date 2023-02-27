#***************************************************************************
# adj_sigmas ----
#***************************************************************************
# function to adjust variances and/or covariances
#' Adjust Covariance Matrizes
#' @description Helper function that adjusts only variances or both variances and covariances.
#' @param adjust_sigmas Should variances and/or covariances be adjusted? Default is 1 = adjust only variances, 2 = adjust both, NULL = no adjustments.
#' @param fit LDA fit with covariance functions
#' @param n number of cases
#' @param l index of variable for which covariance matrix should be adjusted
#' @param SIGMA_list List containing covariance matrizes for L
#' @param trace Controls console output. NULL = no output, 1 = basic output,
#'  2 = detailed output.
adj_sigmas <- function(adjust_sigmas = 1, fit, n, l, SIGMA_list, trace = 1) {


  # print(length(diag(SIGMA_list[[l]])))

  if (adjust_sigmas == 1) {
    # ADJUST just the variances

    if (fit$variograms[[l]]$cov.pars[1] < 1) {
      diag(SIGMA_list[[l]]) <- sapply(rep(1, n),
                                      FUN = f_sigmas,
                                      simplify = T, covs = fit$covs, j = l
      )

      if (trace > 1) {
        cat("Variable =", l, "Sigma < 1 \n")
      }
    } else {
      SIGMA_list[[l]] <- Matrix::Diagonal(n, rep(1, n))

      if (trace > 1) {
        cat("Variable =", l, "Sigma > 1 \n")
      }
    }

  } else {
    if (adjust_sigmas == 2) {
      # ADJUST variances and covariances

      hilf <- as.matrix(SIGMA_list[[l]])

      sigmas <- sqrt(sapply(rep(1, n),
                            FUN = f_sigmas, simplify = T,
                            covs = fit$covs, j = l
      )) / sqrt(diag(hilf))

      # print(dim(SIGMA_list[[l]]))

      # hence we must multiply by sigma again
      SIGMA_list[[l]] <- Matrix::Diagonal(n, sigmas) %*% SIGMA_list[[l]] %*% Matrix::Diagonal(n, sigmas)

    } else {
      cat("No sigma adjustments applied.For adjusting only variances set
          adjust_sigmas = 1, for adjusting both covariances and variance
          set adjust_sigmas = 2")
    }
  }
  return(SIGMA_list)
}

#***************************************************************************
#  calc_SIGMA_cond ----
#***************************************************************************
#' Conditional Covariance Matrix
#' @description Computes a covariance matrix conditional on another subset of
#' the same multivariate normal distribution.
#' @param ind_new Index in combined block diagonal Matrix to identify new, i.e.
#' test cases
#' @param n_all_j Number of all included cases: test cases to be predicted,
#' training cases included in conditional prediction
#' @param L Number of latent variables
#' @param SIGMA Combined block diagonal covariance Matrix
#' @param SIGMA_nn Subset of SIGMA, contains only values for test cases to be
#' predicted
calc_SIGMA_cond <- function(n_new_j, n_all_j, SIGMA) {
  ind_old <- (n_new_j + 1):n_all_j
  ind_new <- 1:n_new_j
  SIGMA_nn <- SIGMA[ind_new, ind_new]
  SIGMA_no <- SIGMA[ind_new, ind_old]
  SIGMA_oo <- SIGMA[ind_old, ind_old]
  SIGMA_oo_inv_SIGMA_on <- try(Matrix::solve(SIGMA_oo, Matrix::t(SIGMA_no)))

  hilf <- try(SIGMA_no %*% SIGMA_oo_inv_SIGMA_on)
  hilf <- Matrix::forceSymmetric(hilf)
  SIGMA_cond <- SIGMA_nn - hilf
  return(tibble::lst(SIGMA_oo_inv_SIGMA_on, SIGMA_cond,SIGMA_no))
}
#***************************************************************************
# small functions ----
#***************************************************************************
# helper functions to retrieve covariances in different levels of a fit
#' Get Covariances1
#' @description Retrieves covariances
#' @param x index
#' @param covs Object to retrieve covariances from
f_covs <- function(x, covs) {
  return(covs[[x]])
}
#' Get Covariances2
#' @description Retrieves covariances for a certain variable
#' @param j Index of block/case to be predicted
#' @inheritParams f_covs
f_sigmas <- function(x, covs, j) {
  return(covs[[x]][j, j])
}
