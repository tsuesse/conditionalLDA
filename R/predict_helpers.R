#***************************************************************************
# get_cutoffdist ----
#***************************************************************************
#' Covariance cutoff distance
#' @description Retrieves distance up to which covariances are > 0
#' @param fit LDA fit providing covariance functions of training data.
#' @param L Number of latent variables.
#' @param eps Tolerance up to which values are regarded as 0.

get_cutoffdist <- function(fit, L, eps) {
  cut_offs <- rep(0, times = L)

  # determine the distance tresholds for all models
  for (l in 1:L) {

    # retrieve covariances with corresponding separation distance
    f <- function(x) {
      return(abs(geoR::cov.spatial(x,
                                   cov.model = fit$variograms[[l]]$cov.model,
                                   cov.pars = fit$variograms[[l]]$cov.pars,
                                   kappa = fit$variograms[[l]]$kappa
      ) - eps))
    }

    # find distance with minimal covariance within 0 up to 20% of the max. distance
    cut_offs[l] <- optimise(f, interval = c(0, fit$variograms[[l]]$max.dist
                                            / 5))$minimum
  } # end for (l in 1:L)

  cut_off <- max(cut_offs)
}

#***************************************************************************
# get_nearby_obs ----
#***************************************************************************
#' Get nearby observations
#' @description Identify training observations that lie within a given distance
#' treshold.
#' @param fit LDA fit providing coords of training data.
#' @param L Number of latent variables.
#' @param index indicates the cases belonging to the block to be predicted.
#' @param n_new_j number of cases belonging to the block to be predicted.
#' @param Y_new Latent variables/predictors of test data.
#' @param coords_old Coordinates of training data.
#' @param coords_new Coordinates of test set/cases to be classified.
#' @param cut_off Distance threshold up to which training observations are included.

get_nearby_obs <- function(fit, L, index, n_new_j, Y_new, coords_old, coords_new, cut_off) {
  if (n_new_j > 1) {
    Y_new_j <- Y_new[index, ]
    coords_new_j <- coords_new[index, ]
  } else {
    Y_new_j <- matrix(Y_new[index, ], 1, L) # L dimension of transformed Y
    coords_new_j <- matrix(coords_new[index, ], 1, 2) # always 2 coordinates
  }


  # calculate distance between block and training set
  dist_old_new <- proxy::dist(coords_old, coords_new_j, method = "euclidean")
  min_dist <- apply(dist_old_new, 1, min)

  # retrieve training data within distance threshold
  index_old <- min_dist < cut_off
  n_old_j <- sum(index_old)

  classes_old_j <- fit$classes[index_old]


  #print(n_old_j)
  #print( index_old)
  #print(min_dist)

  if (n_old_j > 1) {
    coords_old_j <- coords_old[index_old, ]
    Y_old_j <- fit$Y_old[index_old, ]
  }

  if (n_old_j == 1) {
    coords_old_j <- matrix(coords_old[index_old, ], 1, 2) # coordinates are always 2 dimensional, x and y
    Y_old_j <- matrix(fit$Y_old[index_old, ], 1, L) # of dimension L
  }

  if (n_old_j == 0) {
    coords_old_j <- NULL
    Y_old_j <- NULL
  }


  return(tibble::lst(n_old_j, Y_new_j, Y_old_j, coords_new_j, coords_old_j, classes_old_j))
}


#***************************************************************************
#  Predict classification of a block NON CONDITIONAL
#***************************************************************************
# predict_in ----
#' Non conditional Prediction
#' @description Helper function used within predict_cond for NON conditional predictions.
#' @param fit LDA fit
#' @param y_new latent variable values for all cases belonging to the block to be predicted
#' @param trace Controls console output. NULL = no output, 1 = basic output,
#'  2 = detailed output.

predict_in <- function(fit, y_new, trace = FALSE, K) {

  hilf <- dim(y_new)
  n_new <- hilf[1]
  L <- hilf[2]
  #K <- L + 1

  y_n <- c(t(as.matrix(y_new)))

  ll <- rep(NA, K)

  # estimate loglikelihood for each class
  for (k in 1:K) {

    # SIGMA just refers to the observation belonging to the block
    lst <- lapply(rep(k, n_new), FUN = f_covs, fit$covs)
    SIGMA <- Matrix::bdiag(lst) # without spatial correlations
    mu <- rep(as.numeric(fit$means[k, ]), times = n_new)
    logdet_SIGMA <- NULL

    ll[k] <- loglik(y = y_n, mu = mu, C = SIGMA, logdet = logdet_SIGMA)
  } # end for(k in 1:K)

  # predict class with max loglikelihood
  hilf <- list(class = which.max(ll), class_pr = which.max(ll + log(fit$prior)),
               ll = max(ll),ll_pr = max(ll + log(fit$prior)))
  return(hilf)
} # end function

#***************************************************************************
#  Predict classification of a block CONDITIONAL
#***************************************************************************
# predict_in_cond ----
# just implement for equal.cov - Check for it??
#' Conditional Prediction
#' @description Helper function used within predict_cond for conditional predictions.
#' @param fit LDA fit providing covariance functions of training data
#' @param l component of transformed variable
#' @param Y_new_j latent variable values for all cases belonging to the block to be predicted
#' @param Y_old_j latent variable values for all training data cases within distance threshold
#' @param coords_new_j coords of all cases belonging to the block to be predicted
#' @param coords_old_j coord of all training data cases within distance threshold
#' @param classes_old_j class memberships of all training data cases within distance threshold
#' @param trace Controls console output. NULL = no output, 1 = basic output,
#'  2 = detailed output.
#' @param adjust_sigmas Should variances and/or covariances be adjusted? Default is 1 = adjust only variances, 2 = adjust both, NULL = no adjustments.
#' @export

calculate_LL <- function(fit, l ,Y_new_j, Y_old_j, coords_new_j, coords_old_j,
                            classes_old_j, trace) {

  n_new_j <- length(Y_new_j)
  K<-length(fit$lev) # number classes


  # create variables for included training observations if they are >0
  if (!is.null(Y_old_j)) {
    n_old_j <- length(Y_old_j)
    y_o <- Y_old_j
  } else {
    n_old_j <- 0
  }

  y_n <- Y_new_j
  coords_all_j <- rbind(coords_new_j, coords_old_j)
  n_all_j <- n_old_j + n_new_j


  #print(coords_all_j[duplicated(coords_all_j),]


    SIGMA <- geoR::varcov.spatial(
      coords = coords_all_j,
      cov.model = fit$variograms[[l]]$cov.model,
      kappa = fit$variograms[[l]]$kappa,
      nugget = fit$variograms[[l]]$nugget,
      cov.pars = fit$variograms[[l]]$cov.pars
    )$varcov
    # SIGMA <- Matrix::drop0(SIGMA,tol=1e-30)
    # make covariance Matrix sparse
    # loss of information but faster processing? For small dataset in 20 repetitions 3 secs difference in mean

  # retrieve covariance matrix for cases of block j (Matrix of size L*n_new_j x L*n_new_j)
  ind_new <- 1:n_new_j
  SIGMA_nn <- SIGMA[ind_new, ind_new]

  # estimate conditional multivariate normal distribution
  # (cases to be predicted conditional on nearby training observations)
  if (n_old_j > 0) {
    SIGMA_cond_lst <- try(calc_SIGMA_cond(n_new_j, n_all_j, SIGMA))

    if(!inherits(SIGMA_cond_lst, "try-error")){
    SIGMA_cond <- SIGMA_cond_lst$SIGMA_cond
    }else{
      # else use non-conditional covariance matrix
    SIGMA_cond <- SIGMA_nn
    }
  } else {
    # else use non-conditional covariance matrix
    SIGMA_cond <- SIGMA_nn
  }




  chol_SIGMA_cond <- try(chol(SIGMA_cond))

  #  print(head(as.matrix(SIGMA_cond))); print(head(as.matrix(chol_SIGMA_cond)))

  if (!inherits(chol_SIGMA_cond, "try-error")) {
    # cat("calculate chol Sigma-cond\n")
    logdet_SIGMA_cond <- 2 * sum(log(diag(chol_SIGMA_cond)))
  } else {
    warning("Covariance matrix not positive definite or other issues. \n")
  }

  ll <- ll1 <- rep(NA, K)

  for (k in 1:K) {

    mu_cond1 <- rep(as.numeric(fit$means[k,l]), times = n_new_j)
    logdet_SIGMA_cond1 <- NULL

    SIGMA_cond1 <- Matrix::Diagonal(n_new_j) # without spatial correlations

    # estimate loglikelihood for conditional distribution parameters
    ll1[k] <- loglik(y = y_n, mu = mu_cond1, C = SIGMA_cond1,
                     logdet = logdet_SIGMA_cond1)

    if (!inherits(chol_SIGMA_cond, "try-error")) {

      # here SIGMA refers to new cases (block j) and old cases (training data)
      mu_o <- c(t(as.matrix(fit$means[classes_old_j,l])))
      mu_n <- rep(as.numeric(fit$means[k,l]), times = n_new_j)

      # estimate conditional mean
      if (n_old_j > 0) {
        mu_cond <- as.numeric(mu_n + Matrix::t(SIGMA_cond_lst$SIGMA_oo_inv_SIGMA_on) %*% (y_o - mu_o))
      } else {
        mu_cond <- as.numeric(mu_n)
      }

      # estimate loglikelihood for conditional distribution parameters
      ll[k] <- loglik(y = y_n, mu = mu_cond, C = SIGMA_cond,
                      logdet = logdet_SIGMA_cond)

    }else{

      #if there was an issue with spatial model

      ll[k]<-ll1[k]


    }#end

  }#end   for (k in 1:K) {

    return(list(ll=ll,ll_ind=ll1))
  } # end function









