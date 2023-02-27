#***************************************************************************
# lda_sp----
#***************************************************************************
# Function to fit a spatial LDA model that (optionally) includes a covariance
# function for each latent variable.
# (could include an option without formula but separate dataframes with class
# membership and predictors; also for coords)
#' Spatial LDA
#' @description Fits a spatial LDA model that includes covariance
#' functions for each latent variable.
#' @param formula of the form \emph{class} ~ x1 + x2 + ... that describes the response,
#' i.e. the class membership and the predictors on the right hand side.
#' @param data data frame containing the variables in the formula.
#' @param covs_equal Do the classes have equal covariance matrices? (LDA assumption)
#' OR Should class covariances be treated as equal?
#' @param coords vector of length 2 defining the variables in data that contain
#' x and y coordinates of sample locations. The default is c("x", "y"). Must be
#' given to perform spatial LDA.
#' @param trace Controls console output. NULL = no output, 1 = basic output,
#'  2 = detailed output.
#' @param cov_models Character vector defining covariance function models which
#' are to be compared for each latent variable. Default is c("matern", "exponential",
#' "gaussian", "spherical"). Any combination of those as well as "linear",
#' "circular", "cubic", "wave" and ??"power"?? is possible. If set to NULL non
#' spatial/regular LDA will be performed.
#' @export

lda_sp <- function(formula, data, covs_equal = TRUE,  coords = c("x", "y"), trace = 1,
                   cov_models = c("matern", "exponential", "gaussian", "spherical"),max_dist=NULL,phi_max=NULL,
                   showplot=TRUE,ini=NULL,fitLDA =NULL,fixed.model=NULL)
 {


  #cov_models<-c("exponential", "spherical")

  response <- as.character(formula.tools::lhs(formula))
  table_class <- table(data[, response])
  classes <- data[, response]

  if (min(table_class) == 0) {
    stop("One class has sample size ZERO\n")
  }

  if(is.null(fitLDA)){
  # regular LDA fit
  fit <- MASS::lda(formula, data)
  #fit <- mda::fda(formula, data)
  }else{
  fit<-fitLDA
  }
  fit$predictors <- attr(fit$terms, "term.labels")
  # dimension reduction with linear coefficients for higher efficiency
  Y_old <- as.matrix(data[, fit$predictors]) %*% fit$scaling


  hilf <- dim(Y_old)
  L <- hilf[2] # L = K-1
  K <- L + 1 # number of classes

  # get class means of dimension reduced observations
  means <- aggregate(Y_old ~ classes, FUN = mean)
  means <- means[, -1] # remove identifier of cluster, not needed

  # attach training data information to fit for use in conditional prediction
  fit$classes <- classes # class membership
  fit$Y_old <- Y_old # dimension reduced/latent data
  fit$data <- data # original data
  fit$means <- means # class means

  # optional covariance adjustment to meet LDA condition
  if (covs_equal) {
    for (k in 1:K) {
      fit$covs[[k]] <- diag(L)
    }
  } else {
    covs <- by(Y_old, classes, cov)
    fit$covs <- covs
  }

  # check if covariance models are given to proceed with spatial modelling
  if (!is.null(cov_models)) {
    sp_corr <- TRUE
  } else {
    sp_corr <- FALSE
    if (trace > 0) {
      cat("No covariance models are given. Non spatial fit will be returned.\n")
    }
  }

  #*******************************************
  # START SPATIAL CORRELATION
  if (sp_corr) {

    if (is.null(coords)) {
      stop("For spatial fit coords must be given. If you want a non-spatial fit set cov_models = NULL.")
    }

    if (trace > 0) {
      cat("Spatial correlation modelling starts...\n", "\n")
    }

    coords <- data[, coords]

    # get distances among training observations
    if(is.null(max_dist)){
    distances <- dist(coords)
    max_dist <- max(distances)

    if (trace > 1) {
      cat("Max distance between observations", max_dist, "\n", "\n")
    }

    }else{

      if (trace > 1) {
        cat("Max distance used", max_dist, "\n", "\n")
      }


    }#end



    # retrieve residuals and covariances for each predictor l (latent variables!)
    rs <- list()
    omega0 <- rep(0, L)

    for (l in 1:L) {
      rs[[l]] <- (Y_old[, l] - means[classes, l])
      omega0[l] <- fit$covs[[1]][l, l]
    }





    # comparison of covariance models and selection of best model for each predictor
    results <- list()
    for (l in 1:L) {

       cat("pars=",ini," (phi,tausq)\n")

    # fit models if nothing is given
    if(is.null(fixed.model)){
      results[[l]] <- get_best_covmodel(coords, l, rs, omega0,ini, max_dist, trace,
                                        cov_models,showplot=showplot,phi_max=phi_max)




      # Plotting only possible if variogram is returned by get_best_covmodel
      # Is plot nice to have or doesn't matter?

    # if (trace > 1) {
    #   plot(results[[l]], main = paste('geoR', results[[l]]$cov.model)); lines(results[[l]], col='red')
     # }
     hilf<-
       try(practicalRange(cov.model=results[[l]]$cov.model,phi = results[[l]]$cov.pars[2],
                      kappa = results[[l]]$kappa,correlation=1e-2))

    }else{
      results[[l]]<-list()
      # fixed.model
      #cat("hello0\n")
      results[[l]]$cov.model <- fixed.model$cov.model
      results[[l]]$cov.pars  <- fixed.model$cov.pars #  sigma^2,phi
      results[[l]]$kappa     <- fixed.model$kappa    #kappa
      results[[l]]$nugget     <- fixed.model$nugget   # nugget
      hilf<-fixed.model$practicalRange
      #cat("hello1\n")
    }


     if(inherits(hilf,"try-error")){
       results[[l]]$practicalRange <- max_dist
     }else{
       results[[l]]$practicalRange <- hilf
     }


     if (trace > 0) {
       pars<-c(results[[l]]$nugget,results[[l]]$cov.pars[2],results[[l]]$cov.pars[1],results[[l]]$kappa)
       names(pars)<-c("tausq","sigmasq","phi","kappa")
       cat("**************************************************************************\n")
       cat("Model chosen for Variable l =", l, " is ", results[[l]]$cov.model, "\n")
       cat("with Parameters ", pars, " (tausq,phi,sigmasq,kappa)\n")
       cat("Practical Range with Corr<0.01 is",results[[l]]$practicalRange,"\n")
       cat("**************************************************************************\n")
     }




    } # end for l

    fit$coords <- coords
    fit$variograms <- results




  } # end if(sp_corr)
  # END SPATIAL CORRELATION
  #*******************************************

  fit$max_dist<-max_dist

  return(fit)
} # end lda_sp function

