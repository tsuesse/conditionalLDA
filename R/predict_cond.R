#***************************************************************************
# predict_cond----
#***************************************************************************
#' Conditional LDA Prediction
#' @description  This function predicts class membership for an entire block (or a single observation).
#' The conditional prediction approach includes nearby training observations.
#' Therefore coords must be given and LDA fit must provide covariance functions.
#' @param object LDA fit (object of class lda?) Has to be provide covariance functions and training data information for conditional prediction.
#' @param newdata  data frame of cases to be classified.
#' @param coords vector of length 2 defining the variables in newdata that contain
#' x and y coordinates of sample locations. The default is c("x", "y"). Must be
#' given to perform conditional prediction.
#' @param adjust_sigmas Should variances and/or covariances be adjusted? Default is 1 = adjust only variances, 2 = adjust both, NULL = no adjustments.
#' @param trace Controls console output: 0 = no output, 1 = basic output , 2 = detailed output.
#' @param LL Should loglikelihood for training data be estimated?
#' @param blocks character string defining variable of newdata containing block identifiers.
#' @param prior If TRUE (default) predictions including prior class probabilities are returned.
#' @return Predicted class memberships of given data.
#' @export


predict_cond <- function(object, newdata, coords = c("x", "y"),
                         trace = 1, LL = TRUE, blocks = "field", max_dist=NULL, corr_min= 0.01) {

  fit <- object
  if(is.null(max_dist)){ max_dist<- fit$max_dist}

  # fit<-fit4;blocks = block;newdata<-data_test;corr_min= 0.01

  # dimension reduction of original values for higher computational efficiency
  Y_new <- as.matrix(newdata[, fit$predictors]) %*% fit$scaling

  response <- as.character(fit$terms[[2]])

  hilf <- dim(Y_new)
  n_new <- hilf[1]
  L <- hilf[2] # number of latent variables
  K <- L + 1 # number of classes





  # check if cases are grouped in blocks and define unit to be predicted
  if (!is.null(blocks)) {
    block <- as.integer(newdata[, blocks])
  } else {
    block <- as.integer(row.names(newdata)) # ??? or should ID column be delivered? Is that it or are more adjustments necessary?
  }

  unique_block<-sort(unique(block))
  n_cond<- matrix(0,length(unique_block),L)  # number of pixels of training set within distance of field


  ll_un <- ll_un_pr <-  ll<- ll_pr<- pred_pr <- pred  <- pred_un <- pred_un_pr<- rep(NA, n_new)
  # corr_mean<-matrix(0,n_new,9)
  # corr_sum<- rep(0, n_new)
  # corr_sum_block<- rep(0,length(unique_block))

  #*************************************************
  # Conditional approach
  #*************************************************
  if (!is.null(coords)) {

    if (is.null(fit$variograms) || is.null(fit$coords)){
      warning("For conditional predictions fit must provide covariance
              functions and training coords." + "\n" +
                "Will proceed with non conditional predictions...")
      coords <- NULL
    }




    # get coords of test (old) and training (new) set
    coords_new <- as.matrix(newdata[, coords])
    coords_old <- as.matrix(fit$coords) # could also be derived within get_nearby_obs,
    # because no more need afterwards, but then it has to be derived j times.
    # This way only has to be given j times to function...


    cat("Spatial prediction starts\n")

    cut_off<-rep(0,L)

    for(l in 1:L){

      if(corr_min==0.01){
    # if default just use the already practical range computed

     cut_off[l]<-fit$variograms[[l]]$practicalRange


    }else{



      res<-fit$variograms[[l]]
      hilfPR<- try(practicalRange(cov.model=res$cov.model,phi = res$cov.pars[2],
                           kappa = res$kappa,correlation=corr_min))
      if(inherits(hilf,"try-error")){

        cut_off[l]<-res$cov.pars[2] # use range parameter instead
      }else{
        cut_off[l]<-hilfPR
      }# end if error


    }#end if(corr_min==0.01){
    }#end for(l in 1:L){

    cat("For corr=",corr_min," the practical ranges are",round(cut_off,3))
    cat(" (max_dist",max_dist,")\n")

    cut_off<-pmin(cut_off,max_dist)


    if (trace > 0) {
      #cat("cut_off",cut_off,"\n")
      cat("Conditional prediction starts...\n", "\n")
    }


    #print(cut_off)
    #cut_off<-300
    #cut_off<-cut_off*1.5

    j_field<-0

    #print(unique_blocks)



    st_pred <- system.time({
      for (j in unique_block) {
        index <- block == j
        n_new_j <- sum(index)

        j_field<-j_field+1

      ll<-ll1<-matrix(0,L,K)
      for( l in 1:L){
        # get training observations within distance threshold to block
        j_old <- get_nearby_obs(fit, L, index, n_new_j, Y_new, coords_old,
                                coords_new, cut_off[l])


        #cat("Start Prediction\n")

        if (trace > 1) {
          cat("Block", j, "of size", n_new_j, "with", j_old$n_old_j,
              #"points within the distance of", round(cut_off, 1),"sum of corr",corr_sum[index][1], "\n")
              "points within the distance of", round(cut_off[l], 1), "for component ",l," \n")
        }

        n_cond[j_field,l]<-j_old$n_old_j


        # calculate for each class k of K classes the log-likelihood for each dimension l
        hilf <- calculate_LL(
          fit = fit, l=l, Y_new_j = j_old$Y_new_j[,l], Y_old_j = j_old$Y_old_j[,l],
          coords_new_j = j_old$coords_new_j,
          coords_old_j = j_old$coords_old_j,
          classes_old_j = j_old$classes_old_j, trace = trace
        )
        ll[l,]<-hilf$ll
        ll1[l,]<-hilf$ll_ind

        }#end for( l in 1:L){

      ll     <- apply(ll,2,sum)
      ll1    <- apply(ll1,2,sum)

      hilf <- list(class = which.max(ll), class_pr = which.max(ll + log(fit$prior)),
                   # dist_corr=dist_corr,corr_sum=corr_sum,
                   ll = max(ll),ll_pr = max(ll + log(fit$prior)),
                   class_un = which.max(ll1), class_un_pr = which.max(ll1 + log(fit$prior)),
                   ll_un = max(ll1),ll_un_pr = max(ll1 + log(fit$prior))
      )





        pred[index] <- hilf$class
        pred_pr[index] <- hilf$class_pr
        pred_un[index] <- hilf$class_un
        pred_un_pr[index] <- hilf$class_un_pr
        ll[index]  <- hilf$ll
        ll_pr[index]<- hilf$ll_pr
        ll_un[index]  <- hilf$ll_un
        ll_un_pr[index]<- hilf$ll_un_pr

        #corr_mean[index,]<- apply(hilf$dist_corr,2,mean)
        #corr_sum[index]<-hilf$corr_sum
        #corr_sum_block[j_field]<-hilf$corr_sum



        # retrieve true classes of testdata for output
        cases <- newdata[index, ]

        if (trace > 1) {
          cat("True type", cases[ , response][1],  # has to be removed for display in RSAGA
              "Predicted type", hilf$class,
              "(no prior) /", hilf$class_pr, "(prior)", "\n")
        }
      } # end for (j in unique(block))
    })

    if (trace > 0) {
      cat("Time needed in secs for prediction", st_pred[3], "\n")
    }

    if (LL) {
      #*******************************************
      # calculate Log-likelihood for Training set
      #*******************************************
      if (trace > 0) {
        cat("Calculation Log-likelihood of training set starts... \n")
      }

      st_LL <- system.time({
        # loglikelihood calculation with spatial SIGMA (optionally adjusted sigmas)
        LLmodel <- loglik_train_sp(fit, L)
      })

      if (trace > 0) {
        cat("Time needed in secs for Log-likelihood of training set",
            st_LL[3], "\n")
      }

      if (trace > 1) {
        cat("Log-likelihood of training set:", LLmodel, "\n")
      }

    } else {
      LLmodel <- NULL
    }

    # END conditional approach

    #*************************************************
    # NON - Conditional approach
    #*************************************************
  } else { # coords = NULL

    st_pred_ind<-system.time({

    if (trace > 0) {
      cat("Non - Conditional spatial (independent) prediction starts...\n", "\n")
    }

    for (j in unique_block) {
      index <- block == j

      if (sum(index) > 1) {
        y_new_j <- Y_new[index, ]
      } else {
        y_new_j <- matrix(Y_new[index, ], 1) #, dim(newdata)[2])
      }

      # non conditional prediction only based on testdata
      hilf <- predict_in(fit = fit, y_new = y_new_j, trace = trace, K = K)

      #print(hilf$ll)

      pred[index] <- hilf$class
      pred_pr[index] <- hilf$class_pr
      ll[index]  <- hilf$ll
      ll_pr[index]<- hilf$ll_pr




      # retrieve true classes of testdata for output
      cases <- newdata[index, ]

      if (trace > 1) {
        cat("Block", j, "of size", sum(index), "True type",
            cases[ , response][1], "Predicted type", hilf$class, "(no prior) /",
            hilf$class_pr, "(prior) \n")
      }
    } # end for (j in unique(block))

    if (LL) {
      # compute Loglikelihood of model for training data set
      lst <- lapply(fit$classes, FUN = f_covs, fit$covs)
      SIGMA_oo <- Matrix::bdiag(lst) # without spatial correlations
      mu_o <- c(t(fit$means[fit$classes, ]))
      y_o <- c(t(as.matrix(fit$Y_old)))

      LLmodel <- try(loglik(y_o, mu_o, SIGMA_oo))
      if (trace > 1) {
        cat("Log-likelihood of training set:", LLmodel, "\n")
      }

    } else {
      LLmodel <- NULL
    } # end if LL


    }) #end system.time

    st_pred_ind

  } # end non conditional



  # retrieve actual classnames out of fit
  labs <- fit$lev


  # return both prior and non prior predictions for prior = T
  pred_pr <- factor(pred_pr, levels = 1:K, labels = labs)
  pred <- factor(pred, levels = 1:K, labels = labs)

  pred_un_pr <- factor(pred_un_pr, levels = 1:K, labels = labs)
  pred_un <- factor(pred_un, levels = 1:K, labels = labs)


  #if (prior) {
    if (LL) {
      return(tibble::lst(pred_pr,pred_un_pr,pred,pred_un,LLmodel,ll_pr,ll_un_pr,n_cond))
    } else {
      return(tibble::lst(pred_pr,pred_un_pr,pred,pred_un,ll_pr,ll_un_pr,n_cond))
     }
  #}


} #**************************************************************************
