
#***************************************************************************
#  get_best_covmodel ----
#***************************************************************************
#' Best covariance model
#' @description Fits different given covariance function to the empirical semivariogram and chooses the best regarding sum of squares.
#' @param coords Dataframe with two rows containing coordinates
#' @param l The index of the used latent variables
#' @param rs List containing the residuals for each latent variable
#' @param omega0 Covariance matrizes for all latent variables
#' @param max_dist Maximum distance between all points
#' @param trace Controls console output. NULL = no output, 1 = basic output,
#'  2 = detailed output.
#' @param cov_models Covariance function models which are to be compared.


get_best_covmodel <- function(coords, l, rs, omega0,ini=NULL, max_dist, trace, cov_models,showplot=TRUE,phi_max=NULL) {

  # jitter duplicated coordinates
 #  coords <- geoR::jitter2d(coords, max=0.01)

  geodat <- geoR::as.geodata(
    obj = cbind(coords, rs[[l]]),
    coords.col = 1:2, data.col = 3
  )

  breaks<- c(seq(0,100,50),seq(200,max_dist,100))
  # empirical semivariogram
  # "modulus" is Hawkins and Cressie (see Cressie, 1993, pg 75)
  geodat_v1 <- geoR::variog(geodat,
                            # breaks = "default", max.dist = max_dist,
                            breaks = breaks, max.dist = max_dist,
                            option = "bin", estimator.type = "modulus",
                            messages = ifelse(trace > 1, TRUE, FALSE)
  )
  cat("variogram finished\n")

  SSs <- rep(Inf, length(cov_models))

  hilf <- list()

  # comparison of sums of square (SSs) for fitting different covariance models
  # to empirical semivariogram
  for (j in 1:length(cov_models)) {
  cat("Fit Model=",cov_models[j],"\n")
    # max.dist: by default the max.dist from object geodat_v1 is used, so no need
    # to set it again
    #geoExp <- try(geoR::variofit(geodat_v1,
    #                             ini = c(omega0[l], 100), fix.nugget = FALSE,
    #                             nugget = 0.1, fix.kappa = F, kappa = 0.5,
    #                             cov.model = cov_models[j], weights = "cressie",
     #                            messages = ifelse(trace > 1, TRUE, FALSE)
   if(is.null(ini)){ini<-c(0.2,100)}

    geoExp <- try(variofit1(geodat_v1,      #  sigmasq,    phi
                                 ini = c(omega0[l]-ini[2],ini[1]), fix.nugget = FALSE,
                                 nugget = ini[2], fix.kappa = F, kappa = 1,phi_max=phi_max,
                                 cov.model = cov_models[j], weights ="npairs", #cressie", #"cressie",
                                 messages = ifelse(trace > 1, TRUE, FALSE)

    ))
    #cat("Output\n")
    print(geoExp)
    #print(str(geoExp))

    if (!inherits(geoExp, "try-error")) {
      SSs[j] <- geoExp$value
      hilf[[j]] <- geoExp
    } else {
      SSs[j] <- Inf
    }
  } # for j

  k <- which.min(SSs)
  fin <- hilf[[k]]
  if(showplot){
  plot(geodat_v1, main = paste('geoR', fin$cov.model),ylim=c(0,1.2)); lines(fin, col='red')
  }

  return(fin)
}


########################################################################################

#.geoR.env <- new.env()

########################################################################################


    .geoR.env <- new.env()


variofit1 <-
      function (vario, ini.cov.pars, cov.model,
                fix.nugget = FALSE, nugget = 0,
                fix.kappa = TRUE, kappa = 0.5,
                simul.number = NULL,  max.dist = vario$max.dist, phi_max=NULL,
                weights, minimisation.function,
                limits = pars.limits(), messages, ...)
      {
        call.fc <- match.call()
        if(missing(messages))
          messages.screen <- as.logical(ifelse(is.null(getOption("geoR.messages")), TRUE, getOption("geoR.messages")))
        else messages.screen <- messages
        if(length(class(vario)) == 0 || all(class(vario) != "variogram"))
          warning("object vario should preferably be of the geoR's class \"variogram\"")
        if(!missing(ini.cov.pars)){
          if(any(class(ini.cov.pars) == "eyefit"))
            cov.model <- ini.cov.pars[[1]]$cov.model
          if(any(class(ini.cov.pars) == "variomodel"))
            cov.model <- ini.cov.pars$cov.model
        }
        if(missing(cov.model)) cov.model <- "matern"
        #cov.model <- match.arg(cov.model, choices =  .geoR.cov.models)
        if(cov.model == "stable") cov.model <- "powered.exponential"
        if(cov.model == "powered.exponential")
          if(limits$kappa["upper"] > 2) limits$kappa["upper"] <- 2
        ##  if(cov.model == "matern" | cov.model == "    powered.exponential" |
        ##     cov.model == "cauchy" | cov.model == "gneiting.matern")
        ##    fix.kappa <- TRUE
        if(missing(weights)){
          if(vario$output.type == "cloud") weights <- "equal"
          else weights <- "npairs"
        }
        else
          weights <- match.arg(weights, choices = c("npairs", "equal", "cressie"))
        if(messages.screen){
          cat(paste("variofit: covariance model used is", cov.model, "\n"))
          cat(paste("variofit: weights used:", weights, "\n"))
        }
        #  if(missing(minimisation.function)){
        #    if(weights == "equal") minimisation.function <- "nls"
        #    else minimisation.function <- "optim"
        #  }
        if(missing(minimisation.function))
          minimisation.function <- "optim"
        if(any(cov.model == c("linear", "power")) & minimisation.function == "nls"){
          cat("warning: minimisation function nls can not be used with given cov.model.\n          changing for \"optim\".\n")
          minimisation.function <- "optim"
        }
        if(minimisation.function == "nls" & weights != "equal"){
          warning("variofit: minimisation function nls can only be used with weights=\"equal\".\n          changing for \"optim\".\n")
          minimisation.function <- "optim"
        }
        if (is.matrix(vario$v) & is.null(simul.number))
          stop("object in vario$v is a matrix. This function works for only 1 empirical variogram at once\n")
        if (!is.null(simul.number))
          vario$v <- vario$v[, simul.number]
        ##
        ## Setting maximum distance
        ##
        if(mode(max.dist) != "numeric" || length(max.dist) > 1)
          stop("a single numerical value must be provided in the argument max.dist")
        if (max.dist == vario$max.dist)
          XY <- list(u = vario$u, v = vario$v, n=vario$n)
        else
          XY <- list(u = vario$u[vario$u <= max.dist],
                     v = vario$v[vario$u <= max.dist],
                     n = vario$n[vario$u <= max.dist])
        if(cov.model == "pure.nugget"){
          ##
          ## parameter estimation for model which does not require numerical minimisation
          ##
          minimisation.function <- "not used"
          message <- "correlation function does not require numerical minimisation"
          if(weights == "equal") lm.wei <- rep(1, length(XY$u))
          else lm.wei <- XY$n
          if(cov.model == "pure.nugget"){
            if(fix.nugget){
              temp <- lm((XY$v-nugget) ~ 1, weights = lm.wei)
              cov.pars <- c(temp$coef, 0)
            }
            else{
              temp <- lm(XY$v ~ 1, weights = lm.wei)
              nugget <- temp$coef
              cov.pars <- c(0,0)
            }
          }
          value <- sum((temp$residuals)^2)
        }
        else{
          if(messages.screen)
            cat(paste("variofit: minimisation function used:", minimisation.function, "\n"))
          ##
          ## setting things for numerical minimisation
          ##
          ##  Checking initial values
          ##
          umax <- max(vario$u)
          vmax <- max(vario$v)
          if(missing(ini.cov.pars)){
            ini.cov.pars <- as.matrix(expand.grid(c(vmax/2, 3*vmax/4, vmax),
                                                  seq(0, 0.8*umax, len=6)))
            if(!fix.nugget)
              nugget <- unique(c(nugget, vmax/10, vmax/4, vmax/2))
            if(!fix.kappa)
              kappa <- unique(c(kappa, 0.25, 0.5, 1, 1.5, 2))
            if(messages.screen)
              warning("initial values not provided - running the default search")
          }
          else{
            if(any(class(ini.cov.pars) == "eyefit")){
              init <- nugget <- kappa <- NULL
              for(i in 1:length(ini.cov.pars)){
                init <- drop(rbind(init, ini.cov.pars[[i]]$cov.pars))
                nugget <- c(nugget, ini.cov.pars[[i]]$nugget)
                if(cov.model == "gneiting.matern")
                  kappa <- drop(rbind(kappa, ini.cov.pars[[i]]$kappa))
                else
                  kappa <- c(kappa, ini.cov.pars[[i]]$kappa)
              }
              ini.cov.pars <- init
            }
            if(any(class(ini.cov.pars) == "variomodel")){
              nugget <- ini.cov.pars$nugget
              kappa <- ini.cov.pars$kappa
              ini.cov.pars <- ini.cov.pars$cov.pars
            }
          }
          if(is.matrix(ini.cov.pars) | is.data.frame(ini.cov.pars)){
            ini.cov.pars <- as.matrix(ini.cov.pars)
            if(nrow(ini.cov.pars) == 1)
              ini.cov.pars <- as.vector(ini.cov.pars)
            else{
              if(ncol(ini.cov.pars) != 2)
                stop("\nini.cov.pars must be a matrix or data.frame with 2 components: \ninitial values for sigmasq (partial sill) and phi (range parameter)\n")
            }
          }
          else
            if(length(ini.cov.pars) > 2)
              stop("\nini.cov.pars must provide initial values for sigmasq and phi\n")
          ##
          ## Preparing grid of initial values and choosing the best
          ##
          if(is.matrix(ini.cov.pars) | (length(nugget) > 1) | (length(kappa) > 1)) {
            if(messages.screen)
              cat("variofit: searching for best initial value ...")
            ini.temp <- matrix(ini.cov.pars, ncol=2)
            grid.ini <- as.matrix(expand.grid(sigmasq=unique(ini.temp[,1]),
                                              phi=unique(ini.temp[,2]),
                                              tausq=unique(nugget), kappa=unique(kappa)))
            ##  loss function:
            ###########################################################################
            v.loss <- function(parms, u, v, n, cov.model, weights){

              sigmasq <- parms[1]
              sigmasq <- max(min(sigmasq,1),0) #between 0 and 1

              phi <- parms[2]


              if(cov.model == "power") phi <- 2 * exp(phi)/(1+exp(phi))

              tausq <- 1-parms[1]

              # tausq<- min(tausq,0.6)
              # tausq<- max(tausq,0)

              kappa <- parms[3]

              if(cov.model == "power")
                v.mod <- tausq +
                cov.spatial(u, cov.pars=c(sigmasq, phi), cov.model="power", kappa=kappa)
              else
                #v.mod <- (sigmasq + tausq) -
                #cov.spatial(u, cov.pars=c(sigmasq, phi), cov.model = cov.model,
                #           kappa = kappa)
                v.mod <- 1 - cov.spatial(u, cov.pars=c(sigmasq, phi), cov.model = cov.model,
                            kappa = kappa)
              if(weights == "equal")
                loss <- sum((v - v.mod)^2)
              if (weights == "npairs")
                loss <- sum(n * (v - v.mod)^2)
              if (weights == "cressie")
                #loss <- sum((n/(v.mod^2)) * (v - v.mod)^2)
                loss <- sum((n/(v^2)) * (v - v.mod)^2)
              return(loss)
            }
            #####################################################################################


            grid.loss <- apply(grid.ini, 1, v.loss, u=XY$u, v=XY$v, n=XY$n, cov.model = cov.model, weights = weights)
            ini.temp <- grid.ini[which(grid.loss == min(grid.loss))[1],, drop=FALSE]
            if(is.R()) rownames(ini.temp) <- "initial.value"
            if(messages.screen){
              cat(" selected values:\n")
              print(rbind(round(ini.temp, digits=2), status=ifelse(c(FALSE, FALSE, fix.nugget, fix.kappa), "fix", "est")))
              cat(paste("loss value:", min(grid.loss), "\n"))
            }
            names(ini.temp) <- NULL
            ini.cov.pars <- ini.temp[1:2]
            nugget <- ini.temp[3]
            kappa <- ini.temp[4]
            grid.ini <- NULL
          }
          ##
          ## checking for unreasonable initial values
          ##
          if(ini.cov.pars[1] > 2*vmax)
            warning("unreasonable initial value for sigmasq (too high)")
          if(ini.cov.pars[1] + nugget > 3*vmax)
            warning("unreasonable initial value for sigmasq + nugget (too high)")
          if(vario$output.type != "cloud"){
            if(ini.cov.pars[1] + nugget < 0.3*vmax)
              warning("unreasonable initial value for sigmasq + nugget (too low)")
          }
          if(nugget > 2*vmax)
            warning("unreasonable initial value for nugget (too high)")
          if(ini.cov.pars[2] > 1.5*umax)
            warning("unreasonable initial value for phi (too high)")
          ##
          ## transforming kappa for constraint minimisation
          ##
          if(!fix.kappa){
            if(cov.model == "powered.exponential")
              Tkappa.ini <- log(kappa/(2-kappa))
            else
              Tkappa.ini <- log(kappa)
          }
          ##
          ## minimisation using "nls"
          ##
          if (minimisation.function == "nls") {
            if(ini.cov.pars[2] == 0) ini.cov.pars <- max(XY$u)/10
            if(kappa == 0) kappa <- 0.5
            if(cov.model == "power")
              Tphi.ini <- log(ini.cov.pars[2]/(2-ini.cov.pars[2]))
            else Tphi.ini <- log(ini.cov.pars[2])
            XY$cov.model <- cov.model
            ##
            if (fix.nugget) {
              XY$nugget <- as.vector(nugget)
              if(fix.kappa){
                XY$kappa <- as.vector(kappa)
                res <- nls((v-nugget) ~ matrix((1-cov.spatial(u,cov.pars=c(1,exp(Tphi)),
                                                              cov.model=cov.model, kappa=kappa)),
                                               ncol=1),
                           start=list(Tphi=Tphi.ini), data=XY, algorithm="plinear", ...)
              }
              else{
                if(cov.model == "powered.exponential")
                  res <- nls((v-nugget) ~ matrix((1-cov.spatial(u,cov.pars=c(1,exp(Tphi)),
                                                                cov.model=cov.model,
                                                                kappa=(2*exp(Tkappa)/(1+exp(Tkappa))))),
                                                 ncol=1),
                             start=list(Tphi=Tphi.ini, Tkappa = Tkappa.ini),
                             data=XY, algorithm="plinear", ...)
                else
                  res <- nls((v-nugget) ~ matrix((1-cov.spatial(u,cov.pars=c(1,exp(Tphi)),
                                                                cov.model=cov.model,
                                                                kappa=exp(Tkappa))), ncol=1),
                             start=list(Tphi=Tphi.ini, Tkappa = Tkappa.ini),
                             data=XY, algorithm="plinear", ...)
                kappa <- exp(coef(res)["Tkappa"])
                names(kappa) <- NULL
              }
              cov.pars <- coef(res)[c(".lin", "Tphi")]
              names(cov.pars) <- NULL
            }
            else{
              if(fix.kappa){
                XY$kappa <- kappa
                res <- nls(v ~ cbind(1,(1- cov.spatial(u, cov.pars=c(1,exp(Tphi)),
                                                       cov.model = cov.model, kappa=kappa))),
                           start=list(Tphi=Tphi.ini), algorithm="plinear", data=XY, ...)
              }
              else{
                if(cov.model == "powered.exponential")
                  res <- nls(v ~ cbind(1, (1-cov.spatial(u, cov.pars=c(1, exp(Tphi)),
                                                         cov.model = cov.model,
                                                         kappa=(2*exp(Tkappa)/(1+exp(Tkappa)))))),
                             start=list(Tphi=Tphi.ini, Tkappa = Tkappa.ini),
                             algorithm="plinear", data=XY, ...)
                else
                  res <- nls(v ~ cbind(1, (1-cov.spatial(u, cov.pars=c(1, exp(Tphi)),
                                                         cov.model = cov.model,
                                                         kappa=exp(Tkappa)))),
                             start=list(Tphi=Tphi.ini, Tkappa = Tkappa.ini),
                             algorithm="plinear", data=XY, ...)
                kappa <- exp(coef(res)["Tkappa"]);names(kappa) <- NULL
              }
              nugget <- coef(res)[".lin1"];names(nugget) <- NULL
              cov.pars <- coef(res)[c(".lin2", "Tphi")]
              names(cov.pars) <- NULL
            }
            if(cov.model == "power")
              cov.pars[2] <- 2 * exp(cov.pars[2])/(1+exp(cov.pars[2]))
            else cov.pars[2] <- exp(cov.pars[2])
            if(nugget < 0 | cov.pars[1] < 0){
              warning("\nvariofit: negative variance parameter found using the default option \"nls\".\n        Try another minimisation function and/or fix some of the parameters.\n")
              temp <- c(sigmasq=cov.pars[1], phi=cov.pars[2], tausq=nugget, kappa=kappa)
              print(rbind(round(temp, digits=4),
                          status=ifelse(c(FALSE, FALSE, fix.nugget, fix.kappa), "fix", "est")))
              return(invisible())
            }
            value <- sum(resid(res)^2)
            message <- "nls does not provides convergence message"
          }
          ##
          ## minimisation using "optim" or "nlm"
          ##
          if (minimisation.function == "nlm" | minimisation.function == "optim") {
            ##
            ## Preparing lists for the minimiser
            ##
            .global.list <- list(u = XY$u, v = XY$v, n=XY$n, fix.nugget = fix.nugget,
                                 nugget = nugget, fix.kappa = fix.kappa, kappa = kappa,
                                 cov.model = cov.model, m.f = minimisation.function,
                                 weights = weights)
            ##
            ## Preparing initial value
            ##
            ini <- ini.cov.pars
            if(cov.model == "power") ini[2] <- log(ini[2]/(2-ini[2]))
            if(cov.model == "linear") ini <- ini[1]
            if(fix.nugget == FALSE) ini <- c(ini, nugget)
            ## setting kappa > 0 for both methods
            if(!fix.kappa) ini <- c(ini, Tkappa.ini)
            names(ini) <- NULL
            if(minimisation.function == "nlm"){
              result <- nlm(.loss.vario, ini, g.l = .global.list, ...)
              result$par <- result$estimate
              result$value <- result$minimum
              result$convergence <- result$code
              if(!is.null(get(".temp.theta", pos=.geoR.env)))
                result$par <- get(".temp.theta", pos=.geoR.env)
            }
            else{
              #        if(fix.kappa == FALSE) ini <- c(ini, kappa)
              #        names(ini) <- NULL
              lower.l <- sapply(limits, function(x) x[1])
              upper.l <- sapply(limits, function(x) x[2])
              if(fix.kappa == FALSE){
                if(fix.nugget){
                  lower <- lower.l[c("sigmasq.lower", "phi.lower","kappa.lower")]
                  upper <- upper.l[c("sigmasq.upper", "phi.upper","kappa.upper")]
                }
                else{
                  lower <- lower.l[c("sigmasq.lower", "phi.lower",
                                     "tausq.rel.lower", "kappa.lower")]
                  upper <- upper.l[c("sigmasq.upper", "phi.upper",
                                     "tausq.rel.upper", "kappa.upper")]
                }
              }
              else{
                if(cov.model == "power"){
                  if(fix.nugget){
                    lower <- lower.l[c("sigmasq.lower", "phi.lower")]
                    upper <- upper.l[c("sigmasq.upper", "phi.upper")]
                  }
                  else{
                    lower <- lower.l[c("sigmasq.lower", "phi.lower", "tausq.rel.lower")]
                    upper <- upper.l[c("sigmasq.upper", "phi.upper", "tausq.rel.upper")]
                  }
                }
                else{
                  lower <- lower.l["phi.lower"]
                  upper <- upper.l["phi.upper"]
                }
              }
              #cat("ini\n")

              ini<- ini[c(1,2,4)]
              # kappa
              ini[3]<- log(1)  # starts at 1
              #print(ini)
              #sigmsq
              upper[1]<-1  # sigmasq bounded by 0.1 and 1
              lower[1]<-0.5
              #
              if(is.null(phi_max)){upper[2]<-vario$max.dist}else{upper[2]<-phi_max}

              lower[2]<- 20

              #m log(kappa)
              upper[4]<-  4
              lower[4]<- -4

              lower1<-lower[-3]  # remove tausq
              upper1<-upper[-3]  # remove tausq

              if(0){
              cat("ini,lower,upper\n")
              print(ini)
              print(lower)
              print(upper)
              print(lower1)
              print(upper1)
              }
              # print(.global.list)


              cat("############### BEGIN SANN ################\n")
              cat("Cov-model",.global.list$cov.model,"\n")
              result<- GenSA(par=ini, fn=.loss.vario, lower=lower1, upper=upper1,
                             control=list(temperature=100,max.time=2000,maxit=1e4),g.l = .global.list)

              cat("final value SANN",result$value,"\n")
              #if(is.na(result$par[4])){result$par[4]<-1}
              if(0){
               cat("first SANN\n")
               print(result)
               print(result$par)
              }

              cat("############### END SANN - BEGIN OPTIM ################\n")

              # get more precise estimates with local optimiser
              result <- optim(result$par, .loss.vario, method = "L-BFGS-B",
                              hessian = TRUE, lower = lower1,
                              upper = upper1, g.l = .global.list,control=list(trace=2,maxit=200))


              cat("############### END  OPTIM ################\n")
              #cat("2nd L-BFGS-B\n")
              #print(result)
              #print(result$par)

              #        require(methods)
              #        if(exists("trySilent"))
              #          hess <- trySilent(solve(as.matrix(result$hessian)))
              #        else{
              #          op.sem <- options()$show.error.messages
              #          options(show.error.messages = FALSE)
              #          hess <- try(solve(as.matrix(result$hessian)))
              #          options(show.error.messages = op.sem)
              #        }
              #        if(!inherits(hess, "try-error")) hess <- sqrt(diag(hess))
              #        else print("WARNING: unable to compute the hessian")
            }
            value <- result$value
            if(0){
            cat("optim-par - result\n")
            print(result$par)

            }
            hilf<-result$par
            result$par<-c(hilf[1:2],1-hilf[1],hilf[3])
            #print(result$par)


            message <- paste(minimisation.function, "convergence code:", result$convergence)
            if(cov.model == "linear")
              result$par <- c(result$par[1],1,result$par[-1])
            cov.pars <- as.vector(result$par[1:2])
            if(cov.model == "power")
              cov.pars[2] <- 2 * exp(cov.pars[2])/(1+exp(cov.pars[2]))
            if(!fix.kappa){
              if (fix.nugget)
                kappa <- result$par[3]
              else{
                nugget <- result$par[3]
                kappa <- result$par[4]
              }
              ## kappa now > 0 for both nlm() and optim()
              ##        if(minimisation.function == "nlm"){
              if(.global.list$cov.model == "powered.exponential")
                kappa <- 2*(exp(kappa))/(1+exp(kappa))
              else kappa <- exp(kappa)
              ##        }
            }
            else
              if(!fix.nugget)
                nugget <- result$par[3]
          }
        }
        ##
        ## Estimating implicity beta
        ##

        ##
        ## Preparing output
        ##
        cat("prac-range start\n")
        estimation <- list(nugget = nugget, cov.pars = cov.pars,
                           cov.model = cov.model, kappa = kappa, value = value,
                           trend = vario$trend, beta.ols = vario$beta.ols,practicalRange =  NULL,
                          practicalRange = try(practicalRange(cov.model=cov.model,
                                                           phi = cov.pars[2], kappa = kappa)),
                           max.dist = max.dist,
                           minimisation.function = minimisation.function)

        cat("prac-range end\n")
        #  if(exists("hess")) estimation$hessian <- hess
        estimation$weights <- weights
        if(weights == "equal") estimation$method <- "OLS"
        else estimation$method <- "WLS"
        estimation$fix.nugget <- fix.nugget
        estimation$fix.kappa <- fix.kappa
        estimation$lambda <- vario$lambda
        estimation$message <- message
        estimation$call <- call.fc
        oldClass(estimation) <- c("variomodel", "variofit")
        return(estimation)
      }


######################################################################################
    ".loss.vario" <-
      function (theta,g.l)
      {

       #print(theta)
       #print(g.l)
       #break
       #stop


       penalty <- 0
        ##
        ## reading parameters
        ##
        if(!g.l$fix.kappa){
        # if not fixed
            kappa <- exp(theta[3])
            # tau is fixed

        }else{
          # if fixed, take kappa from global list
          kappa <- g.l$kappa
        }


        tausq <- 1- theta[1]

        ##
        sigmasq <- theta[1]
        phi <- theta[2]


        sill.total <- sigmasq + tausq
        ##
        ## Computing values for the theoretical variogram
        ##

          gammaU <- sill.total - cov.spatial(g.l$u, cov.model = g.l$cov.model,
                                             kappa = kappa, cov.pars = c(sigmasq, phi))
        ##
        ## Computing loss function
        ##
        if(g.l$weight == "equal")
          loss <- sum((g.l$v - gammaU)^2)
        if (g.l$weights == "npairs")
          loss <- sum(g.l$n * (g.l$v - gammaU)^2)
        if (g.l$weights == "cressie")
          # loss <- sum((g.l$n/(gammaU^2)) * (g.l$v - gammaU)^2)
          loss <- sum((g.l$n/(g.l$v^2)) * (g.l$v - gammaU)^2)
        if(loss > (.Machine$double.xmax^0.5) | loss == Inf | loss == -Inf | is.nan(loss))
          loss <- .Machine$double.xmax^0.5
        return(loss + penalty)
      }#end  .loss.vario

