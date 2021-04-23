#' @title Fit linear ERR model
#' @description Fits the linear ERR model on a dataset
#' @param data data frame containing matched case-control data, with a number of columns for doses to different locations, a column containing matched set numbers, a column containing the case's tumor location (value between 1 and the number of locations, with location \eqn{x} corresponding to the \eqn{x}-th column index in \code{doses}) and a column serving as a case-control indicator. Other covariates can also be included, in this case a parameter for each covariate column will be estimated. Hence factor variables need to be converted to dummy variables using \code{model.matrix}. If using \code{ccmethod='meandose'}, a column for tumor location is still required but in this case the column can be a vector of ones.
#' @param doses vector containing the indices of columns containing dose information.
#' @param set column index containing matched set numbers.
#' @param status column index containing case status.
#' @param loc column index containing the location of the matched set's case's second tumor.
#' @param corrvars vector containing the indices of columns containing variables to be corrected for.
#' @param repar reparametrize to \eqn{\beta=exp(\xi)}? Defaults to \code{FALSE}
#' @param ccmethod choice of method of analysis: one of meandose, CCML, CCAL or CL. Defaults to CCAL
#' @param initpars initial values for parameters, default is 0 for all parameters. If supplying a different vector, use a vector with an initial value for \eqn{\beta} or \eqn{\xi}, one for all of the other location effects and one for each other covariate (in that order). Note that if \code{repar=TRUE}, the initial value is used for \eqn{\xi}.
#' @param fitopt list with options to pass to \code{control} argument of optimizer
#' @param fitNull boolean: also fit model without dose effect? Defaults to \code{TRUE}. Note: the same optimization algorithm that was used for the MLE will be used for the null model, even if the null model only has one parameter (see details)
#' @param useOld if TRUE, a previous (slower) implementation of the log-likelihood function will be used. Defaults to \code{FALSE}
#' @param uplimBeta upper limit for \eqn{\beta=exp(\xi)}, default value 5. This is used for constraining the MLE estimation in some settings and for the jackknife inclusion criteria, and can be infinite except when Brent optimization is used (see help for \code{linearERR})
#' @return Object with components:
#' \item{fit}{object produced by \code{mle2}}
#' \item{nullfit}{fit without dose effect produced by \code{mle2}}
#' \item{proflik}{profile likelihood: one-dimensional function of \eqn{\beta} or \eqn{\xi}. Note that the optimization used is the same as for the MLE, leading to one-dimensional Nelder-Mead optimization in certain cases (see details of \code{linearERR})}
#' @importFrom numDeriv hessian
#' @importFrom stats constrOptim
#' @details This is a stripped down version of \code{linearERR}, and should only be used when that function does not suffice. For more details refer to the help of \code{linearERR}.
#' @seealso linearERR
#' @examples
#' data(linearERRdata1)
#'
#' fitmeandose <- linearERRfit(data=linearERRdata1, set=1, doses=2:6,
#' status=8, loc=7, corrvars=9, repar=FALSE, ccmethod="meandose")
#'
#' fitCCML <- linearERRfit(data=linearERRdata1, set=1, doses=2:6,
#' status=8, loc=7, corrvars=9, repar=FALSE, ccmethod="CCML")
#'
#' fitCCAL <- linearERRfit(data=linearERRdata1, set=1, doses=2:6,
#' status=8, loc=7, corrvars=9, repar=FALSE, ccmethod="CCAL")
#'
#' fitCL <- linearERRfit(data=linearERRdata1, set=1, doses=2:6,
#' status=8, loc=7, corrvars=NULL, repar=FALSE, ccmethod="CL")
#'
#' fitmeandose$fit$par
#' fitCCML$fit$par
#' fitCCAL$fit$par
#' fitCL$fit$par
#' @export


linearERRfit <- function(data, doses, set, status, loc, corrvars=NULL, repar=FALSE, ccmethod="CCAL", initpars=rep(0,length(doses)+length(corrvars)), fitopt=list(maxit=5000), fitNull=TRUE, useOld=FALSE, uplimBeta=5){

  if(ccmethod=="CL") corrvars <- NULL

  if(ccmethod %in% c("CCML", "meandose") & length(corrvars)==0){
    opt_method <- "Brent"
    if(is.infinite(uplimBeta)) stop("Please provide a finite value for uplimBeta to use in Brent optimization")
  } else{
    opt_method <- "constrOptim"
    if(repar & is.infinite(uplimBeta)) stop("uplimBeta needs to be finite for the current settings")
  }

  if(is.null(fitopt$maxit)){
    fitopt <- c(fitopt, list(maxit=5000))
  }
  if(is.null(fitopt$reltol)){
    fitopt <- c(fitopt, list(reltol=1e-10))
  }
  # if (is.null(fitopt$pgtol) & opt_method=="L-BFGS-B"){
  #   fitopt <- c(fitopt, list(pgtol=1e-8))
  # }
  # if (is.null(fitopt$factr) & opt_method=="L-BFGS-B"){
  #   fitopt <- c(fitopt, list(factr=1e4))
  # }
  # if (is.null(fitopt$ndeps) & opt_method=="L-BFGS-B"){
  #   parlen <- 1+ifelse(ccmethod %in% c("CCAL","CL"),length(doses)-1,0)+length(corrvars)
  #   fitopt <- c(fitopt, list(ndeps=rep(1e-5, parlen)))
  # }

  likfun <- function(params){
    if(repar) params[1] <- exp(params[1])
    if(useOld){
      linERRloglikold(params, data=data, set=set, doses=doses,status=status,loc=loc,corrvars=corrvars, ccmethod=ccmethod)
    } else {
      linERRloglik(params, data=data, set=set, doses=doses,status=status,loc=loc,corrvars=corrvars, ccmethod=ccmethod)
    }
  }
  scorefun <- function(params){
    linERRscore(params, data=data, set=set, doses=doses,status=status,loc=loc,corrvars=corrvars, ccmethod=ccmethod, repar=repar)
  }

  if(opt_method!="Brent"){
    names <- c(ifelse(repar,"xi","beta"), paste0("alpha",2:length(doses)),names(data)[corrvars])
    names(initpars) <- names
    #parnames(likfun) <- names
  } else {
    names <- ifelse(repar, "xi","beta")
    initpars <- list(params=as.numeric(initpars[1]))
    names(initpars) <- names
    #parnames(likfun) <- names
  }

  strt <- initpars

  if(ccmethod=="CCAL"){
    fxd <- NULL

    lw <- sapply(names,FUN=function(x) -Inf)
    if(!repar) lw[1] <- -1/max(data[,doses])+.0001
    if(repar){
      up <- list(xi=log(uplimBeta))
    } else {
      up <- list(beta=uplimBeta)
    }

  } else if (ccmethod=="CL"){
    #strt <- strt[names(strt)%in% c("xi","beta", paste0("alpha",2:length(doses)))]
    fxd <- NULL

    lw <- sapply(names[!names %in% names(data)[corrvars]],FUN=function(x) -Inf)
    if(!repar) lw[1] <- -1/max(data[data[,status]==1,doses])+.0001
    if(repar){
      up <- list(xi=log(uplimBeta))
    } else {
      up <- list(beta=uplimBeta)
    }

  } else if(ccmethod=="CCML"){
    if(opt_method=="constrOptim"){
      fxd <- sapply(paste0("alpha",2:length(doses)), function(i) 0)
      lw <- sapply(c(ifelse(repar,"xi","beta"),tail(names, length(corrvars))),FUN=function(x) -Inf)
      if(!repar) {lw[1] <- -1/max(data[,doses][cbind(1:nrow(data), data[,loc])])+.0001}
      #if(length(corrvars)==0){
      if(repar){
        up <- list(xi=log(uplimBeta))
      } else {
        up <- list(beta=uplimBeta)
      }
    } else if(opt_method=="Brent") {
      fxd <- NULL
      lw <- list(params=ifelse(repar, -10,  -1/max(data[,doses][cbind(1:nrow(data), data[,loc])])+.0001))
      up <- list(params=ifelse(repar,log(uplimBeta),uplimBeta))
    } #else { # Nelder-Mead
    #   fxd <- sapply(paste0("alpha",2:length(doses)), function(i) 0)
    #   lw <-NULL
    #   up <- NULL
    #}
    #} else{
    #  up <- NULL
    #}
  } else if(ccmethod=="meandose"){
    if(opt_method=="constrOptim"){
      fxd <- sapply(paste0("alpha",2:length(doses)), function(i) 0)
      lw <- sapply(c(ifelse(repar,"xi","beta"),tail(names, length(corrvars))),FUN=function(x) -Inf)
      if(!repar) {lw[1] <- -1/max(rowMeans(data[,doses]))+.0001}
      #if(length(corrvars)==0){
      if(repar){
        up <- list(xi=log(uplimBeta))
      } else {
        up <- list(beta=uplimBeta)
      }
    } else if(opt_method=="Brent") {
      fxd <- NULL
      lw <- list(params=ifelse(repar, -10,  -1/max(data[,doses][cbind(1:nrow(data), data[,loc])])+.0001))
      up <- list(params=ifelse(repar,log(uplimBeta),uplimBeta))
    }# else { # Nelder-Mead
    #   fxd <- sapply(paste0("alpha",2:length(doses)), function(i) 0)
    #   up <- NULL
    #   lw <- NULL
    # }
    #} else{
    #  up <- NULL
    #}
  }



  if(opt_method=="constrOptim"){
    strt0 <- strt[!names(strt)%in%names(fxd)]
    strt[names(fxd)] <- unlist(fxd)
    if(repar){
      #fit <- mle2(likfun, start=strt,fixed=fxd, method=opt_method,control=fitopt, vecpar=TRUE, parnames=names)
      ui <- matrix(c(-1,rep(0,length(strt0)-1)),nrow=1)
      ci <- -1*log(uplimBeta)
    } else {
      ui <- matrix(c(1,rep(0,length(strt0)-1)),nrow=1)
      ci <- lw[1]
    }
    betathresh <- ci

    fit <- constrOptim(strt0,function(x){
      y <- strt
      y[names(x)] <- x
      likfun(y)
    }, ui=ui, ci=ci, grad=function(x){
      y <- strt
      y[names(x)] <- x
      if(!is.null(fxd)){
        res <-  -scorefun(y[-which(names(y) %in% names(fxd))])$U
        if(any(names(res) %in% names(fxd))) res <- res[-which(names(res) %in% names(fxd))]
      } else {
        res <-  -scorefun(y)$U
      }
      res
    }, control=fitopt, outer.iterations=200, hessian=TRUE)
    tmpcoef <- strt
    tmpcoef[names(fit$par)] <- fit$par
    #fit$par <- tmpcoef

    fit$fullcoef <- c(tmpcoef[1],tmpcoef[-1])
    names(fit$fullcoef) <- names
    fit$betathresh <- betathresh
  } else {
    betathresh <- lw
    fit <- optim(initpars, function(x){
      likfun(c(x, rep(0, length(doses)+length(corrvars)-1)))
    }, method=opt_method, lower=lw, upper=up, control=fitopt, hessian=TRUE)

    fit$fullcoef <- fit$par
    fit$betathresh <- betathresh
  }


  #fit <- mle2(likfun, start=initpars, vecpar=(opt_method!="Brent"),fixed=fxd,upper=up, lower=lw, control=fitopt, method=opt_method)

  if(fitNull){
    repar0 <- repar
    repar <- FALSE
    names[1] <- "beta"
    if(opt_method!="Brent" ){


      fxd2 <- c(fxd, list(beta=0))


      strt <- sapply(names,FUN=function(x) 0)
      strt0 <- strt[!names(strt)%in%names(fxd2)]
      strt[names(fxd2)] <- unlist(fxd2)
      ui <- matrix(c(-1,rep(0,length(strt0)-1)),nrow=1)
      ci <- -20

      nullfit <- constrOptim(strt0,function(x){
        y <- strt
        y[names(x)] <- x
        likfun(y)
      }, ui=ui, ci=ci, grad=function(x){
        y <- strt
        y[names(x)] <- x
        if(!is.null(fxd)){
          res <-  -scorefun(y[-which(names(y) %in% names(fxd))])$U
          res <- res[-which(names(res) %in% names(fxd2))]
        } else {
          res <-  -scorefun(y)$U
          res <- res[-which(names(res) %in% c("beta"))]
        }

      }, control=fitopt, outer.iterations=200, hessian=TRUE)




    } else {
      nullfit <- list(par=0, value=likfun(c(0, rep(0, length(doses)+length(corrvars)-1))), hessian=hessian(function(x) likfun(c(x, rep(0, length(doses)+length(corrvars)-1))), 0))

    }

    repar <- repar0
    names[1] <- ifelse(repar, "xi","beta")

    #
    # lw2 <- lw[-1]
    # up2 <- up[-1]
    # if(length(lw2)==0) lw2 <- NULL
    # if(length(up2)==0) up2 <- NULL
    # nullfit <- mle2(likfun, start=initpars, vecpar=(opt_method!="Brent"),fixed=fxd2,upper=up2, lower=lw2, control=fitoptnull, method=opt_method)

  } else {
    nullfit <- NULL
  }

  proflik <- function(x){
    if(opt_method!="Brent"){

      if(repar){
        fxd3 <- c(fxd, list(xi=x))
      } else {
        fxd3 <- c(fxd, list(beta=x))
      }
      strt <- sapply(names,FUN=function(x) 0)
      strt0 <- strt[!names(strt)%in%names(fxd3)]
      strt[names(fxd3)] <- unlist(fxd3)
      ui <- matrix(c(-1,rep(0,length(strt0)-1)),nrow=1)
      ci <- -20

      proffit <- constrOptim(strt0,function(xx){
        y <- strt
        y[names(xx)] <- xx
        likfun(y)
      }, ui=ui, ci=ci, grad=function(xx){
        y <- strt
        y[names(xx)] <- xx
        if(!is.null(fxd)){
          res <-  -scorefun(y[-which(names(y) %in% names(fxd))])$U
          res <- res[-which(names(res) %in% names(fxd3))]
        } else {
          res <-  -scorefun(y)$U
          res <- res[-which(names(res) %in% c("xi","beta"))]
        }

      }, control=fitopt, outer.iterations=200, hessian=TRUE)
    } else {
      proffit <- list(par=x, value=likfun(c(x, rep(0, length(doses)+length(corrvars)-1))), hessian=hessian(function(xxx) likfun(c(xxx, rep(0, length(doses)+length(corrvars)-1))), x))

      # fitoptprof <- fitopt
      # fxd3 <- list(params=x)
    }
    # lw3 <- lw[-1]
    # up3 <- up[-1]
    # if(length(lw3)==0) lw3 <- NULL
    # if(length(up3)==0) up3 <- NULL


    proffit$value
    #mle2(likfun, start=initpars, vecpar=(opt_method!="Brent"),fixed=fxd3,upper=up3, lower=lw3, control=fitoptprof, method=opt_method)@min
  }


  list(fit=fit, nullfit=nullfit, proflik=proflik)

}
