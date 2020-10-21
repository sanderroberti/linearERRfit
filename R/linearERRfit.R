#' @title Fit linear ERR model
#' @description Fits the linear ERR model on a dataset
#' @param data data frame containing matched case-control data, with a number of columns for doses to different locations, a column containing matched set numbers, a column containing the case's tumor location and a column serving as a case-control indicator. Other covariates can also be included.
#' @param doses vector containing the indices of columns containing dose information, in the desired order.
#' @param set column index containing matched set numbers.
#' @param status column index containing case status.
#' @param loc column index containing the location of the matched set's case's second tumor.
#' @param corrvars vector containing the indices of columns containing variables to be corrected for.
#' @param repar reparametrize to beta=exp(xi)
#' @param ccmethod choice of method of analysis; one of meandose, CCML, CCAL or CL.
#' @param initpars initial values for parameters
#' @param fitopt list with options to pass to 'control' argument of optimizer
#' @param fitNull boolean: also fit model without dose effect?
#' @param useOld if TRUE, a previous (slower) implementation of the log-likelihood function will be used
#' @param uplimBeta upper limit for beta=exp(xi)
#' @return object with components 'MLE' and 'jackknife'. The 'jackknife' component in turn has components 'firstorder' and 'secondorder'
#' @importFrom bbmle mle2 parnames parnames<- coef
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
#' status=8, loc=7, corrvars=9, repar=FALSE, ccmethod="CL")
#'
#' bbmle::coef(fitmeandose$fit, exclude.fixed=TRUE)
#' bbmle::coef(fitCCML$fit, exclude.fixed=TRUE)
#' bbmle::coef(fitCCAL$fit, exclude.fixed=TRUE)
#' bbmle::coef(fitCL$fit, exclude.fixed=TRUE)
#' @export


linearERRfit <- function(data, doses, set, status, loc, corrvars=NULL, repar=FALSE, ccmethod="CCAL", initpars=rep(0,length(doses)+length(corrvars)), fitopt=list(maxit=5000), fitNull=TRUE, useOld=FALSE, uplimBeta=5){

  if(ccmethod %in% c("CCML", "meandose") & length(corrvars)==0){
    opt_method <- "Brent"
  } else{
    opt_method <- ifelse(repar, "Nelder-Mead", "L-BFGS-B")
  }

  if(is.null(fitopt$maxit)){
    fitopt <- c(fitopt, list(maxit=5000))
  }
  if(is.null(fitopt$reltol) & opt_method!="L-BFGS-B"){
    fitopt <- c(fitopt, list(reltol=1e-10))
  }
  if (is.null(fitopt$pgtol) & opt_method=="L-BFGS-B"){
    fitopt <- c(fitopt, list(pgtol=1e-8))
  }
  if (is.null(fitopt$factr) & opt_method=="L-BFGS-B"){
    fitopt <- c(fitopt, list(factr=1e4))
  }

  likfun <- function(params){
    if(repar) params[1] <- exp(params[1])
    if(useOld){
      linERRloglikold(params, data=data, set=set, doses=doses,status=status,loc=loc,corrvars=corrvars, ccmethod=ccmethod)
    } else {
      linERRloglik(params, data=data, set=set, doses=doses,status=status,loc=loc,corrvars=corrvars, ccmethod=ccmethod)
    }

  }

  if(opt_method!="Brent"){
    names <- c(ifelse(repar,"xi","beta"), paste0("alpha",2:length(doses)),names(data)[corrvars])
    names(initpars) <- names
    parnames(likfun) <- names
  } else {
    names <- "params"
    initpars <- list(params=as.numeric(initpars[1]))
    parnames(likfun) <- names
  }

  if(ccmethod=="CCAL"){
    fxd <- NULL
    if(opt_method!="Nelder-Mead"){
      lw <- sapply(names,FUN=function(x) -Inf)
      if(!repar) lw[1] <- -1/max(data[,doses])+.0001
      if(repar){
        up <- list(xi=log(uplimBeta))
      } else {
        up <- list(beta=uplimBeta)
      }
    } else {
      lw <- NULL
      up <- NULL
    }
  } else if (ccmethod=="CL"){
    fxd <- sapply(names(data)[corrvars], function(i) 0)
    if(opt_method!="Nelder-Mead"){
      lw <- sapply(names[!names %in% names(data)[corrvars]],FUN=function(x) -Inf)
      if(!repar) lw[1] <- -1/max(data[data[,status]==1,doses])+.0001
      if(repar){
        up <- list(xi=log(uplimBeta))
      } else {
        up <- list(beta=uplimBeta)
      }
    } else {
      lw <- NULL
      up <- NULL
    }
  } else if(ccmethod=="CCML"){
    if(opt_method=="L-BFGS-B"){
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
    } else { # Nelder-Mead
      fxd <- sapply(paste0("alpha",2:length(doses)), function(i) 0)
      lw <-NULL
      up <- NULL
    }
    #} else{
    #  up <- NULL
    #}
  } else if(ccmethod=="meandose"){
    if(opt_method=="L-BFGS-B"){
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
    } else { # Nelder-Mead
      fxd <- sapply(paste0("alpha",2:length(doses)), function(i) 0)
      up <- NULL
      lw <- NULL
    }
    #} else{
    #  up <- NULL
    #}
  }


  fit <- mle2(likfun, start=initpars, vecpar=(opt_method!="Brent"),fixed=fxd,upper=up, lower=lw, control=fitopt, method=opt_method)

  if(fitNull){
    if(opt_method!="Brent"){
      if(repar){
        fxd2 <- c(fxd, list(xi=-999999))
      } else {
        fxd2 <- c(fxd, list(beta=0))
      }
    } else {
      fxd2 <- list(params=ifelse(repar, -999999, 0))
    }

    lw2 <- lw[-1]
    up2 <- up[-1]
    if(length(lw2)==0) lw2 <- NULL
    if(length(up2)==0) up2 <- NULL
    nullfit <- mle2(likfun, start=initpars, vecpar=(opt_method!="Brent"),fixed=fxd2,upper=up2, lower=lw2, control=fitopt, method=opt_method)
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
    } else {
      fxd3 <- list(params=x)
    }
    lw3 <- lw[-1]
    up3 <- up[-1]
    if(length(lw3)==0) lw3 <- NULL
    if(length(up3)==0) up3 <- NULL

    mle2(likfun, start=initpars, vecpar=(opt_method!="Brent"),fixed=fxd3,upper=up3, lower=lw3, control=fitopt, method=opt_method)@min
  }


  list(fit=fit, nullfit=nullfit, proflik=proflik)

}
