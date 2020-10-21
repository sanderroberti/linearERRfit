#' @title Fit linear ERR model and perform jackknife correction
#' @description Fits the linear ERR model on a dataset and performs first and second order jackknife correction
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
#' @param uplimBeta upper limit for beta=exp(xi)
#' @param profCI boolean: compute 95\% profile likelihood confidence interval for beta/xi?
#' @param doJK1 perform first order jackknife correction?
#' @param doJK2 perform second order jackknife correction?
#' @param jkscorethresh square L2 norm threshold for leave-one-out and leave-two-out estimates to be included in the computation of the first and second order jackknife corrected estimate, respectively
#' @param jkvalrange range of leave-one-out and leave-two-out beta/xi estimates to be allowed in the computation of the first and second order jackknife corrected estimate, respectively
#' @return object with components 'MLE' and 'jackknife'. The 'jackknife' component in turn has components 'firstorder' and 'secondorder'
#' @importFrom stats model.matrix pchisq uniroot
#' @importFrom utils combn
#' @examples
#' data(linearERRdata1)
#'
#' fitCCML <- linearERR(data=linearERRdata1, set=1, doses=2:6, status=8,
#' loc=7, corrvars=9, repar=FALSE, ccmethod="CCML", doJK1=TRUE)
#'
#' fitCCML$MLE
#' fitCCML$jackknife$firstorder$coef

#' @export
#'


linearERR <- function(data, doses, set, status, loc, corrvars=NULL, ccmethod="CCAL", repar=FALSE, initpars=rep(0,length(doses)+length(corrvars)), fitopt=NULL, uplimBeta=5,profCI=TRUE,doJK1=FALSE,doJK2=FALSE,jkscorethresh=.01,jkvalrange=c(-Inf, Inf)){

  if(doJK2) doJK1 <- TRUE

  mainfit <- linearERRfit(data=data, doses=doses, set=set, status=status, loc=loc, corrvars=corrvars, repar=repar, initpars=initpars, fitopt=fitopt, ccmethod=ccmethod, uplimBeta=uplimBeta)

  MLEscore <- linERRscore(params=mainfit$fit@coef, data=data, doses=doses, set=set, status=status, loc=loc, corrvars=corrvars, repar=repar,ccmethod=ccmethod)$U

  pval <- pchisq(2*(mainfit$nullfit@min-mainfit$fit@min), df=1, lower.tail=FALSE)


  if(profCI){


    g <- function(para){
      1-pchisq(2*(mainfit$proflik(para)-mainfit$fit@min),df=1)-.05
    }
    lowLim <- tryCatch(uniroot(g, lower=ifelse(repar,-20,mainfit$fit@call$lower[1]), upper=mainfit$fit@coef[1], extendInt="no")$root, error=function(e) NA)
    upLim <- tryCatch(uniroot(g, lower=mainfit$fit@coef[1],upper=ifelse(repar,log(100),100), extendInt="no", maxiter=150)$root, error=function(e) NA)

  } else{
    lowLim <- NULL
    upLim <- NULL
  }

  MLE <- list(coef=mainfit$fit@coef,sd=sqrt(diag(mainfit$fit@vcov)), vcov=mainfit$fit@vcov, score=MLEscore, convergence=mainfit$fit@details$convergence, message=mainfit$fit@details$message, dosepval=pval, profCI=c(lo=lowLim, up=upLim))

  # Jackknife

  setnrs <- unique(data[,set])


  if(doJK1){

    outfunjk1 <- function(exclset) {
      jk1fit <- linearERRfit(data=data[!(data[,set]== exclset),], doses=doses, set=set, status=status, loc=loc, corrvars=corrvars, repar=repar, initpars=initpars, fitopt=fitopt, ccmethod=ccmethod, fitNull=FALSE, uplimBeta=uplimBeta)
      jk1score <- tryCatch(linERRscore(params=jk1fit$fit@coef, data=data[!(data[,set]== exclset),], doses=doses, set=set, status=status, loc=loc, corrvars=corrvars, repar=repar,ccmethod=ccmethod)$U, error=function(e) rep(NA, length(jk1fit$fit@fullcoef)))
      list(fit=jk1fit,score=jk1score)
    }
    jk1out <- lapply(setnrs, function(k) outfunjk1(k))

    jk1coefs <- sapply(jk1out, function(x) x$fit$fit@coef)
    jk1coefs <- as.data.frame(matrix(jk1coefs, nrow=length(setnrs),byrow=TRUE, dimnames=list(NULL,names(jk1out[[1]]$fit$fit@coef))))
    jk1scores <- sapply(jk1out, function(x) x$score)
    jk1scores <- as.data.frame(matrix(jk1scores, nrow=length(setnrs),byrow=TRUE, dimnames=list(NULL,names(jk1out[[1]]$fit$fit@coef))))
    jk1conv <- sapply(jk1out, function(x) x$fit$fit@details$convergence)

    jk1included <- (1-1*(jk1coefs[,1]==ifelse(repar,log(uplimBeta),uplimBeta)))*(rowSums(jk1scores^2)<jkscorethresh)*(1-(jk1conv==1))*(jk1coefs[,1]>=jkvalrange[1])*(jk1coefs[,1]<=jkvalrange[2])

    jk1coef <- length(setnrs)*mainfit$fit@coef-(length(setnrs)-1)*colMeans(jk1coefs[jk1included==1,])

    jk1details <- data.frame(cbind(setnrs, included=jk1included,conv=jk1conv,coef=jk1coefs,score=jk1scores))

    jackknife1 <- list(coef=jk1coef, details=jk1details)

  } else {
    jackknife1 <- NULL
  }




  if(doJK2){
    allpairs <- t(combn(setnrs,2))

    outfunjk2 <- function(pair) {
      jk2fit <- linearERRfit(data=data[!(data[,set]%in% pair),], doses=doses, set=set, status=status, loc=loc, corrvars=corrvars, repar=repar, initpars=initpars, fitopt=fitopt, ccmethod=ccmethod, fitNull=FALSE, uplimBeta=uplimBeta)
      jk2score <- tryCatch(linERRscore(params=jk2fit$fit@coef, data=data[!(data[,set]%in% pair),], doses=doses, set=set, status=status, loc=loc, corrvars=corrvars, repar=repar,ccmethod=ccmethod)$U, error=function(e) rep(NA, length(jk2fit$fit@fullcoef)))
      list(fit=jk2fit, score=jk2score)
    }
    jk2out <- lapply(1:nrow(allpairs), function(k) outfunjk2(allpairs[k,]))

    jk2coefs <- sapply(jk2out, function(x) x$fit$fit@coef)
    jk2coefs <- as.data.frame(matrix(jk2coefs, nrow=nrow(allpairs),byrow=TRUE, dimnames=list(NULL,names(jk2out[[1]]$fit$fit@coef))))
    jk2scores <- sapply(jk2out, function(x) x$score)
    jk2scores <- as.data.frame(matrix(jk2scores, nrow=nrow(allpairs),byrow=TRUE, dimnames=list(NULL,names(jk2out[[1]]$fit$fit@coef))))
    jk2conv <- sapply(jk2out, function(x) x$fit$fit@details$convergence)

    jk2included <- (1-1*(jk2coefs[,1]==ifelse(repar,log(uplimBeta),uplimBeta)))*(rowSums(jk2scores^2)<jkscorethresh)*(1-(jk2conv==1))*(jk2coefs[,1]>=jkvalrange[1])*(jk2coefs[,1]<=jkvalrange[2])

    jk2coef <- (length(setnrs)^3*mainfit$fit@coef-(2*length(setnrs)^2-2*length(setnrs)+1)*(length(setnrs)-1)*colMeans(jk1coefs[jk1included==1,, drop=FALSE])+(length(setnrs)-1)^2*(length(setnrs)-2)*colMeans(jk2coefs[jk2included == 1,, drop=FALSE]))/(2*length(setnrs)-1)

    allpairs <- as.data.frame(allpairs)
    names(allpairs) <- c("set1","set2")
    jk2details <- data.frame(cbind(allpairs, included=jk2included,conv=jk2conv,coef=jk2coefs,score=jk2scores))

    jackknife2 <- list(coef=jk2coef, details=jk2details)

  } else {
    jackknife2 <- NULL
  }

  jackknife <- list(firstorder=jackknife1, secondorder=jackknife2)

  list(MLE=MLE, jackknife=jackknife)

}

