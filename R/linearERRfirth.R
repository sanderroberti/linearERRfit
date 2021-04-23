#' @title Derive the Firth-corrected estimate for the linear ERR model
#' @description Finds roots to the Firth-corrected score equations for the linear ERR model using a matched case-control study.
#' @param data data frame containing matched case-control data, with a number of columns for doses to different locations, a column containing matched set numbers, a column containing the case's tumor location (value between 1 and the number of locations, with location \eqn{x} corresponding to the \eqn{x}-th column index in \code{doses}) and a column serving as a case-control indicator. Other covariates can also be included, in this case a parameter for each covariate column will be estimated. Hence factor variables need to be converted to dummy variables using \code{model.matrix}. If using \code{ccmethod='meandose'}, a column for tumor location is still required but in this case the column can be a vector of ones.
#' @param doses vector containing the indices of columns containing dose information.
#' @param set column index containing matched set numbers.
#' @param status column index containing case status.
#' @param loc column index containing the location of the matched set's case's second tumor.
#' @param corrvars vector containing the indices of columns containing variables to be corrected for. Not used with \code{ccmethod='CL'}
#' @param repar reparametrize to \eqn{\beta=exp(\xi)}? It is recommended to reparametrize when using CL or CCAL or when using additional covariates. Defaults to \code{FALSE}
#' @param ccmethod choice of method of analysis: one of meandose, CCML, CCAL or CL. Defaults to CCAL
#' @param initpars initial values for parameters, default is 0 for all parameters. If supplying a different vector, use a vector with an initial value for all free parameters (\eqn{\beta} or \eqn{\xi}, one for each location effect (except the reference) when using CL or CCAL, and for each other covariate if applicable, in that order). Note that if \code{repar=TRUE}, the first initial value is used for \eqn{\xi}.
#' @param lowerlim lower bound for model parameters, in the same order as \code{initpars}. At least one upper or lower limit needs to be finite. Note that when \code{repar=TRUE}, the first entry is the lower limit for \eqn{\xi}. When \code{repar=FALSE}, the lower limit for \eqn{\beta} cannot be smaller than \eqn{-1/max(d)} where the maximum is taken among all relevant doses for the chosen \code{ccmethod}. If this is the case, the limit will automatically be changed to that value
#' @param upperlim upper bound for model parameters, in the same order as \code{initpars}. At least one upper or lower limit needs to be finite. Note that when \code{repar=TRUE}, the first entry is the upper limit for \eqn{\xi}. When \code{repar=TRUE}, if no other lower or upper limit is given as input, an upper limit of \eqn{\log(5)} will be used for \eqn{\xi}
#' @param fitopt list with options to pass to \code{control} argument of optimizer (see details)
#' @return \code{optim} object with fit results.
#' @references David Firth, Bias reduction of maximum likelihood estimates, Biometrika, Volume 80, Issue 1, March 1993, Pages 27–38, \href{https://doi.org/10.1093/biomet/80.1.27}{https://doi.org/10.1093/biomet/80.1.27}
#' @importFrom stats optim constrOptim
#' @details This function looks for roots of the Firth-corrected score functions.
#'
#' The underlying model is HR=\eqn{\sum(1+\beta d_l)exp(\alpha_l+X^T\gamma)}, where the sum is over organ locations. Here \eqn{\beta} is the dose effect, \eqn{\alpha} are the location effects and \eqn{\gamma} are other covariate effects. The model can be reparametrized to HR=\eqn{\sum(1+exp(\xi) d_l)exp(\alpha_l+X^T\gamma)} using \code{repar=TRUE}. In the original parametrization, \eqn{\beta} is constrained such that HR cannot be negative. There are different choices for the design used to estimate the parameters: mean organ dose, CCML, CL, and CCAL. Mean organ dose (\code{ccmethod='meandose'}) uses the mean of the supplied location doses and compares that mean dose between case and matched controls. The other choices (CCML, CL and CCAL) use the tumor location for the case and compare either only between patients (CCML), only within patients (CL) or both between and within patients (CCAL). CCML only compares the same location between patients, and hence cannot be used to estimate location effects. Similarly, CL compares within patients and cannot be used to estimate covariate effects other than dose, meaning \code{corrvars} should not be supplied for CL. For this model, the Firth correction (Firth 1993) is used as a method for bias correction, or for obtaining an estimate when there is separation in the data.
#'
#' To avoid using unstable multidimensional root finders, this function minimizes the square L2 norm of the modified score instead. This is done using the \code{optim} function. If desired, it is possible to use \code{linERRscore} and optimize or search for roots directly. For one-dimensional models (i.e., mean dose or CCML without additional covariates), the Brent algorithm is used with the user-supplied search interval (\code{lowerlim},\code{upperlim}). Note that the choice for search interval is crucial as this determines convergence. For this reason, there is no default setting in this case. For other optimizations, the L-BFGS-B algorithm (with constraints \code{lowerlim} and \code{upperlim}) is used. For details refer to the function optim, also for \code{fitopt} settings. When \code{repar=FALSE}, if the lower bound for \eqn{\beta} is set too small, it is automatically changed according to the positivity constraint for HR.
#'
#' It is advisable to interpret the results with caution. It was found that the modified score function sometimes has multiple roots, which makes setting initial values and search intervals crucial. It is recommended to try different settings for these inputs. Further, it seemed that reparametrizing improved the performance for multidimensional models.
#' @examples
#' data(linearERRdata1)
#'
#' fitMLE <- linearERR(data=linearERRdata1,doses=2:6,set=1,status=8,loc=7,
#' corrvars=9,repar=TRUE,ccmethod="CCAL",profCI=FALSE)
#'
#' fitfirth <- linearERRfirth(data=linearERRdata1,doses=2:6,set=1,status=8,loc=7,
#' corrvars=9,repar=TRUE,ccmethod="CCAL",initpars=fitMLE$MLE$coef)
#'
#' data.frame(MLE=fitMLE$MLE$coef, Firth=fitfirth$par)
#'
#' @export


linearERRfirth <- function(data, doses, set, status, loc, corrvars=NULL, repar=FALSE, ccmethod="CCAL", initpars=NULL,lowerlim=NULL, upperlim=NULL, fitopt=list(maxit=5000)){

  if(ccmethod=="CL" & !is.null(corrvars)) stop("corrvars needs to be set to NULL when using CL")
  #if(ccmethod=="CL") corrvars <- NULL


  if(is.null(initpars)){
    if(ccmethod %in%c("CCML","meandose")) initpars <- rep(0,1+length(corrvars))
    if(ccmethod %in%c("CCAL","CL")) initpars <- rep(0,length(doses)+length(corrvars))
  }

  if(is.null(lowerlim)) lowerlim <- rep(-Inf, length(initpars))
  if(is.null(upperlim)) {
    upperlim <- rep(Inf, length(initpars))
    if(repar) upperlim[1] <- log(5)
  }


  if(ccmethod%in%c("CCML","meandose") & length(initpars)!=(length(corrvars)+1)) stop("Length of initpars incorrect. Please provide an initial value for the dose effect and one for each of the other covariates.")
  if(ccmethod%in%c("CCML","meandose") & length(lowerlim)!=(length(corrvars)+1)) stop("Length of lowerlim incorrect. Please provide a lower bound for the dose effect and one for each of the other covariates.")
  if(ccmethod%in%c("CCML","meandose") & length(upperlim)!=(length(corrvars)+1)) stop("Length of upperlim incorrect. Please provide an upper bound for the dose effect and one for each of the other covariates.")

  if(ccmethod%in%c("CL","CCAL") & length(initpars)!=(length(corrvars)+length(doses))) stop("Length of initpars incorrect. Please provide an initial value for the dose effect, for all but one of the alphas, and one for each of the other covariates if used.")
  if(ccmethod%in%c("CL","CCAL") & length(lowerlim)!=(length(corrvars)+length(doses))) stop("Length of lowerlim incorrect. Please provide a lower bound for the dose effect, for all but one of the alphas, and one for each of the other covariates if used.")
  if(ccmethod%in%c("CL","CCAL") & length(upperlim)!=(length(corrvars)+length(doses))) stop("Length of upperlim incorrect. Please provide an upper bound for the dose effect, for all but one of the alphas, and one for each of the other covariates if used.")




  scorefun <- function(params){
    linERRscore(params=params,data=data, doses=doses, set=set, status=status, loc=loc, ccmethod=ccmethod, corrvars=corrvars, repar=repar)
  }


  if(ccmethod %in% c("CCML", "meandose") & length(corrvars)==0){
    opt_method <- "Brent"
  } else{
    opt_method <- "constrOptim"
  }

  if(opt_method=="Brent" & length(initpars)>1) initpars=initpars[1]

  lw <- lowerlim
  up <- upperlim
  if(opt_method=="Brent" & any(is.infinite(c(lw,up)))) stop("Please provide finite search bounds for Brent")

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
  #   fitopt <- c(fitopt, list(ndeps=rep(1e-5, length(initpars))))
  # }



  if(ccmethod=="CCAL"){
    if(!repar) lw[1] <- max(lw[1],-1/max(data[,doses])+.0001)
  } else if (ccmethod=="CL"){
    if(!repar) lw[1] <- max(lw[1],-1/max(data[data[,status]==1,doses])+.0001)
  } else if(ccmethod=="CCML"){
    if(!repar) {lw[1] <- max(lw[1],-1/max(data[,doses][cbind(1:nrow(data), data[,loc])])+.0001)}

  } else if(ccmethod=="meandose"){
    if(!repar) {lw[1] <- max(lw[1],-1/max(rowMeans(data[,doses]))+.0001)}
  }

  if(opt_method=="Brent"){

    fit <- optim(initpars,
                 function(x){
                   tmp <- scorefun(x)
                   sum((tmp$U+tmp$A)^2)
                 },method=opt_method, control=fitopt, lower=lw, upper=up)
  } else {

    ui <- rbind(-1*diag(length(initpars)), diag(length(initpars)))
    ci <- c(-1*up, lw)

    ui <- ui[which(!is.infinite(ci)),, drop=FALSE]
    ci <- ci[which(!is.infinite(ci))]

    fit <- constrOptim(initpars,
                       function(x){
                         tmp <- scorefun(x)
                         sum((tmp$U+tmp$A)^2)
                       }, grad=NULL, control=fitopt, ci=ci,ui=ui, outer.iterations = 200)
  }


  return(fit)

}
