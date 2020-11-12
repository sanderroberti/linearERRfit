#' @title Negative log-likelihood
#' @description Compute the negative log-likelihood for the linear ERR model in a matched case-control dataset. This is an outdated function, replaced by the substantially faster linERRloglik
#' @param params vector of parameter values (\eqn{\beta}, \eqn{\alpha_2}, ... , \eqn{\alpha_L}, \eqn{\gamma_1}, ... , \eqn{\gamma_p})  to evaluate the log-likelihood at
#' @param data data frame containing matched case-control data, with a number of columns for doses to different locations, a column containing matched set numbers, a column containing the case's tumor location (value between 1 and the number of locations, with location \eqn{x} corresponding to the \eqn{x}-th column index in \code{doses}) and a column serving as a case-control indicator. Other covariates can also be included, in this case a parameter for each covariate column will be estimated. Hence factor variables need to be converted to dummy variables using \code{model.matrix}. If using \code{ccmethod='meandose'}, a column for tumor location is still required but in this case the column can be a vector of ones.
#' @param doses vector containing the indices of columns containing dose information, in the desired order.
#' @param set column index containing matched set numbers.
#' @param status column index containing case status.
#' @param loc column index containing the location of the matched set's case's second tumor.
#' @param corrvars vector containing the indices of columns containing variables to be corrected for.
#' @param ccmethod choice of method of analysis: one of meandose, CCML, CCAL or CL. Defaults to CCAL
#' @return Minus log likelihood in params
#' @examples
#' data(linearERRdata1)
#'
#' #log-likelihood in the truth
#' -linERRloglikold(params=c(.3,rep(log(.25),4),log(2)),
#' data=linearERRdata1,set=1, doses=2:6, status=8, loc=7, corrvars=9)
#'
#' #log-likelihood in 0
#' -linERRloglikold(params=c(0,rep(0,4),0),
#' data=linearERRdata1,set=1, doses=2:6, status=8, loc=7, corrvars=9)

#' @importFrom utils tail
#' @export

linERRloglikold <- function(params, data, doses, set, status, loc, corrvars=NULL, ccmethod="CCAL"){ # -log(likelihood)
  beta0 <- params[1]

  alpha <- c(0,params[2:(length(doses))])
  delta <- tail(params, length(corrvars))

  sp <- split(data,data[,set])
  sp <- lapply(sp, as.matrix)

  byset <- unlist(lapply(sp, function(x){
    if(ccmethod=="meandose"){
      num <- as.numeric((1+rowMeans(x[,doses, drop=FALSE])*beta0)[x[,status]==1]*ifelse(length(corrvars)>0,exp(x[x[,status]==1, corrvars]%*%delta),1))
      den <- sum(as.numeric((1+rowMeans(x[,doses, drop=FALSE])*beta0)*ifelse(length(corrvars)>0,list((exp(x[,corrvars,drop=FALSE]%*%delta))),1)[[1]]))
    } else if(ccmethod=="CCML"){
      num <- as.numeric((1+(x[x[,status]==1,doses, drop=FALSE]*beta0)[x[x[,status]==1, loc]])*ifelse(length(corrvars)>0,exp(x[x[,status]==1, corrvars]%*%delta),1))
      den <- sum(as.matrix((1+(x[,doses, drop=FALSE]*beta0)[,x[x[,status]==1, loc]])*ifelse(length(corrvars)>0,list((exp(x[,corrvars,drop=FALSE]%*%delta))),1)[[1]]))
    } else if(ccmethod=="CCAL"){
      num <- (as.numeric((1+(x[x[,status]==1,doses]*beta0))*exp(alpha)*ifelse(length(corrvars)>0,exp(x[x[,status]==1, corrvars]%*%delta),1))[x[x[,status]==1, loc]])
      den <- sum(as.matrix(1+(x[,doses]*beta0), nrow=nrow(x))*matrix(rep(exp(alpha),nrow(x)),nrow=nrow(x),byrow=TRUE)*ifelse(length(corrvars)>0,list(matrix(rep(exp(x[,corrvars,drop=FALSE]%*%delta),length(doses)), nrow=nrow(x))),1)[[1]])
    } else if(ccmethod=="CL"){
      num <- (as.numeric((1+(x[x[,status]==1,doses]*beta0))*exp(alpha))[x[x[,status]==1, loc]])
      den <- sum(as.numeric(1+(x[x[,status]==1,doses]*beta0))*exp(alpha))
    } else if(ccmethod=="OED"){
      num <- mean((as.numeric((1+(x[x[,status]==1,doses]*beta0))*exp(alpha)*ifelse(length(corrvars)>0,exp(x[x[,status]==1, corrvars]%*%delta),1))))
      den <- sum((as.numeric(rowMeans((1+(x[,doses]*beta0))*matrix(rep(exp(alpha),nrow(x)),nrow=nrow(x),byrow=TRUE))*ifelse(length(corrvars)>0,list((exp(x[,corrvars,drop=FALSE]%*%delta))),1)[[1]])))
    }

    return(ifelse(num!=0,log(num/den),0))

  }))

  -sum(byset)
}
