#' @title Negative log-likelihood
#' @description Compute the negative log-likelihood for the linear ERR model in a matched case-control dataset
#' @param params vector of parameter values (\eqn{\beta}, \eqn{\alpha_2}, ... , \eqn{\alpha_L}, \eqn{\gamma_1}, ... , \eqn{\gamma_p}) to evaluate the log-likelihood at
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
#' -linERRloglik(params=c(.3,rep(-1.386294,4),log(2)),
#' data=linearERRdata1,set=1, doses=2:6, status=8, loc=7, corrvars=9)
#'
#' #log-likelihood in 0
#' -linERRloglik(params=c(0,rep(0,4),0),
#' data=linearERRdata1,set=1, doses=2:6, status=8, loc=7, corrvars=9)

#' @importFrom utils tail
#' @importFrom stats reshape
#' @export

linERRloglik <- function(params, data, doses, set, status, loc, corrvars=NULL, ccmethod="CCAL"){ # -log(likelihood)
  #if(!(length(params)==length(doses)+length(corrvars))) stop("Mismatch between length of params and specified columns in data")

  data <- data[data[,set] %in% as.numeric(names(which(table(data[,set],data[,status])[,2]==1))),]
  if(ccmethod!="CL") data <- data[data[,set] %in% as.numeric(names(which(table(data[,set])>1))),]

  beta0 <- params[1]

  # alpha are location effects
  if(ccmethod %in% c("CCAL","CL")){
    alpha <- c(0,params[2:(length(doses))])
  } else {
    alpha <- rep(0,length(doses))
  }

  # gamma are other patient-level effect parameters
  gamma <- tail(params, length(corrvars))

  if(!ccmethod=="meandose"){

    # rrmat is a matrix with relative risk for each location (columns) from each patient (rows)
    rrmat <- matrix(rep(exp(as.matrix(data[,corrvars,drop=FALSE])%*%gamma), length(doses)), ncol=length(doses))*(1+beta0*data[,doses])*matrix(rep(exp(alpha), nrow(data)), ncol=length(doses), byrow=TRUE)

    if(ccmethod=="CCAL"){
      # CCAL

      rrcases <- rrmat[data[,status]==1,][cbind(1:sum(data[,status]), data[data[,status]==1,loc])]

      # reshape rrmat to wide format (with each row representing an entire matched set of locations) for efficient summation
      rrdf <- data.frame(set=data[,set], id=1,rrmat)
      rrdf$id <- with(rrdf, ave(as.numeric(set), set, FUN = seq_along))
      rrdfwide <- reshape(rrdf, direction="wide", idvar="set",timevar="id")

      minuslogl <- -(sum(log(rrcases)) - sum(log(rowSums(rrdfwide[,-1], na.rm=TRUE))))
    } else if(ccmethod=="CL"){

      # CL

      rrcases <- rrmat[data[,status]==1,][cbind(1:sum(data[,status]), data[data[,status]==1,loc])]
      rrdf <- data.frame(set=data[,set],rrmat)
      rrdf <- rrdf[data[status]==1,]

      minuslogl <- -(sum(log(rrcases)) - sum(log(rowSums(rrdf[,-1]))))
    } else if(ccmethod=="CCML"){
      # CCML
      rrcases <- rrmat[data[,status]==1,][cbind(1:sum(data[,status]), data[data[,status]==1,loc])]

      # reshape rrmat to wide format (with each row representing an entire matched set of locations) for efficient summation
      rrdf <- data.frame(set=data[,set], id=1, rrmat[cbind(1:nrow(rrmat), data[,loc])])
      rrdf$id <- with(rrdf, ave(as.numeric(set), set, FUN = seq_along))
      rrdfwide <- reshape(rrdf, direction="wide", idvar="set",timevar="id")

      minuslogl <- -(sum(log(rrcases)) - sum(log(rowSums(rrdfwide[,-1], na.rm=TRUE))))
    }
  } else {
    # Mean dose method, rrmat based on mean dose and patient-level effects

    rrmat <- (1+beta0*rowMeans(data[,doses]))*exp(as.matrix(data[,corrvars,drop=FALSE])%*%gamma)
    rrdf <- data.frame(set=data[,set], id=1, rrmat)
    rrdf$id <- with(rrdf, ave(as.numeric(set), set, FUN = seq_along))
    rrdfwide <- reshape(rrdf, direction="wide", idvar="set",timevar="id")

    minuslogl <- -(sum(log(rrmat[data[,status]==1,]))-sum(log(rowSums(rrdfwide[,-1], na.rm=TRUE))))

  }

  minuslogl
}
