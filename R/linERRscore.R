#' @title Compute the Firth-corrected score function
#' @description Compute the Firth-corrected score for a matched case-control dataset
#' @param params vector of parameters (\eqn{\beta/\xi}, \eqn{\alpha_2}, ... , \eqn{\alpha_L}, \eqn{\gamma_1}, ... , \eqn{\gamma_p}) for which to compute the modified score. Parameters \eqn{\alpha} only need to be supplied when \code{ccmethod="CL"} or \code{ccmethod="CCAL"}. Note that when \code{repar=TRUE}, \eqn{\xi} needs to be supplied
#' @param data data frame containing matched case-control data, with a number of columns for doses to different locations, a column containing matched set numbers, a column containing the case's tumor location (value between 1 and the number of locations, with location \eqn{x} corresponding to the \eqn{x}-th column index in \code{doses}) and a column serving as a case-control indicator. Other covariates can also be included, in this case a parameter for each covariate column will be estimated. Hence factor variables need to be converted to dummy variables using \code{model.matrix}. If using \code{ccmethod='meandose'}, a column for tumor location is still required but in this case the column can be a vector of ones.
#' @param doses vector containing the indices of columns containing dose information.
#' @param set column index containing matched set numbers.
#' @param status column index containing case status.
#' @param loc column index containing the location of the matched set's case's second tumor.
#' @param corrvars vector containing the indices of columns containing variables to be corrected for.
#' @param ccmethod choice of method of analysis: one of meandose, CCML, CCAL or CL. Defaults to CCAL
#' @param repar reparametrize to \eqn{\beta=exp(\gamma)}? Note that this only affects the Firth modified score, not the original score. Defaults to FALSE.
#' @return List object with components:
#' \item{U}{the original score}
#' \item{A}{the modification to the score}
#' The Firth corrected score is equal to \code{U+A}.
#' @examples
#' data(linearERRdata1)
#'
#' # score in the truth
#' score1 <- linERRscore(params=c(.3,rep(-1.386294,4),log(2)),
#' data=linearERRdata1,set=1, doses=2:6, status=8, loc=7, corrvars=9,ccmethod="CCAL")
#'
#' # score in the truth, reparametrized
#' score1_repar <- linERRscore(params=c(log(.3),rep(-1.386294,4),log(2)),
#' data=linearERRdata1,set=1, doses=2:6, status=8, loc=7, corrvars=9,ccmethod="CCAL", repar=TRUE)
#'
#' score1$U # Original score
#' score1$U+score1$A # Firth score
#'
#' score1_repar$U+score1_repar$A # Firth score under reparametrization
#'
#' # score in 0
#' score2 <- linERRscore(params=c(0,rep(0,4),0),
#' data=linearERRdata1,set=1, doses=2:6, status=8, loc=7, corrvars=9, ccmethod="CCAL")
#'
#' score2$U # Original score
#' score2$U+score2$A # Firth score
#'
#' @export



linERRscore <- function(params,data, doses, set, status, loc, ccmethod, corrvars=NULL, repar=FALSE){

  data <- data[data[,set] %in% as.numeric(names(which(table(data[,set],data[,status])[,2]==1))),]
  if(ccmethod!="CL") data <- data[data[,set] %in% as.numeric(names(which(table(data[,set])>1))),]

  sp <- split(data,data[,set])
  sp <- lapply(sp, as.matrix)

  if(repar==FALSE){
    byset <- lapply(sp, function(x){
      # x is the x-matrix for a single matched set
      x <- as.data.frame(x)

      # transform such that each row corresponds to an organ location
      longx <- reshape(x, direction="long", varying=list(doses), v.names="dose", timevar="loc")

      # construct X: a matrix with D*=d/(1+beta*d), dummies for location effects and all variables in the exponential part of the model
      if(ccmethod=="CCAL"){
        X <- as.matrix(cbind(D=longx$dose/(1+params[1]*longx$dose), model.matrix(~factor(loc)-1, data=longx)[,-1], longx[, names(x)[corrvars], drop=FALSE]))

      } else if(ccmethod=="CL"){
        # CL only uses locations from the case patient
        longx <- longx[longx[,names(data)[status]]==1,]
        X <- as.matrix(cbind(D=longx$dose/(1+params[1]*longx$dose), model.matrix(~factor(loc)-1, data=longx)[,-1]))

      } else if(ccmethod%in%c("CCML")){
        # CCML only uses the tumor location
        longx <- longx[longx$loc==longx[,names(data)[loc]],]
        X <- as.matrix(cbind(D=longx$dose/(1+params[1]*longx$dose), longx[, names(x)[corrvars], drop=FALSE]))

      } else if(ccmethod=="meandose"){
        # Mean dose only uses the mean organ dose
        longx <- cbind(x[, set, drop = FALSE], data.frame(dose = rowMeans(x[, doses])), x[, -c(set, doses), drop = FALSE])
        X <- as.matrix(cbind(D = rowMeans(x[, doses])/(1 + params[1] * rowMeans(x[, doses])), x[, names(x)[corrvars], drop = FALSE]))
      }

      # vector of tumor location probabilities (RR_l/sum(RR))
      p <- (1+params[1]*longx$dose)*exp(c(X[,-1,drop=FALSE]%*%params[-1]))/sum((1+params[1]*longx$dose)*exp(c(X[,-1,drop=FALSE]%*%params[-1])))

      # compute S0 (scalar), SX (vector, one component for each variable), SXX (matrix) and SXXX (array) as written in manuscript. The first component is always for beta
      S0 <- sum((1+params[1]*longx$dose)*exp(c(X[,-1,drop=FALSE]%*%params[-1])))
      SX <- ((1+params[1]*longx$dose)*exp(c(X[,-1,drop=FALSE]%*%params[-1])))%*%X

      SXX <- Reduce('+',lapply(1:nrow(longx), function(j) X[j,]%*%t(X[j,])*(1+params[1]*longx$dose[j])*exp(c(X[,-1,drop=FALSE]%*%params[-1]))[j]))

      SXSX <- t(SX)%*%SX

      SXXX <- array(dim=c(ncol(X),ncol(X),ncol(X)))
      for(k in 1:ncol(X)){
        for(s in 1:ncol(X)){
          for(t in 1:ncol(X)){
            SXXX[k,s,t] <- sum(X[,k]*X[,s]*X[,t]*(1+params[1]*longx$dose)*exp(c(X[,-1,drop=FALSE]%*%params[-1])))
          }
        }
      }

      # Expected value for X and x^2
      EX <- p%*%X
      EX2 <- p%*%X^2

      # Components for EUU, expected value of product of 2 U's
      tmp1 <- EX2[1]-EX[1]^2
      tmp2 <- SXX[1,-1]/S0-SX[1]*SX[-1]/S0^2
      tmp3 <- SXX[-1,-1]/S0- SXSX[-1,-1]/S0^2

      EUU <- rbind(c(tmp1, tmp2),cbind(tmp2,tmp3))


      # Compute EUUU (expected values of product of 3 U's) and EUU2 (expected values of product of U_a and U_{b,c})
      # Definitions of components are given in the manuscript
      EUUU <- array(dim=c(ncol(X),ncol(X),ncol(X)))
      EUU2 <- array(dim=c(ncol(X),ncol(X),ncol(X)))
      for(k in 1:ncol(X)){
        for(s in 1:ncol(X)){
          for(t in 1:ncol(X)){
            if(k+s+t == 3){
              EUUU[1,1,1] <- SXXX[1,1,1]/S0-3*(SXX[1,1]/S0)*SX[1]/S0+2*SX[1]^3/S0^3
              EUU2[1,1,1] <- -SXXX[1,1,1]/S0+SX[1]/S0*SXX[1,1]/S0
            } else if(k+s == 2){
              EUUU[k,s,t] <- SXXX[k,s,t]/S0 - SXX[1,1]/S0*SX[t]/S0-2*SXX[1,t]/S0*SX[1]/S0+2*SX[1]^2*SX[t]/S0^3
              EUU2[k,s,t] <- 0
            } else if(k+t == 2){
              EUUU[k,s,t] <- SXXX[k,s,t]/S0 - SXX[1,1]/S0*SX[s]/S0-2*SXX[1,s]/S0*SX[1]/S0+2*SX[1]^2*SX[s]/S0^3
              EUU2[k,s,t] <- 0
            } else if(s+t == 2){
              EUUU[k,s,t] <- SXXX[k,s,t]/S0 - SXX[1,1]/S0*SX[k]/S0-2*SXX[1,k]/S0*SX[1]/S0+2*SX[1]^2*SX[k]/S0^3
              EUU2[k,s,t] <- -SXXX[k,1,1]/S0+SX[k]/S0*SXX[1,1]/S0
            } else if(k==1){
              EUUU[k,s,t] <- SXXX[k,s,t]/S0-SXX[1,s]/S0*SX[t]/S0-SXX[s,t]/S0*SX[1]/S0+2*SX[1]*SX[s]*SX[t]/S0^3-SXX[1,t]/S0*SX[s]/S0
              EUU2[k,s,t] <- 0
            } else if(s==1){
              EUUU[k,s,t] <- SXXX[k,s,t]/S0-SXX[1,k]/S0*SX[t]/S0-SXX[k,t]/S0*SX[1]/S0+2*SX[1]*SX[k]*SX[t]/S0^3-SXX[1,t]/S0*SX[k]/S0
              EUU2[k,s,t] <- 0
            } else if(t==1){
              EUUU[k,s,t] <- SXXX[k,s,t]/S0-SXX[1,k]/S0*SX[s]/S0-SXX[k,s]/S0*SX[1]/S0+2*SX[1]*SX[k]*SX[s]/S0^3-SXX[1,s]/S0*SX[k]/S0
              EUU2[k,s,t] <- 0
            } else{
              EUUU[k,s,t] <- SXXX[k,s,t]/S0-SX[s]/S0*SXX[k,t]/S0-SX[k]/S0*SXX[s,t]/S0-SX[t]/S0*SXX[k,s]/S0+2*SX[k]*SX[s]*SX[t]/S0^3
              EUU2[k,s,t] <- 0
            }
          }
        }
      }

      # Compute the score U and second derivative U2
      if(!(ccmethod=="meandose")){
        U <- c(X[longx[,names(data)[status]]==1 & longx[,names(data)[loc]] == longx$loc,1]-SX[1]/S0, X[longx[,names(data)[status]]==1 & longx[,names(data)[loc]] == longx$loc,-1]-SX[-1]/S0)
        tmp1 <- -X[longx[,names(data)[status]]==1 & longx[,names(data)[loc]] == longx$loc,1]^2+SX[1]^2/S0^2
      } else {
        U <- c(X[longx[,names(data)[status]]==1,1]-SX[1]/S0, X[longx[,names(data)[status]]==1 ,-1]-SX[-1]/S0)
        tmp1 <- -X[longx[,names(data)[status]]==1,1]^2+SX[1]^2/S0^2
      }


      tmp2 <- SX[1]*SX[-1]/S0^2-SXX[1,-1]/S0
      tmp3 <- SXSX[-1,-1]/S0^2-SXX[-1,-1]/S0

      U2 <- rbind(c(tmp1, tmp2),cbind(tmp2,tmp3))

      list(S0=S0, SX=SX, SXX=SXX, SXXX=SXXX, SXSX=SXSX, EUU=EUU, EUUU=EUUU, EUU2=EUU2, U=U, U2=U2)
    })
  } else {
    # Computation for the reparametrized setting structured equivalently to the original parametrization. See documentation above for more details
    params[1] <- exp(params[1])
    byset <- lapply(sp, function(x){
      x <- as.data.frame(x)
      longx <- reshape(x, direction="long", varying=list(doses), v.names="dose", timevar="loc")

      if(ccmethod=="CCAL"){
        X <- as.matrix(cbind(D=longx$dose/(1+params[1]*longx$dose), model.matrix(~factor(loc)-1, data=longx)[,-1], longx[, names(x)[corrvars], drop=FALSE]))

      } else if(ccmethod=="CL"){
        longx <- longx[longx[,names(data)[status]]==1,]
        X <- as.matrix(cbind(D=longx$dose/(1+params[1]*longx$dose), model.matrix(~factor(loc)-1, data=longx)[,-1]))

      } else if(ccmethod%in%c("CCML")){
        longx <- longx[longx$loc==longx[,names(data)[loc]],]
        X <- as.matrix(cbind(D=longx$dose/(1+params[1]*longx$dose), longx[, names(x)[corrvars], drop=FALSE]))

      } else if(ccmethod=="meandose"){
        longx <- cbind(x[, set, drop = FALSE], data.frame(dose = rowMeans(x[, doses])), x[, -c(set, doses), drop = FALSE])
        X <- as.matrix(cbind(D = rowMeans(x[, doses])/(1 + params[1] * rowMeans(x[, doses])), x[, names(x)[corrvars], drop = FALSE]))
      }


      p <- (1+params[1]*longx$dose)*exp(c(X[,-1,drop=FALSE]%*%params[-1]))/sum((1+params[1]*longx$dose)*exp(c(X[,-1,drop=FALSE]%*%params[-1])))

      S0 <- sum((1+params[1]*longx$dose)*exp(c(X[,-1,drop=FALSE]%*%params[-1])))
      SX <- ((1+params[1]*longx$dose)*exp(c(X[,-1,drop=FALSE]%*%params[-1])))%*%X

      SXX <- Reduce('+',lapply(1:nrow(longx), function(j) X[j,]%*%t(X[j,])*(1+params[1]*longx$dose[j])*exp(c(X[,-1,drop=FALSE]%*%params[-1]))[j]))

      SXSX <- t(SX)%*%SX

      SXXX <- array(dim=c(ncol(X),ncol(X),ncol(X)))
      for(k in 1:ncol(X)){
        for(s in 1:ncol(X)){
          for(t in 1:ncol(X)){
            SXXX[k,s,t] <- sum(X[,k]*X[,s]*X[,t]*(1+params[1]*longx$dose)*exp(c(X[,-1,drop=FALSE]%*%params[-1])))
          }
        }
      }

      EX <- p%*%X
      EX2 <- p%*%X^2

      tmp1 <- params[1]^2*(EX2[1]-EX[1]^2)
      tmp2 <- params[1]*(SXX[1,-1]/S0-SX[1]*SX[-1]/S0^2)
      tmp3 <- SXX[-1,-1]/S0- SXSX[-1,-1]/S0^2

      EUU <- rbind(c(tmp1, tmp2),cbind(tmp2,tmp3))



      EUUU <- array(dim=c(ncol(X),ncol(X),ncol(X)))
      EUU2 <- array(dim=c(ncol(X),ncol(X),ncol(X)))
      for(k in 1:ncol(X)){
        for(s in 1:ncol(X)){
          for(t in 1:ncol(X)){
            if(k+s+t == 3){
              EUUU[1,1,1] <- params[1]^3*(SXXX[1,1,1]/S0-3*(SXX[1,1]/S0)*SX[1]/S0+2*SX[1]^3/S0^3)
              EUU2[1,1,1] <- params[1]^3*(-SXXX[1,1,1]/S0+SX[1]/S0*SXX[1,1]/S0)+params[1]^2*(SXX[1,1]/S0-SX[1]^2/S0^2)
            } else if(k+s == 2){
              EUUU[k,s,t] <- params[1]^2*(SXXX[k,s,t]/S0 - SXX[1,1]/S0*SX[t]/S0-2*SXX[1,t]/S0*SX[1]/S0+2*SX[1]^2*SX[t]/S0^3)
              EUU2[k,s,t] <- 0
            } else if(k+t == 2){
              EUUU[k,s,t] <- params[1]^2*(SXXX[k,s,t]/S0 - SXX[1,1]/S0*SX[s]/S0-2*SXX[1,s]/S0*SX[1]/S0+2*SX[1]^2*SX[s]/S0^3)
              EUU2[k,s,t] <- 0
            } else if(s+t == 2){
              EUUU[k,s,t] <- params[1]^2*(SXXX[k,s,t]/S0 - SXX[1,1]/S0*SX[k]/S0-2*SXX[1,k]/S0*SX[1]/S0+2*SX[1]^2*SX[k]/S0^3)
              EUU2[k,s,t] <- params[1]^2*(-SXXX[k,s,t]/S0 + SXX[1,1]/S0*SX[k]/S0)+params[1]*(SXX[1,k]/S0-SX[1]*SX[k]/S0^2)
            } else if(k==1){
              EUUU[k,s,t] <- params[1]*(SXXX[k,s,t]/S0-SXX[1,s]/S0*SX[t]/S0-SXX[s,t]/S0*SX[1]/S0+2*SX[1]*SX[s]*SX[t]/S0^3-SXX[1,t]/S0*SX[s]/S0)
              EUU2[k,s,t] <- 0
            } else if(s==1){
              EUUU[k,s,t] <- params[1]*(SXXX[k,s,t]/S0-SXX[1,k]/S0*SX[t]/S0-SXX[k,t]/S0*SX[1]/S0+2*SX[1]*SX[k]*SX[t]/S0^3-SXX[1,t]/S0*SX[k]/S0)
              EUU2[k,s,t] <- 0
            } else if(t==1){
              EUUU[k,s,t] <- params[1]*(SXXX[k,s,t]/S0-SXX[1,k]/S0*SX[s]/S0-SXX[k,s]/S0*SX[1]/S0+2*SX[1]*SX[k]*SX[s]/S0^3-SXX[1,s]/S0*SX[k]/S0)
              EUU2[k,s,t] <- 0
            } else{
              EUUU[k,s,t] <- SXXX[k,s,t]/S0-SX[s]/S0*SXX[k,t]/S0-SX[k]/S0*SXX[s,t]/S0-SX[t]/S0*SXX[k,s]/S0+2*SX[k]*SX[s]*SX[t]/S0^3
              EUU2[k,s,t] <- 0
            }
          }
        }
      }

      if(!(ccmethod=="meandose")){
        U <- c(params[1]*(X[longx[,names(data)[status]]==1 & longx[,names(data)[loc]] == longx$loc,1]-SX[1]/S0), X[longx[,names(data)[status]]==1 & longx[,names(data)[loc]] == longx$loc,-1]-SX[-1]/S0)
        tmp1 <- params[1]^2*(-X[longx[,names(data)[status]]==1 & longx[,names(data)[loc]] == longx$loc,1]^2+SX[1]^2/S0^2)+params[1]*(X[longx[,names(data)[status]]==1 & longx[,names(data)[loc]] == longx$loc,1]-SX[1]/S0)

      } else {
        U <- c(params[1]*(X[longx[,names(data)[status]]==1,1]-SX[1]/S0), X[longx[,names(data)[status]]==1,-1]-SX[-1]/S0)
        tmp1 <- params[1]^2*(-X[longx[,names(data)[status]]==1,1]^2+SX[1]^2/S0^2)+params[1]*(X[longx[,names(data)[status]]==1,1]-SX[1]/S0)

      }

      tmp2 <- params[1]*(SX[1]*SX[-1]/S0^2-SXX[1,-1]/S0)
      tmp3 <- SXSX[-1,-1]/S0^2-SXX[-1,-1]/S0

      U2 <- rbind(c(tmp1, tmp2),cbind(tmp2,tmp3))

      list(S0=S0, SX=SX, SXX=SXX, SXXX=SXXX, SXSX=SXSX, EUU=EUU, EUUU=EUUU, EUU2=EUU2, U=U, U2=U2)
    })
  }

  # Sum all components
  EUU <- Reduce('+',lapply(byset, function(x) x$EUU))
  EUUU <- Reduce('+',lapply(byset, function(x) x$EUUU))
  EUU2 <- Reduce('+',lapply(byset, function(x) x$EUU2))

  U <- Reduce('+',lapply(byset, function(x) x$U))
  U2 <- Reduce('+',lapply(byset, function(x) x$U2))

  infomat <- solve(EUU)

  # A will contain the Firth modification terms
  A <- rep(0, length(params))

  # Compute modifications using equation (B9) in manuscript
  for(r in 1:length(params)){
    res <- 0
    for(s in 1:length(params)){
      for(t in 1:length(params)){
        for(u in 1:length(params)){
          for(v in 1:length(params)){
            res <- res + U2[r,s]*infomat[s,t]*infomat[u,v]*(EUUU[t,u,v]+EUU2[t,u,v])
          }
        }
      }
    }
    A[r] <- -res/2

  }

  if(ccmethod == "CCAL"){
    names(U) <- names(A) <- c(ifelse(repar,"xi","beta"), paste0("alpha",2:length(doses)),names(data)[corrvars])
  } else if (ccmethod=="CL"){
    names(U) <- names(A) <- c(ifelse(repar,"xi","beta"), paste0("alpha",2:length(doses)))
  } else {
    names(U) <- names(A) <- c(ifelse(repar,"xi","beta"),names(data)[corrvars])
  }

  list(U=U, A=A)
}

