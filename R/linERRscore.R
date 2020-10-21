#' @title Compute the Firth-corrected score function
#' @description Compute the Firth-corrected score
#' @param params vector of parameters (beta, alpha, delta) for which to compute the modified score
#' @param data data frame containing matched case-control data, with a number of columns for doses to different locations, a column containing matched set numbers, a column containing the case's tumor location and a column serving as a case-control indicator.
#' @param doses vector containing the indices of columns containing dose information, in the desired order.
#' @param set column index containing matched set numbers.
#' @param status column index containing case status.
#' @param loc column index containing the location of the matched set's case's second tumor.
#' @param corrvars vector containing the indices of columns containing variables to be corrected for.
#' @param ccmethod choice of method of analysis; one of meandose, CCML, CCAL or CL.
#' @param repar reparametrize beta = exp(gamma)
#' @export



linERRscore <- function(params,data, doses, set, status, loc, ccmethod, corrvars=NULL, repar=FALSE){

  sp <- split(data,data[,set])
  sp <- lapply(sp, as.matrix)

  if(repar==FALSE){
    byset <- lapply(sp, function(x){
      x <- as.data.frame(x)
      longx <- reshape(x, direction="long", varying=list(doses), v.names="dose", timevar="loc")

      if(ccmethod=="CCAL"){
        X <- as.matrix(cbind(D=longx$dose/(1+params[1]*longx$dose), model.matrix(~factor(loc)-1, data=longx)[,-1], longx[, names(x)[corrvars], drop=FALSE]))

      } else if(ccmethod=="CL"){
        longx <- longx[longx[,names(data)[status]]==1,]
        X <- as.matrix(cbind(D=longx$dose/(1+params[1]*longx$dose), model.matrix(~factor(loc)-1, data=longx)[,-1], longx[, names(x)[corrvars], drop=FALSE]))

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

      tmp1 <- EX2[1]-EX[1]^2
      tmp2 <- SXX[1,-1]/S0-SX[1]*SX[-1]/S0^2
      tmp3 <- SXX[-1,-1]/S0- SXSX[-1,-1]/S0^2

      EUU <- rbind(c(tmp1, tmp2),cbind(tmp2,tmp3))



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
    params[1] <- exp(params[1])
    byset <- lapply(sp, function(x){
      x <- as.data.frame(x)
      longx <- reshape(x, direction="long", varying=list(doses), v.names="dose", timevar="loc")

      if(ccmethod=="CCAL"){
        X <- as.matrix(cbind(D=longx$dose/(1+params[1]*longx$dose), model.matrix(~factor(loc)-1, data=longx)[,-1], longx[, names(x)[corrvars], drop=FALSE]))

      } else if(ccmethod=="CL"){
        longx <- longx[longx[,names(data)[status]]==1,]
        X <- as.matrix(cbind(D=longx$dose/(1+params[1]*longx$dose), model.matrix(~factor(loc)-1, data=longx)[,-1], longx[, names(x)[corrvars], drop=FALSE]))

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

  EUU <- Reduce('+',lapply(byset, function(x) x$EUU))
  EUUU <- Reduce('+',lapply(byset, function(x) x$EUUU))
  EUU2 <- Reduce('+',lapply(byset, function(x) x$EUU2))

  U <- Reduce('+',lapply(byset, function(x) x$U))
  U2 <- Reduce('+',lapply(byset, function(x) x$U2))

  infomat <- solve(EUU)
  A <- rep(0, length(params))


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



  # if(ccmethod %in% c("CCML","meandose")){
  #   U <- c(U[1], rep(0, length(doses)-1),U[-1])
  #   A <- c(A[1], rep(0, length(doses)-1),A[-1])
  # }
  if(ccmethod %in% c("CCAL","CL")){
    names(U) <- names(A) <- c(ifelse(repar,"xi","beta"), paste0("alpha",2:length(doses)),names(data)[corrvars])
  } else {
    names(U) <- names(A) <- c(ifelse(repar,"xi","beta"),names(data)[corrvars])
  }

  list(U=U, A=A)
}

