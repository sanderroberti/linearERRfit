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
#' @param repar reparametrize to \eqn{\beta=exp(\xi)}? Note that this only affects the Firth modified score, not the original score. Defaults to FALSE.
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
#' @importFrom data.table data.table ':=' .SD
#' @export



linERRscore <- function(params,data, doses, set, status, loc, ccmethod, corrvars=NULL, repar=FALSE){

  data <- data[data[,set] %in% as.numeric(names(which(table(data[,set],data[,status])[,2]==1))),] # take only the matched sets with one case patient
  if(ccmethod!="CL") data <- data[data[,set] %in% as.numeric(names(which(table(data[,set])>1))),] # for methods other than CL, take only matched sets with at least 2 patients


  if(ccmethod=="CL") corrvars <- NULL

  # alpha are location effects, which are not used for CCML and mean dose methods
  if(ccmethod %in% c("CCAL","CL")){
    alpha <- c(0,params[2:(length(doses))])
  } else {
    alpha <- rep(0,length(doses))
  }

  if(repar) params[1] <- exp(params[1])

  # gamma are other patient-level effect parameters
  gamma <- tail(params, length(corrvars))

  # reshape data frame to long format with different organ locations in different rows
  datalong <- reshape(data, direction="long", varying=list(doses), v.names="dose", timevar="loc")
  datalong <- datalong[order(datalong[,names(data)[set]]),] # order the long dataframe based on matched set nr

  # Construct the object X: a sort of 'design matrix' based on the long data frame. The first column is D=d/(1+beta*d), which is followed by dummy variables for all locations except the first, and finally by the columns indicated by corrvars in the original data frame
  if(ccmethod=="CCAL"){
    X <- data.frame(D=datalong$dose/(1+params[1]*datalong$dose), model.matrix(~factor(loc)-1, data=datalong)[,-1], datalong[, names(data)[corrvars], drop=FALSE])
  } else if(ccmethod=="CCML"){
    datalong <- datalong[datalong$loc == datalong[,names(data)[loc]],]
    X <- data.frame(D=datalong$dose/(1+params[1]*datalong$dose), datalong[, names(data)[corrvars], drop=FALSE])
  } else if(ccmethod=="CL"){
    datalong <- datalong[datalong[names(data)[status]]==1,]
    X <- data.frame(D=datalong$dose/(1+params[1]*datalong$dose), model.matrix(~factor(loc)-1, data=datalong)[,-1], datalong[, names(data)[corrvars], drop=FALSE])
  } else if(ccmethod=="meandose"){
    datalong <- data[,-c(doses)]
    datalong$dose <- rowMeans(data[,doses])
    X <- data.frame(D=datalong$dose/(1+params[1]*datalong$dose), datalong[, names(data)[corrvars], drop=FALSE])
  }

  # Construct a long version of S0 (see Appendix in the paper). The proper S0 for a matched set is obtained by summing S0long within the matched set
  # S0long is a vector
  S0long <- (1+params[1]*datalong$dose)*exp(as.matrix(datalong[,names(data)[corrvars], drop=FALSE])%*%gamma)
  if(ccmethod!="meandose") S0long <- S0long*exp(model.matrix(~factor(loc)-1, data=datalong)%*%alpha)

  # Create S0 by aggregating S0long within each matched set
  S0 <- as.data.frame(data.table(S0long)[, list(S0=sum(V1)), by=list(set=datalong[,names(data)[set]])]) #aggregate(list(S0=S0long),by=list(set=datalong[,names(data)[set]]), FUN=sum)

  #SX <- lapply(X, function(x) as.data.frame(data.table(x*S0long)[, list(SX=sum(V1)), by=list(set=datalong[,names(data)[set]])]))#aggregate(list(SX=x*S0long), by=list(set=datalong[,names(data)[set]]),FUN=sum))

  # Construct SX. The object is a matrix with columns S_d, S_X1, S_X2, ... Note that the location effect variables are part of the S_X. Each row corresponds with a matched set
  SX <- sapply(X, function(x) {
    tmp <- as.data.frame(data.table(x*S0long)[, list(SX=sum(V1)), by=list(set=datalong[,names(data)[set]])])
    tmp$SX[match(datalong[,names(data)[set]], tmp$set)]
  })


  # matlist <- split(X, datalong[,names(data)[set]])
  # matlist <- lapply(matlist, as.matrix)
  #
  # matlist2 <- split(S0long, datalong[,names(data)[set]])
  # matlist2 <- lapply(matlist2, function(x) matrix(rep(x, each=ncol(X)), ncol=ncol(X),byrow=TRUE))
  # diagmat2 <- bdiag(matlist2)
  #
  # diagmat <- bdiag(matlist)
  #
  # diagmat <- t(diagmat)%*%(diagmat*diagmat2)
  #
  #
  #
  # Avec <- c(1,rep(0, ncol(X)-1))
  # Amat0 <- bdiag(lapply(1:length(matlist), function(x) matrix(Avec, nrow=1)))

  # Construct SXX, which is a list object. Each list element is a matrix, for example the first list element is a matrix with columns S_dd, S_dX1, S_dX2, ...
  SXX <- lapply(1:ncol(X),  function(x) matrix(nrow=nrow(datalong),ncol=ncol(X)))
  #setnrlong <- data.table(set=datalong[,names(data)[set]])
  dt1 <- data.table(X*S0long)
  for(c1 in 1:ncol(X)){ # For each column in the matrix X, I construct a new list element

    #Amat <- cbind(matrix(0, nrow=nrow(Amat0), ncol=c1-1), Amat0[,1:(ncol(Amat0)-c1+1)])

    #for(c2 in c1:ncol(X)){
    #SXX[[c1]][[c2]] <- SXX[[c2]][[c1]] <- as.data.frame(data.table(X[,c1]*X[,c2]*S0long)[, list(SXX=sum(V1)), by=list(set=datalong[,names(data)[set]])])
    #SXX[[c1]][[c2]] <- SXX[[c2]][[c1]] <- data.frame(set=unique(datalong[,names(data)[set]]), SXX=diagmat[cbind(c1+(1:length(matlist)-1)*ncol(X),c2+(1:length(matlist)-1)*ncol(X))])

    #tmp <- as.data.frame(data.table(X[,c1]*X[,c2]*S0long)[, list(SXX=sum(V1)), by=list(set=datalong[,names(data)[set]])])
    #SXX[[c1]][,c2] <- SXX[[c2]][,c1] <- tmp$SXX[match(datalong[,names(data)[set]], tmp$set)]

    # create a data table by multiplying column c1 of X with S0long and with each column of X, and then aggregate within each matched set
    tmp <- data.table(set=datalong[,names(data)[set]],X[,c1]*dt1)[,(names(X)):=lapply(.SD, sum), by="set"]
    #
    #

    # Transform the data table to a matrix and store in the list
    SXX[[c1]] <- as.matrix(tmp[,-1])

    # Bvec <- rep(0, ncol(X))
    # Bvec[c2] <- 1
    # Bmat <- matrix(rep(Bvec, length(matlist)), ncol=1)
    #
    # #SXX[[c1]][[c2]] <- SXX[[c2]][[c1]] <- data.frame(set=unique(datalong[,names(data)[set]]), SXX= as.numeric(Amat%*%diagmat%*%Bmat))
    #
    # tmp <- data.frame(set=unique(datalong[,names(data)[set]]), SXX= as.numeric(Amat%*%diagmat%*%Bmat))
    #
    # SXX[[c1]][,c2] <- SXX[[c2]][,c1] <- tmp$SXX[match(datalong[,names(data)[set]], tmp$set)]
    #
    # #SXX[[c1]][[c2]] <- SXX[[c2]][[c1]] <- setnrlong[tmp, on='set']$SXX

    #}
  }
  # SXX <- lapply(X, function(x) {
  #   lapply(X, function(y){
  #     as.data.frame(data.table(x*y*S0long)[, list(SXX=sum(V1)), by=list(set=datalong[,names(data)[set]])])
  #     #data.frame(set=datalong[,names(data)[set]],SXX=x*y*S0long) %>% group_by(set) %>% summarise(SXX=sum(SXX), .groups="drop")
  #     #aggregate(list(SXX=x*y*S0long), by=list(set=datalong[,names(data)[set]]),FUN=sum)
  #   })
  # })

  # Take S0long and aggregate within matched sets to obtain the total relative in a matched set. Then use that as the denominator in the computation of conditional probabilities of each location being a case
  S0setlong <- S0$S0[match(datalong[,names(data)[set]], S0$set)]
  probs <- S0long/S0setlong # conditional tumor probabilities for each location, to be used for computing expected values

  # Compute a long version of U, which is a matrix. Each column corresponds to a derivative with respect to one variable. Each value in the matrix is the contribution to the derivative for the ENTIRE MATCHED SET, assuming the location corresponding to the row IS THE TUMOR LOCATION.
  # Later, the proper derivative is extracted from this object by taking the correct subset of rows (i.e., only the observed tumor locations). The reason for this construction is to facilitate easy computation of expected values needed for Firth's bias correction
  Ulong <- sapply(1:ncol(X), function(x) X[,x]-SX[,x]/S0setlong)
  if(repar) Ulong[,1] <- Ulong[,1]*params[1]

  # U2long <- lapply(1:ncol(X), function(x){
  #   sapply(1:ncol(X), function(y) (SX[,x]*SX[,y]-S0setlong*SXX[[x]][,y])/S0setlong^2)
  # })

  # Similar to Ulong, compute a long matrix of second derivatives. This is a list, similar to SXX, with each list element a matrix. Column j of matrix i is the second derivative of log(L) with respect to variables i and j. Again, each value in the matrix is the contribution to the particular second derivative for the ENTIRE MATCHED SET, assuming the location corresponding to the row IS THE TUMOR LOCATION.
  U2long <- lapply(1:ncol(X), function(x){
    (SX[,x]*SX-S0setlong*SXX[[x]])/S0setlong^2
    #sapply(1:ncol(X), function(y) (SX[,x]*SX[,y]-S0setlong*SXX[[x]][,y])/S0setlong^2)
  })


  if(ncol(X)>1){
    for(ucol in 2:ncol(X)){
      U2long[[1]][,ucol] <- ifelse(repar,params[1],1)*U2long[[1]][,ucol] # Needs to be multiplied with exp(xi) when repar=TRUE
      U2long[[ucol]][,1] <- ifelse(repar,params[1],1)*U2long[[ucol]][,1] # Needs to be multiplied with exp(xi) when repar=TRUE
    }
  }
  U2long[[1]][,1] <- ifelse(repar,params[1],0)*(X[,1]-SX[,1]/S0setlong)+ifelse(repar, params[1]^2,1)*(SX[,1]^2/S0setlong^2-X[,1]^2)

  # Now extract the actual contributions to the derivatives U and U2 by taking only tumor location rows
  if(ccmethod!="meandose"){
    #U <- lapply(Ulong, function(x) x[datalong$loc==datalong[,names(data)[loc]] & datalong[,names(data)[status]]==1])
    U <- Ulong[datalong$loc==datalong[,names(data)[loc]] & datalong[,names(data)[status]]==1,,drop=FALSE]
    U2 <- lapply(U2long, function(x) {
      x[datalong$loc==datalong[,names(data)[loc]] & datalong[,names(data)[status]]==1,,drop=FALSE]}
    )
  } else { # For mean dose, take rows corresponding to case patients
    #U <- lapply(Ulong, function(x) x[datalong[,names(data)[status]]==1])
    U <- Ulong[datalong[,names(data)[status]]==1,,drop=FALSE]
    U2 <- lapply(U2long, function(x) {
      x[datalong[,names(data)[status]]==1,,drop=FALSE]}
    )
  }

  # Compute the expected values of second derivatives by multiplying each column of U2long with the conditional probabilities and summing
  EU2 <- sapply(U2long, function(x){
    apply(x,2,function(i) sum(i*probs))
  })

  # Compute the expected values of the product of derivatives using matrix multiplication. The object EUU is a matrix with (i,j)-element the expected second derivative of log(L) wrt variables i and j
  Umat <- Ulong
  Umat2 <- Umat*matrix(rep(probs,each=ncol(X)), ncol=ncol(X),byrow=TRUE)
  EUU <- t(Umat) %*% Umat2



  # EUU2 <- lapply(Ulong, function(x){
  #   sapply(U2long, function(y){
  #     unname(sapply(as.data.frame(y), function(z) sum(x*z*probs)))
  #   })
  # })
  #probslong <- matrix(rep(probs,each=ncol(X)),ncol=ncol(X),byrow=TRUE)

  # EUU2 <- lapply(1:ncol(X), function(x){
  #   sapply(U2long, function(y){
  #     t(Ulong[,x]*probs)%*%(y)
  #   })
  # })

  # Similarly, compute the expected values of all products of a first and second derivative. The output is a list with matrices. Element (j,k) of matrix i is the expected value of U_i*U_jk
  EUU2 <- lapply(U2long, function(x){
    t(Ulong*as.numeric(probs))%*%x
  })

  EUU2 <- lapply(1:ncol(X), function(x){
    sapply(EUU2, function(y) y[x,])
  })



  # Similarly, compute expected values of products of three first derivatives. The output has the same form as EUU2
  EUUU <- lapply(1:ncol(X), function(x){
    t(Umat) %*% (Umat2*Ulong[,x])
  })

  # Now obtain the first and second derivative of log(L) by summing all contributions over matched sets
  U <- colSums(U)
  U2 <- matrix(sapply(U2, function(x) colSums(x)), ncol=length(params))

  # EUU is the information matrix, we also need its inverse
  invinfomat <- solve(EUU)

  # A will contain the Firth modification terms
  A <- rep(0, length(params))

  # Compute modifications using equation (B5) in paper
  for(r in 1:length(params)){
    res <- 0
    for(s in 1:length(params)){
      for(t in 1:length(params)){
        res <- res + U2[r,s]*invinfomat[s,t]*sum(invinfomat*(EUUU[[t]]+EUU2[[t]]))
        #res <- res + U2[r,s]*invinfomat[s,t]*sum(invinfomat*(EUUU[[t]]+sapply(EUU2, function(x) x[,t])))
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

  list(U=U, A=A, infomat=EUU)

}

