# This independent function compares multiple solutions
# August 8, 2017

 

CompareSolutions <- function(lambdas, compare = 'First'){
  nsol <- length(lambdas)
  p <- dim(lambdas[[1]])[1]
  m <- dim(lambdas[[1]])[2]

  # Align sample solutions ==============
  LA <- list()
  LA[[1]] <- lambdas[[1]]
  for(j in 2:nsol) LA[[j]] <- Align.Matrix(LA[[1]], rbind(lambdas[[j]],diag(m)))[1:p,]
  
  Cong <- rep(NA, m)
  CMax <- CMean <- CMin <- matrix(0, nsol, nsol)
  CRaw <- matrix(0, nsol, m*nsol)
  for(i in 1:nsol){
    if(i < nsol){
      for(j in (i+1):nsol){
        for(k in 1:m){
  Cong[k] <- sum(LA[[i]][,k]*LA[[j]][,k])/(sqrt(sum(LA[[i]][,k]^2))*sqrt(sum(LA[[j]][,k]^2)))
        }
        CMin[i,j] <- CMin[j,i] <- min(abs(Cong))
        CMean[i,j] <- CMean[j,i] <- mean(abs(Cong))
        CMax[i,j] <- CMax[j,i] <- max(abs(Cong))
        CRaw[j,(i*m-(m-1)):(i*m)] <- abs(Cong)
      }
    }
  }
  diag(CMax) <- diag(CMean) <- diag(CMin) <- 1
  
  # Clean the output =================================
  sname <- c(paste('Solution',1:nsol,sep=''))
  rownames(CMean) <- colnames(CMean) <- rownames(CMin) <- colnames(CMin) <- sname
  rownames(CRaw) <- rownames(CMax) <- colnames(CMax) <- sname
  
  name <- rep(NA, m*nsol)
  for(i in 1:nsol) {
    for(k in 1:m) name[(i-1)*m+k] <- paste(sname[i], paste('F',k,sep = ""), sep="-")
  }
  colnames(CRaw) <- name
  
  # Output ===========================================
  if(compare == 'First'){
    list(MinimumCongruence = round(CMin[1,-1],3),
         MeanCongruence = round(CMean[1,-1],3),
         MaximumCongruence = round(CMax[1,-1],3),
         RawCongruences = round(CRaw[-1,1:m],3))
  } else if(compare == 'All'){
    list(MinimumCongruence = round(CMin,3),
         MeanCongruence = round(CMean,3),
         MaximumCongruence = round(CMax,3),
         RawCongruences = round(CRaw,3))
  }
} # ComapareSolutions() -------------------------------
