
# The function BootJack compute SE estimates using the bootstrap or jackknife method.



BootJack <- function(X.raw, M.order, BJ.Arg) {
  
  # Two internal functions Get.MR and fact.extract
  # It requires two external functions and a R package.
  
  # source('Alignment20160616.R')
  # source('EFAEstimation20160811.R')
  # library(GPArotation)
  
  
  #----------------------------------------------------------------------------
  
  Get.MR <- function(X.raw, bj='bootstrap', Ib = NULL, dist = 'continuous') {
    
    # Get.MR stands for get manifest variable correlation matrices
    
    if ((bj=='bootstrap') & (is.null(Ib))) Ib = 2000
    
    p = ncol(X.raw)
    n = nrow(X.raw)
    
    if (bj =='jackknife') Ib = n 
    
    Total.MR = array(rep(0, Ib * p * p), dim=c(Ib,p,p))
    
    
    for (i in 1:Ib) {
      
      if (bj=='bootstrap') {
        
        Temp.raw = matrix(0,n,p)
        
        boot.index = sample(1:n, n, replace = TRUE)  
        
        for (j in 1:n) {
          Temp.raw[j,1:p] = X.raw[boot.index[j],1:p]
        }
        
      } else {
        
        Temp.raw = matrix(0,(n-1),p)
        
        if (i==1) {
          Temp.raw = X.raw[2:n, 1:p]
        } else if (i==n) {
          Temp.raw = X.raw[1:(n-1),1:p]  
        } else {
          Temp.raw[1:(i-1),1:p] = X.raw[1:(i-1),1:p]
          Temp.raw[i:(n-1),1:p] = X.raw[(i+1):n, 1:p]
        } #  1 < i < n
        
        
      } # jackknife
      
      if (dist=='ordinal') {
        Total.MR[i,1:p, 1:p] = PolychoricRM(Temp.raw)$correlation
      } else {
        Total.MR[i,1:p, 1:p] = cor(Temp.raw)
      }
      
      
    } # for (i in 1:Ib)
    
    Total.MR
    
  } # Get.MR <- function
  
  #------------------------------------------------------------------------
  
  
  # Dummy arguments
  # bj -> 'bootstrap' | 'jackknife'
  # Ib -> # of bootstrap samples
  # rypte -> 'oblique' | 'orthogonal'
  # Level.Confid -> 
  
  # EF.Arg -> dummy arguments for factor extraction
  # Rotation.Arg -> dummy arguments for factor rotation
  
  
  bj = BJ.Arg$bj
  rtype = BJ.Arg$rtype
  Ib = BJ.Arg$Ib
  Level.Confid = BJ.Arg$Level.Confid
  FE.Arg = BJ.Arg$FE.Arg
  Rotation.Arg = BJ.Arg$Rotation.Arg
  fnames = BJ.Arg$fnames
  dist = BJ.Arg$dist
  
  n = nrow(X.raw)
  p = ncol(X.raw)
  m = ncol(M.order)
  
  
  ### Step 1: get bootstrap / jackknife manifest variable correlation matrices ###
  
  if (bj == 'bootstrap') {
    Boot.R = Get.MR(X.raw,bj='bootstrap',Ib = Ib) ## 
  } else {
    Jack.R = Get.MR(X.raw,bj='jackknife') ## 
  }
  
  ### Step 2: extract factors   ###
  
  Total.A <- array(rep(0,Ib*p*m),dim=c(Ib,p,m))
  Total.Index <- matrix(0,Ib,4)
  
  
  for (i in 1:Ib) {
    
    if (bj == 'bootstrap') {
      Test = do.call('fa.extract',append(list(Boot.R[i,,]),FE.Arg))
    } else {
      Test = do.call('fa.extract',append(list(Jack.R[i,,]),FE.Arg))
    }
    
    Total.A[i,1:p,1:m] = Test$Unrotated
    Total.Index[i,1] = Test$f[1]
    Total.Index[i,2] = Test$convergence
    Total.Index[i,3] = Test$heywood
    
  }
  
  
  ### Step 3: conduct factor rotations   ###
  
  
  Total.LPhi <- array(rep(0,Ib*(p+m)*m),dim=c(Ib,(p+m),m))
  Total.Phi <- array(rep(0,Ib*m*m),dim=c(Ib,m,m))
  
  
  
  for (i in 1:Ib) {
    
    Test = do.call(fnames,append(list(Total.A[i,1:p,1:m]),Rotation.Arg)) 
    Total.LPhi[i,1:p,1:m] = Test$loadings[1:p,1:m]
    Total.LPhi[i,(p+1):(p+m),1:m] = diag(m)
    if (rtype == 'oblique') Total.LPhi[i,(p+1):(p+m),1:m] = Test$Phi[1:m,1:m]           
    if (! Test$convergence) Total.Index[i,4] = 1 
  }
  
  
  ### Step 4, Reflected the rotated factor loading matrices and factor correlation matrices   ###
  
  
  Total.Rotated = array(rep(0,Ib*(p+m+1)*m),dim=c(Ib,(p+m+1),m))
  
  for (i in 1:Ib) {
    Total.Rotated[i,1:(p+m+1),1:m] = Align.Matrix(M.order,Total.LPhi[i,1:(p+m),1:m])
  }
  
  
  ### Step 5, Compute SE and percentile confidence intervals  ###
  
  if (bj=='bootstrap') {
    scale = 1
  } else {
    scale = (n-1) / sqrt(n)
  }
  
  SE.LPhi <- matrix(0, p+m, m)
  
  for (j in 1:m) {
    for (i in 1:(p+m)) {
      SE.LPhi[i,j] = sd(Total.Rotated[1:Ib,i,j]) * scale
    }
  }
  
  
  
  
  
  if (bj=='bootstrap') {
    
    percentile = rep(0,2)
    percentile[1] = ( 1 - Level.Confid ) / 2
    percentile[2] = 1 - percentile[1]
    
    Confid.Boot <- array(rep(0, (p+m)*m*2), dim=c(p+m, m, 2))
    
    for (j in 1:m) {
      for (i in 1:(p+m)) {
        Confid.Boot[i,j,1:2] = quantile(Total.Rotated[1:Ib,i,j], percentile)
      }
    }
    
    
  } # (bj=='bootstrap')
  
  
  ###  Make a list of outputs   ###
  
  if (bj=='bootstrap') {
    
    list(SE.LPhi = SE.LPhi, Total.Index = Total.Index, Confid.Boot = Confid.Boot)
    
  } else { # jackknife
    
    list(SE.LPhi = SE.LPhi, Total.Index = Total.Index)
    
  } # jackknife
  
  
  
} # BootJack





