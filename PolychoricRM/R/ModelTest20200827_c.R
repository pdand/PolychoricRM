### It contains two functions compute.stat and Model.Fit
### The file was made on 2017-08-08

################################################################

### The function Compute.stat implments Proposition 4 in Browne 1984.

Compute.stat <- function(R0,u.r,Unrotated, acm.type) {
  
  # Input variables: 
  # R0 -> the sample correlation matrix
  # u.r -> the asymptotic covariance matrix of sample correlations
  # Unrotated -> the unrotated factor loading matrix
  # acm.type -> 1: u.r = Yhat; 2: u.r is of p2 by p2  
  
  # Output variables: the test statistic and the degrees of freedom
  # statistic 
  # df
  
  
  # Step 0, housekeeping
  #library(MASS)
  p = nrow(Unrotated)
  m = ncol(Unrotated)
  p.star = p * (p-1)/2
  p2 = p * p
  
  # Step 0.1, compute df
  df = ((p-m)**2 - p - m ) /2
  
if (df > 0) {
    
  # Step 1, make a selection matrix
#  M.Select = matrix(0,p2,p.star)
#  ij=0
#  ij.new = 0
#  for (j in 1:p) {
#    for (i in 1:p) {
#      ij = ij + 1
#      if (i<j) {
#        ij.new = ij.new + 1
#        M.Select[ij, ij.new] = 1
#      }
#    } # i
#  } # j
  
  
  # Step 2, compute the Delta matrix and its null matrix
  
  Delta.3d = array( rep(0, p*p*( p*m)), dim=c(p,p,( p*m ) ))
  
  
  ij = 0
  
  for (j in 1:m) {
    for (i in 1:p) {
      ij = ij + 1
      Delta.3d[i,1:p,ij] = Unrotated[1:p,j]
      Delta.3d[1:p,i,ij] = Unrotated[1:p,j] 
    } # i
  } # j
  
  Delta = matrix(0,p.star,(p*m))
  
  ij.new=0
  for (j in 2:p) {
    for (i in 1:(j-1)) {
      ij.new = ij.new + 1
      Delta[ij.new,1:(p*m)] = Delta.3d[i,j,1:(p*m)]
    } # i 
  } # j
  
  Delta.c = Null (Delta)
  
  
  # Step 3, Compute the test statistics
  ### Step 3.1, select non-duplicated elements from the asymptotic covariance matrix
   if (acm.type == 1) {
     
     Y.Hat = u.r
     
   } else {
     
     Y.Hat = matrix(0, p.star, p.star)
     u.r.col = matrix(0, p2, p.star)
     
     ij = 0
     for (j in 2:p) {
       for (i in 1:(j-1)) {
         ij = ij + 1
         u.r.col[,ij] = u.r[,((j-1)*p + i)]
       }
     }
     
     
     ij = 0
     for (j in 2:p) {
       for (i in 1:(j-1)) {
         ij = ij + 1
         Y.Hat[ij,] = u.r.col[((j-1)*p + i) , ]
       }
     }
     
     
   } 
   
  ### Step 3.2, compute the residuals and select the nonduplicated elements
  Residual = R0 - Unrotated %*% t(Unrotated)
  
  r.v = rep(0,p.star)
  
  ij = 0
  for (j in 2:p) {
    for (i in 1:(j-1)) {
      ij = ij + 1
      r.v[ij] = Residual[i,j]
    }
  }
  
  
  
  ### Step 3,3, compute the sample statistic
  temp.Delta.c.e = t(Delta.c) %*% r.v

  temp.middle = t(Delta.c) %*% Y.Hat %*% Delta.c
  
 # ev = eigen(Y.Hat)$values
  ev = eigen(temp.middle,symmetric=TRUE,only.values=TRUE)$values
  
  if (min(ev)<= 1.0e-6) {
   statistic = NaN 
  } else {
  temp.right = solve(temp.middle, temp.Delta.c.e)
  statistic = sum(temp.Delta.c.e * temp.right)
  }
  
} else {
  statistic = 0
}  
  
  
  list (statistic = statistic, df = df)
  
}


##------------------------------------------------------------------------------------------



Model.Fit <- function(statistic.sample,fm,p,m,n,confid.level) {

## Model.Fit computes p values for close fit, perfect fit, RMSEA and its confidence intervals,
## ECVI and its confidence intervals for ML  
    
  ## Housekeeping
  pstar = p * (p-1)/2
  q = p * m - m * (m-1) /2
  df = pstar - q
  q = q + p # unique variances are also paramaters
  
  ## Compute an estimate for ECVI
  ECVI = ( statistic.sample + 2 * q ) / n 
  
  
  ### Compute less biased estimate for the RMSEA
  if (statistic.sample < df) {
    
    RMSEA = 0
    
  } else {
    
    RMSEA = sqrt((statistic.sample - df)/(df*n))
    
  }
  
  ### compute the p values for the test of perfect fit and close fit
  
  p.perfect = 1 - pchisq(statistic.sample,df,ncp=0)
  p.close = 1 - pchisq(statistic.sample,df,ncp=0.0025*n*df)
  
  
  ### Compute the lower end and the upper end of the noncentrality parameter
  
  find.lambda <- function(x,f.dis,d,prob) pchisq(f.dis,d,ncp=x) - prob
  
  lambda.range = statistic.sample
  for (i in 1:10) {
    if (pchisq(statistic.sample,df,ncp=lambda.range)<(1-confid.level)/2) break
    lambda.range = lambda.range * 2
  }
  
  
  
  if (p.perfect > (1 - confid.level)/2) {
    l.lambda = 0
  } else {
    temp = uniroot (find.lambda, interval= c(0,lambda.range),f.dis = statistic.sample, d=df, prob= (0.5 + confid.level/2) )
    l.lambda = temp$root 
  }
  
  if (p.perfect > (0.5 + confid.level)/2) {
    u.lambda = 0
  } else{
    temp = uniroot (find.lambda, c(0,lambda.range),f.dis = statistic.sample, d=df, prob= (0.5 - confid.level/2) )
    u.lambda = temp$root 
  }  
  
  
  ### 
  
  if (fm=='ml') {
    n.ECVI = n - p - 1
  } else {n.ECVI = n} 
  
  
  list (f.stat = statistic.sample, df = df, n = n, RMSEA = RMSEA, p.perfect=p.perfect, p.close = p.close,
        confid.level = confid.level,RMSEA.l = sqrt(l.lambda/(n*df)), RMSEA.u = sqrt(u.lambda/(n*df)), 
        ECVI = ECVI, ECVI.l=(l.lambda + pstar + q)/n.ECVI, ECVI.u=(u.lambda + pstar + q)/n.ECVI)
  
  
} # Model.Fit

#---------------------------------------------------------------------------------------------------
