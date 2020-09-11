AsyCovCorr <- function(X) {

## The function AsyCovCorr implements the asymptotic covariance matrix of correlations for an arbitrary distribution
## described in Browne & Shapiro (1986), the asymptotic covariance matrix of sample correlation coefficients under
## general conditions. Linear Algebra and its applications, 82, 169-176.
## Equation (3.2) on page 171.


p = ncol(X)
N = nrow(X)

Xmean = colMeans(X)

Xcenter = matrix(0,N,p)
for (j in 1:p) {
 Xcenter[,j] = X[,j] - Xmean[j]
}


S0 = crossprod(Xcenter) / N
SD = sqrt(diag(S0))

R0 = S0

for (j in 1:p) {
R0[,j] = R0[,j] / SD[j]
R0[j,] = R0[j,] / SD[j]
}


Theta = matrix(0, p*p, p*p)
ij = 0

for (j in 1:p) {
 for (i in 1:p) {
   ij = ij + 1
   kl = 0
    for (l in 1:p) {
      for (k in 1:p) {
      kl = kl + 1

      Theta[ij,kl] = mean(Xcenter[,i] * Xcenter[,j] * Xcenter[,k] * Xcenter[,l]) / (SD[i]*SD[j]*SD[k]*SD[l])

     } # k
   } # l
  } # i
} # j

Gamma = Theta - outer(array(R0),array(R0))  ## Equation (3.1) on page 171


A = matrix(0,p*p,p)  # Equation (3,3) on page 171

ij = 0
for (j in 1:p) {
  for (i in 1:p) {
    ij = ij + 1
     if (i==j) {  
     A[ij,i] = R0[i,j]
     } else {
     A[ij,i] = R0[i,j] / 2
     A[ij,j] = R0[i,j] /2
     } 
          
 } # i 
} # j


B = matrix(0,p*p,p) # Equation (3.4) on page 171

ij = 0
for (j in 1:p) {
  for (i in 1:p) {
    ij = ij + 1
    for (k in 1:p) {
     kk = (k-1) * p + k
     B[ij,k] = Gamma[ij,kk]
   } # k
          
 } # i 
} # j


G = matrix(0,p,p) # Equation (3.5) on page 171

for (j in 1:p) {
  for (i in 1:p) {
    ii = (i - 1) * p + i
    jj = (j - 1) * p + j
    G[i,j] = Gamma[ii,jj]
   }
}

Temp2 = tcrossprod(A,B)
Temp4 = tcrossprod(tcrossprod(A,G), A)
C.Gamma = Gamma -  Temp2 - t(Temp2) + Temp4 # Equation (3.2) on page 171

list(mean=Xmean,sd = SD, corr = R0, asc = C.Gamma)

} ## AsyCovCorr <- function(X)
