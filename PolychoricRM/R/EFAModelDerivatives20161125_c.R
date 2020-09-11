### The file contains five functions.
### EFA.Hessian is the head function


###########################################################################

EFA.Hessian <- function(Lambda, Phi, RSample, extraction=NULL) {

### SHessianML and SHessianOLS

#..................................................................................

# SHessianML computes the hessian matrix of the ML discrepancy. 
# Lambda <- (p,m)
# Phi <- (m,m), a symmetric matrix
# R.Sample <- (p,p), a symmetric matrix
# Result <- (ntemp,ntemp), a matrix of second order derivatives

SHessianML <- function(Lambda, Phi, R.Sample) {

# It invokes two external functions: DifS2LPhiPsi and Dif2S2LPhiPsi.

# library(MASS)

p = dim(Lambda)[1]
m = dim(Lambda)[2]
ntemp = p*m + m*(m-1)/2 + p

Result = array( rep(0,ntemp*ntemp), dim=c(ntemp, ntemp))

P = Lambda %*% Phi %*% t(Lambda)
for (i in 1:p) {
P[i,i] = 1
}

PInverse = solve(P)

MTemp1 = PInverse %*% (R.Sample - P) %*% PInverse

FirstDerivative = DifS2LPhiPsi (Lambda, Phi)
SecondDerivative = Dif2S2LPhiPsi (Lambda, Phi)


for (j in 1:ntemp) {
  for (i in 1:ntemp) {
    MTemp2 = FirstDerivative[1:p,1:p,i] %*% PInverse %*% FirstDerivative[1:p,1:p,j] 
    Result[i,j] = 2 * sum(MTemp1 * MTemp2) + sum (MTemp2 * PInverse) - sum ( MTemp1 * SecondDerivative[1:p,1:p,i,j] )        
  } # i
} # j

Result

} # SHessianML


#..................................................................................

# SHessianOLS computes the hessian matrix of the OLS discrepancy. 
# Lambda <- (p,m)
# Phi <- (m,m), a symmetric matrix
# R.Sample <- (p,p), a symmetric matrix
# Result <- (ntemp,ntemp), a matrix of second order derivatives

SHessianOLS <- function(Lambda, Phi, R.Sample) {

# It invokes two external functions: DifS2LPhiPsi and Dif2S2LPhiPsi.


# library(MASS)

p = dim(Lambda)[1]
m = dim(Lambda)[2]
ntemp = p*m + m*(m-1)/2 + p

Result = array( rep(0,ntemp*ntemp), dim=c(ntemp, ntemp))

P = Lambda %*% Phi %*% t(Lambda)
for (i in 1:p) {
P[i,i] = 1
}

Residual = R.Sample - P


FirstDerivative = DifS2LPhiPsi (Lambda, Phi)
SecondDerivative = Dif2S2LPhiPsi (Lambda, Phi)


for (j in 1:ntemp) {
  for (i in 1:ntemp) {
    Result[i,j] = - sum ( Residual * SecondDerivative[1:p,1:p,i,j] ) + sum ( FirstDerivative[1:p,1:p,i] * FirstDerivative[1:p,1:p,j] ) 
  } # i
} # j

Result * 2

} # SHessianOLS

#....................................................................


#..................................................................................

# DifS2LPhiPsi computes the partial derivatives of the manifest variable
# covariance matrix WRT factor loadings, factor correlations, and unique variances.
# Lambda <- (p,m)
# Phi <- (m,m), a symmetric matrix
# Result <- (p,p,(p*m + m*(m-1)/2 + p))

DifS2LPhiPsi <- function(Lambda, Phi){

# The function invokes no other external functions.

p = dim(Lambda)[1]
m = dim(Lambda)[2]

Result = array( rep(0, p*p*( p*m + m*(m-1)/2 + p)), dim=c(p,p,( p*m + m*(m-1)/2 + p) ))

## Factor loadings

ij = 0
LPhi = Lambda %*% Phi

for (j in 1:m) {
  for (i in 1:p) {
     ij = ij + 1
     Result[i,1:p,ij] = LPhi[1:p,j]
     Result[1:p,i,ij] = LPhi[1:p,j] 
  } # i
 } # j

## Factor correlations

ij = p*m

if (m>1) {
for (j in 2:m) {
  for (i in 1:(j-1)) {
    ij = ij + 1    
     Temp = Lambda[1:p,i] %*% t(Lambda[1:p,j]) 
     Result[1:p,1:p,ij] = Temp + t(Temp)
   } # i
 } # j
}
  
## unique variances
ij = p*m + m*(m-1)/2

for (i in 1:p) {
ij = ij + 1
Result[i,i,ij] = 1
}
 

## Output

Result

} # DifS2LPhiPsi


#..................................................................................

### Dif2S2LPhiPsi computes second order derivatives of P WRT parameters
# Lambda <- (p,m)
# Phi <- (m,m), a symmetric matrix
# Result <- (p,p,(p*m+m*(m-1)/2+p) , (p*m+m*(m-1)/2+p) )

Dif2S2LPhiPsi <- function(Lambda, Phi){

# The function invokes no other external functions.


p = dim(Lambda)[1]
m = dim(Lambda)[2]

Result = array( rep(0, p*p*(p*m+m*(m-1)/2+p)^2), dim=c(p,p,(p*m+m*(m-1)/2+p) , (p*m+m*(m-1)/2+p) ))

## Step 1, d^2 S / d L^2

kl =0

for (l in 1:m) {
  for (k in 1:p) {
     kl = kl + 1
     ij = 0
    for (j in 1:m) {
       for (i in 1:p) {
         ij = ij + 1
         Result[i,k,ij,kl] = Phi[j,l]
         Result[k,i,ij,kl] = Phi[j,l]           
         Result[i,k,kl,ij] = Phi[j,l]
         Result[k,i,kl,ij] = Phi[j,l]           
      } # i
    } # j
  } # k
} # l

## Step 2, d^2 S / d L d Phi'

kl = p*m 

if (m>1) {

for (l in 2:m) {
  for (k in 1:(l-1)) {
    # k and l refer to factor correlations in Step 2
    kl = kl + 1
    ij=0
    for (j in 1:m) {
      for (i in 1:p) {
        ij = ij + 1
          Temp = array(rep(0, p*p), dim=c(p,p)) 
          if (j==k) Temp[i, 1:p] = Lambda[1:p,l]
          if (j==l) Temp[i, 1:p] = Lambda[1:p,k]
          Result[1:p,1:p,ij,kl] = Temp + t(Temp)
          Result[1:p,1:p,kl,ij] = Temp + t(Temp)          
      } # i
    } # j        
   } # k
 } # l
}

## Step 3, All other second order derivatives are zero. 


## Output

Result

} # Dif2S2LPhiPsi
# ....................................................................................


if (is.null(extraction)) extraction='ml'


if (extraction=='ml') {

Hessian.Analytic = SHessianML(Lambda, Phi, RSample)

} else if (extraction=='ols') {

Hessian.Analytic = SHessianOLS(Lambda, Phi, RSample)

} else {
  stop ("wrong specification for the factor extraction method")
}

Hessian.Analytic

} # EFA.Hessian

#############################################################################


