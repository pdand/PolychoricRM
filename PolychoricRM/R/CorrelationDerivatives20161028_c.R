### Two external functions: EliU, D.g.2.r


EliU <- function(MP, eta = 1) {


## Browne, M. W. & Shapiro, A. (1986). The asymptotic Covariance matrix of 
## sample correlation coefficients under general conditions. Linear Algebra
## and its applications, 82, 169-176.
## Equations (4.1) and (4.3)

  

## EliU -> The asymptotic covariance matrix of sample correlations if manifest variables are
## of an elliptical distribution.

# It does not require any external functions.

p = dim(MP)[1]

Ms= matrix(0,p*p,p*p)

for (j in 1:p) {
 for (i in 1:p)  {
 
   if (j==i) {
   ii = (i-1)*p + i
   Ms[ii,ii] = 1
   } else
   {
   ij = (j-1)*p + i
   ji = (i-1)*p + j
   Ms[ij,ij] = 0.5
   Ms[ij,ji] = 0.5
   }

  } # i
} # j


Kd = matrix(0,p*p,p)
for (i in 1:p) {
 ii = (i-1) * p + i
 Kd[ii,i] = 1 
}


A = Ms %*% (MP %x% diag(p)) %*% Kd

Gamma = 2 * Ms %*% (MP %x% MP)

if (eta != 1 ) {
MP.v = array(MP)
Gamma = eta * Gamma + (eta -1) * outer(MP.v, MP.v)
}

B = Gamma %*% Kd
G = t(Kd) %*% Gamma %*% Kd

Cov.r = Gamma - A %*% t(B) - B %*% t(A) + A %*% G %*% t(A)

} # EliU


#################################################################
 
D.g.2.r <- function(Lambda, Phi, extraction=NULL) {

### use DifS2LPhiPsi


#-----------------------------------------------------------------------------------

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
} # if (m>1)
  
## unique variances
ij = p*m + m*(m-1)/2

for (i in 1:p) {
ij = ij + 1
Result[i,i,ij] = 1
}
 

## Output

Result

} # DifS2LPhiPsi


#-------------------------------------------------------------------------------------


if (is.null(extraction)) extraction='ml'
p = dim(Lambda)[1]
m = dim(Lambda)[2]
Nq = p * m + m * (m - 1) / 2 + p


PM = Lambda %*% Phi %*% t(Lambda)
for (i in 1:p) {
PM[i,i] = 1
}

PInverse = solve(PM)

Gradient = DifS2LPhiPsi (Lambda, Phi) # An external function

D.Gradient.2.R.vector = matrix(0,p*p,Nq)

if (extraction=='ml') {

for (i in 1:Nq) {
D.Gradient.2.R.vector[, i] = array((-PInverse %*% Gradient[1:p,1:p,i] %*% PInverse), dim=c(p*p,1)) 
}

} else if (extraction=='ols') {

for (i in 1:Nq) {
D.Gradient.2.R.vector[, i] = - 2 * array( Gradient[1:p,1:p,i], dim=c(p*p,1)) 
} 

} else {
  stop ("wrong specification for the factor extraction method")
}

D.Gradient.2.R.vector

} # D.g.2.r

###############################################################################################