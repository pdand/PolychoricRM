
Derivative.Orth.Constraints.Numerical <- function (Lambda, rotation=NULL,normalize=FALSE,geomin.delta = NULL, MWeight=NULL, MTarget=NULL) {



vgQ.cf <- function (L, kappa = 0) 
{
    k <- ncol(L)
    p <- nrow(L)
    N <- matrix(1, k, k) - diag(k)
    M <- matrix(1, p, p) - diag(p)
    L2 <- L^2
    f1 <- (1 - kappa) * sum(diag(crossprod(L2, L2 %*% N)))/4
    f2 <- kappa * sum(diag(crossprod(L2, M %*% L2)))/4
    list(Gq = (1 - kappa) * L * (L2 %*% N) + kappa * L * (M %*% 
        L2), f = f1 + f2, Method = paste("Crawford-Ferguson:k=", 
        kappa, sep = ""))
} # vgQ.cf


###-------------------------------------------------------------------------

vgQ.pst <- function(L, W=NULL, Target=NULL){
   if(is.null(W))      stop("argument W must be specified.")
   if(is.null(Target)) stop("argument Target must be specified.")
   Btilde <- W * Target
   list(Gq= 2*(W*L-Btilde), 
        f = sum((W*L-Btilde)^2),
        Method="Partially specified target")
}

###---------------------------------------------------------------------------

vgQ.geomin <- function(L, delta=.01){
  k <- ncol(L)
  p <- nrow(L)
  L2 <- L^2 + delta
  pro <- exp(rowSums(log(L2))/k) 
  list(Gq=(2/k)*(L/L2)*matrix(rep(pro,k),p),
       f= sum(pro), 
       Method="Geomin")
  }

#------------------------------------------------------------------------------------

DCon.Unrotated <- function (A, extraction) {

p = dim(A)[1]
m = dim(A)[2]
Nc = m*(m-1)/2
Nq = p*m + p

if (extraction == "ml") {

psi = 1 - diag(A %*% t(A))
Dpsi.inv = diag (1/psi)
A.star = Dpsi.inv %*% A

} else if (extraction == "ols") {

A.star = A

} else stop (paste("DCon.Unrotated Error:", extraction, "is wrongly specified for extraction."))


d.Con.Parameters = matrix(0,Nc,Nq)

uv=0
if (m>1) {
for (v in 2:m) {
  for (u in 1:(v-1)) {
   uv = uv + 1
   M.temp = matrix(0,p,m)
   M.temp[1:p,u] = A.star[1:p,v]
   M.temp[1:p,v] = A.star[1:p,u]
   d.Con.Parameters [uv, 1:(p*m)] = array(M.temp)
   d.Con.Parameters [uv, (p*m + 1):(p*m+p)] = - A.star[1:p,v] * A.star[1:p,u]
  } # u
} # v
}
  
d.Con.Parameters

} 
#--------------------------------------------------------------------------------------



if (is.null(rotation)) stop ("No rotaton criterion is specified for numberical approximation")


p = dim(Lambda)[1]
m = dim(Lambda)[2]
Nc = m * (m-1) / 2 # the number of constraints
Nq = p * m  + p  # the number of parameters


if (normalize) {
Lambda0 = Lambda
h = rowSums(Lambda0 **2)
h.half.inverse =  1 / sqrt(h)
h.onehalf.inverse = (h.half.inverse)**3
 for (k in 1:p) {
  Lambda[k,] = Lambda[k,] * h.half.inverse[k]
  } # (k in 1:p)

} # (normalize)



d.Con.Parameters = matrix(0, Nc , Nq)


if ((rotation == 'ml') | (rotation == 'ols')) { d.Con.Parameters = DCon.Unrotated(Lambda,rotation)

} else {

d.Con.Loading <- array (rep(0,m*m*p*m), dim=c(m,m,p,m))

Z = matrix(0,p,m)
eps = 1e-4
i = 1
j = 1


for (j in 1:m) {
  for (i in 1:p) {
    dZ = Z
    dZ[i,j] = eps


if (rotation == 'geomin') {
    if (is.null(geomin.delta)) geomin.delta = 0.01

    d.Con.Loading [1:m,1:m,i,j] =
    (crossprod((Lambda + dZ), (vgQ.geomin(Lambda + dZ, geomin.delta) $ Gq )) - 
    crossprod((Lambda - dZ), (vgQ.geomin(Lambda - dZ, geomin.delta) $ Gq ))) / (2*eps)  }

else if (rotation == 'target') {
    if ((is.null(MWeight)) | (is.null(MTarget))) stop ("MWeight or MTarget is not specified for target rotation")

    d.Con.Loading [1:m,1:m,i,j] =
    (crossprod((Lambda + dZ), (vgQ.pst(Lambda + dZ, MWeight, MTarget) $ Gq )) - 
    crossprod((Lambda - dZ), (vgQ.pst(Lambda - dZ, MWeight, MTarget) $ Gq ))) / (2*eps)  }

else if (rotation == 'CF-varimax') {
     cf.kappa = 1 / p
    d.Con.Loading [1:m,1:m,i,j] =
    (crossprod((Lambda + dZ), (vgQ.cf(Lambda + dZ, cf.kappa) $ Gq )) - 
    crossprod((Lambda - dZ), (vgQ.cf(Lambda - dZ, cf.kappa) $ Gq ))) / (2*eps)  }

else if (rotation == 'CF-quartimax') {
     cf.kappa = 0
    d.Con.Loading [1:m,1:m,i,j] =
    (crossprod((Lambda + dZ), (vgQ.cf(Lambda + dZ, cf.kappa) $ Gq )) - 
    crossprod((Lambda - dZ), (vgQ.cf(Lambda - dZ, cf.kappa) $ Gq ))) / (2*eps)  }

    
else if (rotation == 'CF-facparsim') {
      cf.kappa = 1
      d.Con.Loading [1:m,1:m,i,j] =
        (crossprod((Lambda + dZ), (vgQ.cf(Lambda + dZ, cf.kappa) $ Gq )) - 
           crossprod((Lambda - dZ), (vgQ.cf(Lambda - dZ, cf.kappa) $ Gq ))) / (2*eps)  }
    

else if (rotation == 'CF-equamax') {
      cf.kappa = m/(2*p)
      d.Con.Loading [1:m,1:m,i,j] =
        (crossprod((Lambda + dZ), (vgQ.cf(Lambda + dZ, cf.kappa) $ Gq )) - 
           crossprod((Lambda - dZ), (vgQ.cf(Lambda - dZ, cf.kappa) $ Gq ))) / (2*eps)  }
    

else if (rotation == 'CF-parsimax') {
      cf.kappa = (m-1)/(p+m-2)
      d.Con.Loading [1:m,1:m,i,j] =
        (crossprod((Lambda + dZ), (vgQ.cf(Lambda + dZ, cf.kappa) $ Gq )) - 
           crossprod((Lambda - dZ), (vgQ.cf(Lambda - dZ, cf.kappa) $ Gq ))) / (2*eps)  }
    

else {
  stop (paste(rotation," is not for numerical approximation"))
}


   } # (i in 1:p)
 }   # (j in 1:m)





ICon = 0
for (j in 2:m) {
  for (i in 1:(j-1))  {

       
      ICon = ICon + 1
      
      tempc2L = d.Con.Loading[i,j,1:p,1:m] - d.Con.Loading[j,i,1:p,1:m]
      

if (normalize) {

      temp = rowSums(tempc2L * Lambda0)     
      for (k in 1:p) {
      tempc2L[k,] = tempc2L[k,] * h.half.inverse[k] - temp[k] * h.onehalf.inverse[k] * Lambda0[k,]
      } # k

} 


      for (l in 1:m) {
         d.Con.Parameters[ICon,((l-1)*p+1):(l*p)] = tempc2L[1:p,l]            
      } # (l in 1:m)
 

  } # (i in 1:m)
} 

} 


d.Con.Parameters

} 

##############################################################
