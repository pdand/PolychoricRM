
ssem.se.augmt <- function(Lambda, Phi, Rsample, N, extraction, rotation, normalize=FALSE, modelerror, geomin.delta=NULL, MTarget=NULL, MWeight=NULL,
                          BGWeight = NULL, BGTarget = NULL,PhiWeight = NULL, PhiTarget = NULL, u.r = NULL, acm.type, wxt2=1e0) {


# It contains two internal functions D.G.2.Phi and D2.G.2.Phi.N



#----------------------------------------------------------------------

#-------------------------------------------------------------------
# D.BG.2.Phi computes derivatives of BG WRT Phi

D.BG.2.Phi <- function(Phi,m1) {

# Phi -> the factor correlation matrix, m by m, real
# m1, the number of endogenous factors, input

m = dim(Phi)[1]

Result = array(rep(0,m1*m*m*m), dim=c(m1,m,m,m))

# Phi = t(T) %*% T

 for (i in 1:m1) {   # Note that we need to compute the derivative ROW by ROW rather than column by column.

 Phi.inv = solve(Phi[(i+1):m, (i+1):m])
 w.i = Phi[i,(i+1):m] %*% Phi.inv


 Result[i,i,i,(i+1):m] = - w.i * 2
 Result[i,i,(i+1):m, (i+1):m] =  t(w.i) %*% w.i * (matrix(1,m-i,m-i) - diag(m - i)) 


 for (j in (i+1):m) {

 Result[i,j,i,(i+1):m] = Phi.inv[(j-i),1:(m-i)]
 Result[i,j,(i+1):m, (i+1):m] = - t(w.i) %*% Phi.inv[(j-i),1:(m-i)] * (matrix(1,m-i,m-i) - diag(m - i)) 

  } # (j in (i+1):m)


 } # (i in 1:m1)

Result

} # D.BG.2.Phi


#----------------------------------------------------------------------------------------

# D.G.2.Phi differentiates Q WRT Phi.

D.G.2.Phi <- function(Phi, BGWeight , BGTarget, PhiW , PhiTarget) { 

# Phi -> the factor correlation matrix, m by m, real
# BGWeight -> the weight matrix of [B|G], m1 by m, real
# BGTarget -> the target matrix of [B|G], m1 by m, real
# PhiW     -> the weight matrix of Phi_xi, m2 by m2,real
# PhiTarget-> the target matrix of Phi_xi, m2 by m2,real

m = dim(Phi)[1]
m2 = dim(PhiW)[1]
m1 = m - m2


dQ2Phi = matrix(0,m,m)


if (m1 > 0) {

   BGtilde <- BGWeight * BGTarget

   BG = matrix(0,m1,m)

   for (i in 1:m1) {
   BG[i,(i+1):m] =  solve(t(Phi[(i+1):m, (i+1):m]), Phi[i,(i+1):m])
   BG[i,i] = Phi[i,i] - sum(Phi[i,(i+1):m] * BG[i,(i+1):m])
   }
   

   Gq2BG = 2*(BGWeight * BG - BGtilde)

#   f.BG = sum((BGWeight * BG - BGtilde)^2) 


   BG2Phi  = D.BG.2.Phi(Phi,m1)

   

   for (i in 1:m1) { # Row by Row
    for (j in (i+1):m) {
       if (BGWeight[i,j] == 1) {
          dQ2Phi = dQ2Phi + BG2Phi[i,j,1:m,1:m] * Gq2BG[i,j]
       }
     }
   }

} ## (m1 > 0)

## The targets on factor correlations among exogenous variables


if ( m2 > 1 ) {
   
   Phi.xi = Phi[(m1+1):m,(m1+1):m]   
   Phitilde.xi <- PhiW * PhiTarget
   
   Gq2phi.xi = 2*(PhiW * Phi.xi - Phitilde.xi)

#   f.Phi.xi = sum((PhiW * Phi.xi - Phitilde.xi)^2) / 2  

  
   for (j2 in 1:m2) {
      for (i2 in 1:m2) {
        if (PhiW[i2,j2]==1) { 
         i = m1 + i2
         j = m1 + j2
         dQ2Phi[i,j] = dQ2Phi[i,j] + Gq2phi.xi[i2,j2] / 2           
        }
   }
  }

} # ( m2 > 1 )

Result = dQ2Phi

} # D.G.2.Phi

#---------------------------------------------------------------------------------------

D2.G.2.Phi.N <- function(Phi, BGWeight , BGTarget, PhiW , PhiTarget, epsilon = 0.0001) { 

# D2.G.2.Phi.N computes the second derivatives of constraint functions with regard to Phi
# The method is a numerical method.


# Phi -> the factor correlation matrix, m by m, real
# BGWeight -> the weight matrix of [B|G], m1 by m, real
# BGTarget -> the target matrix of [B|G], m1 by m, real
# PhiW     -> the weight matrix of Phi_xi, m2 by m2,real
# PhiTarget-> the target matrix of Phi_xi, m2 by m2,real

m = dim(Phi)[1]

Result = array(rep(0,m**4), dim=c(m,m,m,m))

for (j in 2:m) {
  for (i in 1:(j-1)) {

     IJ = matrix(0,m,m)
     IJ[i,j] = epsilon
     IJ[j,i] = epsilon
     Phi.ij = Phi + IJ

    T.D.G.2.Phi.ij = D.G.2.Phi (Phi.ij, BGWeight , BGTarget, PhiWeight, PhiTarget)               
    M.Positive = T.D.G.2.Phi.ij + t(T.D.G.2.Phi.ij)     


    Phi.ij = Phi - IJ
    T.D.G.2.Phi.ij = D.G.2.Phi (Phi.ij, BGWeight , BGTarget, PhiWeight, PhiTarget)               
    M.Negative = T.D.G.2.Phi.ij + t(T.D.G.2.Phi.ij)     
     
    Result[1:m,1:m,i,j] = (M.Positive - M.Negative) / (2 * epsilon) 
   
  }
}

Result

} # D2.G.2.Phi.N

#--------------------------------------------------------------------------------------------



# We assume that the manifest variables are normally distributed for the time being.

if (is.null(rotation)) stop ("No rotaton criterion is specified for numberical approximation")
if ((rotation=='geomin') & (is.null(geomin.delta))) geomin.delta = 0.01
if ((rotation=='target') & ((is.null(MWeight)) | (is.null(MTarget)))) stop ("MWeight or MTarget is not specified for target rotation")
if ((rotation=='xtarget') & ((is.null(MWeight)) | (is.null(MTarget)) | (is.null(PhiWeight)) | (is.null(PhiTarget)) )) stop ("MWeight or MTarget is not specified for target rotation") # 2016-06-03, GZ
if (is.null(modelerror)) modelerror='YES' 

p = dim(Lambda)[1]
m = dim(Lambda)[2]
m2 = dim(PhiTarget)[1]
m1 = m - m2
Nc = m * (m-1) # the number of constraints
Nq = p * m + m * (m - 1) / 2 + p  # the number of parameters

PM = Lambda %*% Phi %*% t(Lambda) 
PM = PM - diag(diag(PM)) + diag(p)

##

if (modelerror== 'NO') Rsample = PM

if (is.null(u.r)) u.r = EliU(Rsample)

# if (is.null(u.r)) {

# if (modelerror== 'YES') {
#  u.r = EliU(Rsample) 
#} else if (modelerror == 'NO') {
#  u.r = EliU(PM)
#} else {
# stop("Model Error option is inappropriately specified.")
#}
# } ## (u.r == NULL)

##

if (extraction == 'ml') {

dg2r = D.g.2.r (Lambda, Phi, extraction='ml')
Hessian = EFA.Hessian (Lambda, Phi, Rsample, extraction='ml') ####


} else if (extraction == 'ols') {

dg2r = D.g.2.r (Lambda, Phi, extraction='ols')
Hessian = EFA.Hessian (Lambda, Phi, Rsample, extraction='ols') #####


} else {

stop ('Factor extraction method is incorrectly specified.')
}


##


if (acm.type==1) {
  Y.Hat = u.r
} else {
  
  
  Y.Hat = matrix(0, (p*(p-1)/2), (p*(p-1)/2))
  u.r.col = matrix(0, (p*p), (p*(p-1)/2))
  
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

dg2r.upper = matrix(0, (p*(p-1)/2), Nq)


ij = 0
for (j in 2:p) {
  for (i in 1:(j-1)) {
    ij = ij + 1
    dg2r.upper[ij,] = dg2r[((j-1)*p + i) , ]
  }
}

Ham = (t(dg2r.upper) %*% Y.Hat %*% dg2r.upper) * 4 # Multiplying by 4 to accommodate the partial derivatives



# Ham = t(dg2r) %*% u.r %*% dg2r

###

if (rotation=='CF-varimax') {

Olq.Con.Parameters = Extended.CF.Family.c.2.LPhi (Lambda, Phi, 'CF-varimax',normalize)

} else if (rotation=='CF-quartimax') {

Olq.Con.Parameters = Extended.CF.Family.c.2.LPhi (Lambda, Phi, 'CF-quartimax',normalize)

} else if (rotation=='CF-facparsim') {
  
  Olq.Con.Parameters = Extended.CF.Family.c.2.LPhi (Lambda, Phi, 'CF-facparsim',normalize)
  
} else if (rotation=='CF-equamax') {
  
  Olq.Con.Parameters = Extended.CF.Family.c.2.LPhi (Lambda, Phi, 'CF-equamax',normalize)
  
} else if (rotation=='CF-parsimax') {
  
  Olq.Con.Parameters = Extended.CF.Family.c.2.LPhi (Lambda, Phi, 'CF-parsimax',normalize)

} else if (rotation=='geomin') {

if (is.null(geomin.delta)) geomin.delta = 0.01

Olq.Con.Parameters = Derivative.Constraints.Numerical(Lambda,Phi,'geomin',normalize,geomin.delta)

} else if (rotation=='target') {

Olq.Con.Parameters = Derivative.Constraints.Numerical(Lambda,Phi,'target',normalize,MWeight=MWeight, MTarget=MTarget)

} else if (rotation=='xtarget') {

Olq.Con.Parameters = Derivative.Constraints.Numerical(Lambda,Phi,'xtarget',normalize,MWeight=MWeight, MTarget=MTarget,PhiWeight = PhiWeight, PhiTarget = PhiTarget, wxt2 = wxt2) # 2017-11-28, GZ!


} else if (rotation=='semtarget') {

Olq.Con.Parameters = Derivative.Constraints.Numerical(Lambda,Phi,'target',normalize,MWeight=MWeight, MTarget=MTarget)

T.D.G.2.Phi = D2.G.2.Phi.N (Phi, BGWeight , BGTarget, PhiWeight , PhiTarget, epsilon = 0.0001)

ICon = 0
for (j in 1:m) {
  for (i in 1:m) {
    if (i != j) {
    ICon = ICon + 1
    Itemp = p * m
    for (l in 2:m) {
      for (k in 1: (l-1)) {
    Itemp = Itemp + 1
    Olq.Con.Parameters[ICon, Itemp] = Olq.Con.Parameters[ICon, Itemp] - T.D.G.2.Phi[i,j,k,l]
       } # k
     } # l
   } # (i != j)
  } # i
} # j



} else {

  stop ("wrong specification for the factor rotation criterion")
}

### 


Temp.Bigger = matrix(0,(Nq+Nc),(Nq+Nc))
Temp.Bigger[1:Nq,1:Nq] = Hessian
Temp.Bigger[1:Nq,(Nq+1):(Nq+Nc)] = t(Olq.Con.Parameters)
Temp.Bigger[(Nq+1):(Nq+Nc),1:Nq] = Olq.Con.Parameters


Temp.Bigger.inverse = solve(Temp.Bigger)

Sandwich = Temp.Bigger.inverse[1:Nq,1:Nq] %*% Ham %*% Temp.Bigger.inverse[1:Nq,1:Nq]

###

SE = sqrt(diag(Sandwich)/(N-1))

Lambda.se <- array(SE[1: (p*m) ],dim=c(p,m))

Phi.se <- matrix(0,m,m)
ij=0
for (j in 2:m) {
 for (i in 1:(j-1)) {
    ij = ij + 1
    Phi.se[i,j] = SE[p*m+ij]
  }
}

Phi.se = Phi.se + t(Phi.se)

ME.se <- SE[(p*m+m*(m-1)/2+1):Nq]



if (m1>0) {


D.BG.Phi.2 = array(rep(0,m1*m*m*m), dim=c(m1,m,m,m))
BG = matrix(0,m1,m)

 for (i in 1:m1) {   # Note that we need to compute the derivative ROW by ROW rather than column by column.

 Phi.inv = solve(Phi[(i+1):m, (i+1):m])
 w.i = Phi[i,(i+1):m] %*% Phi.inv
 BG[i,(i+1):m] = w.i
 BG[i,i] = Phi[i,i] - sum(w.i * Phi[i,(i+1):m])

 D.BG.Phi.2[i,i,i,(i+1):m] = - w.i * 2
 D.BG.Phi.2[i,i,(i+1):m, (i+1):m] =  t(w.i) %*% w.i * (matrix(1,m-i,m-i) - diag(m - i)) 


 for (j in (i+1):m) {

 D.BG.Phi.2[i,j,i,(i+1):m] = Phi.inv[(j-i),1:(m-i)]
 D.BG.Phi.2[i,j,(i+1):m, (i+1):m] = - t(w.i) %*% Phi.inv[(j-i),1:(m-i)] * (matrix(1,m-i,m-i) - diag(m - i)) 

  } # (j in (i+1):m)


 } # (i in 1:m1)

### Compute standard errors for BG
ACM.phi = Sandwich[(p*m+1):(p*m+(m*(m-1)/2)),(p*m+1):(p*m+(m*(m-1)/2))]
ACM.BG = matrix(0,m1,m)

for (i in 1:m1) {
 for (j in i:m) {
   
     temp.d.phi = rep(0,m*(m-1)/2)
     kl = 0 
     for (l in 2:m) {
       for (k in 1:(l-1)) {
     kl = kl + 1
     temp.d.phi[kl] = D.BG.Phi.2[i,j,k,l] + D.BG.Phi.2[i,j,l,k]
       } # k
     } # l
    ACM.BG[i,j] = sum((ACM.phi %*% temp.d.phi) * (temp.d.phi))

 } # j
} # i


BG.se = sqrt(ACM.BG / (N-1))
  psi = rep(0,m1)
  psi.se = rep(0,m1)

  for (i in 1:m1) {
  psi[i] = BG[i,i]
  psi.se[i] = BG.se[i,i]
  BG[i,i] = 0
  BG.se[i,i] = 0
  } # i

} else {

BG = NULL
BG.se = NULL
psi = NULL
psi.se = NULL


} # (m1>0)



if (m2>1) {
Phi.xi = Phi[(m1+1):m,(m1+1):m]
Phi.xi.se = Phi.se[(m1+1):m,(m1+1):m]
} else {
Phi.xi = NULL
Phi.xi.se = NULL
} # (m2>1)





##
## list(Lambda.se = Lambda.se, Phi.se = Phi.se, Psi.se = Psi.se, Temp.Bigger = Temp.Bigger, Hessian = Hessian)

list(Lambda.se = Lambda.se, Phi.se = Phi.se, ME.se = ME.se, BG = BG, BG.se=BG.se, psi = psi, psi.se=psi.se, Phi.xi = Phi.xi, Phi.xi.se = Phi.xi.se )

} # oblq.se.augmt

