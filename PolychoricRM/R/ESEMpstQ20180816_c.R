### 2018-07-31, Tuesday, Guangjian Zhang
### SEMpstQ standards for Exploratory Structural Equation Modeling with partially specified target, oblique rotation




###########################################################################
### The function "GPFoblq.RCPhi" conducts exploratory SEM as an oblique rotation method. 
### The function is modified from "GPFoblq" in the R package "GPArotation",


ESEMpstQ = function (A, Tmat = diag(ncol(A)), normalize = FALSE, eps = 1e-05, 
    maxit = 1000, method = "quartimin", methodArgs = NULL, BGTarget=NULL, BGWeight=NULL, PhiTarget = NULL, PhiWeight = NULL,  wxt2 = 1e0) 
{


#---------------------------------------------------------------------------------
vgQ.pst <- function(L, W=NULL, Target=NULL){
   if(is.null(W))      stop("argument W must be specified.")
   if(is.null(Target)) stop("argument Target must be specified.")
  
   Btilde <- W * Target
   list(Gq= 2*(W*L-Btilde), 
        f = sum((W*L-Btilde)^2),
        Method="Partially specified target")
}
#---------------------------------------------------------------------------------

#------------------------------------------
# D.Phi.2.T computes derivatives of Phi WRT T

D.Phi.2.T <- function(T) {


m = dim(T)[1]

Result = array(rep(0,m*m*m*m), dim=c(m,m,m,m))

for (j in 1:m) {
  for (i in 1:m) {
   if (i != j) {
     Result[i,j,1:m,i] = T[1:m,j]
     Result[i,j,1:m,j] = T[1:m,i]
   }
  } # (i in 1:m)
} # (j in 1:m)

Result

} # D.Phi.2.T

#-------------------------------------------------------------------
# D.BG.2.Phi computes derivatives of BG WRT Phi

D.BG.2.Phi <- function(T,m1) {

m = dim(T)[1]


Result = array(rep(0,m1*m*m*m), dim=c(m1,m,m,m))

Phi = t(T) %*% T

 for (i in 1:m1) {

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


#----------------------------------------------------------------------
# D.BG.2.T computes the derivatives of BG WRT T using the chain rule

D.BG.2.T <- function(BG.2.Phi,Phi.2.T) {

m1 = dim(BG.2.Phi)[1]
m = dim(BG.2.Phi)[2]

Result = array(rep(0,m1*m*m*m), dim=c(m1,m,m,m))

for (j in 1:m) {
  for (i in 1:m1) {
   for (l in 1:m) {
     for (k in 1:m) {

    Result[i,j,k,l] = sum(BG.2.Phi[i,j,1:m,1:m] * Phi.2.T[1:m,1:m,k,l])

    }
   }
 }
}

Result

} # D.BG.2.T

#---------------------------------------------------------------------------------


#---------------------------------------------------------------------------------

### The function "vgQ.ESEM" conducts exploratory SEM as a rotation method.
### The function is modified from vgQ.pst in the R package "GPArotation"

vgQ.ESEM <- function(Transform, BGWeight=NULL, BGTarget=NULL, PhiW=NULL, PhiTarget=NULL, wxt2 = 1e0){


   if(is.null(PhiW))      stop("argument PhiW must be specified.")
   if(is.null(PhiTarget)) stop("argument PhiTarget must be specified.")


   if (max(abs(PhiTarget - t(PhiTarget)))>1.0e-10) stop(" PhiTarget must be symmetric.")
   if (max(abs(PhiW - t(PhiW)))>1.0e-10) stop(" PhiW must be symmetric.")



   m  = dim(Transform)[1]
   m2 = dim(PhiTarget)[1]
   m1 = m - m2
   
   if((is.null(BGWeight)) & (m1>0) ) stop("argument BGWeight must be specified.")
   if((is.null(BGTarget)) & (m1>0)) stop("argument BGTarget must be specified.")


   Phi = t(Transform) %*% Transform

   dQ2T = matrix(0,m,m)
   f.BG = 0
   f.Phi.xi = 0
    
## BG

   
   

   
   if (m1 > 0) {
  
   BGtilde <- BGWeight * BGTarget

   BG = matrix(0,m1,m)

   for (i in 1:m1) {
   BG[i,(i+1):m] =  solve(t(Phi[(i+1):m, (i+1):m]), Phi[i,(i+1):m])
   BG[i,i] = Phi[i,i] - sum(Phi[i,(i+1):m] * BG[i,(i+1):m])
   }
   

   Gq2BG = 2*(BGWeight * BG - BGtilde)
   f.BG = sum((BGWeight * BG - BGtilde)^2) 


   BG2Phi  = D.BG.2.Phi(Transform,m1)
   Phi2T   = D.Phi.2.T(Transform)
   BG2T    = D.BG.2.T (BG2Phi,Phi2T)


   for (i in 1:m1) {
    for (j in (i+1):m) {
       if (BGWeight[i,j] == 1) {
          dQ2T = dQ2T + BG2T[i,j,1:m,1:m] * Gq2BG[i,j]
       }
     }
   }

  } # (m1 > 0)



   
   Phi.xi = Phi[(m1+1):m,(m1+1):m]   
   Phitilde.xi <- PhiW * PhiTarget
   
   Gq2phi.xi = 2*(PhiW * Phi.xi - Phitilde.xi)
   f.Phi.xi = sum((PhiW * Phi.xi - Phitilde.xi)^2) / 2  
   
if (m2>1) {
  
   for (j2 in 2:m2) {
      for (i2 in 1:j2) {
        if (PhiW[i2,j2]==1) { 
         i = m1 + i2
         j = m1 + j2
         dQ2T[1:m,i] = dQ2T[1:m,i] + Transform[1:m,j] * Gq2phi.xi[i2,j2] 
         dQ2T[1:m,j] = dQ2T[1:m,j] + Transform[1:m,i] * Gq2phi.xi[i2,j2] 
        }
   }
  }

} # (m2>1)

  Method="Saturated SEM"   


   list(dQ2T = dQ2T * wxt2, 
        f.ESEMT = (f.BG + f.Phi.xi) * wxt2,
        Method=Method)
} # vgQ.ESEM


# ---------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------

    if (1 >= ncol(A)) 
        stop("rotation does not make sense for single factor models.")
    if ((!is.logical(normalize)) || normalize) {
        A2 = A * A
        Com = rowSums(A2)
        W = sqrt(Com) %*% matrix(1,1,ncol(A))
        normalize <- TRUE
        A <- A/W
    }
    al <- 1
    L <- A %*% t(solve(Tmat))
    Method <- paste("vgQ", method, sep = ".")
    VgQ <- do.call(Method, append(list(L), methodArgs))
    G1 <- -t(t(L) %*% VgQ$Gq %*% solve(Tmat))
    f1 <- VgQ$f

    VgQ.2 = vgQ.ESEM(Tmat, BGWeight, BGTarget, PhiWeight, PhiTarget, wxt2)
    f = f1 + VgQ.2$f.ESEMT
    G = G1 + VgQ.2$dQ2T
    Table <- NULL
    VgQt <- do.call(Method, append(list(L), methodArgs))
    VgQ.2 = vgQ.ESEM(Tmat, BGWeight, BGTarget, PhiWeight, PhiTarget, wxt2)

    for (iter in 0:maxit) {
        Gp <- G - Tmat %*% diag(c(rep(1, nrow(G)) %*% (Tmat * 
            G)))
        s <- sqrt(sum(diag(crossprod(Gp))))
        Table <- rbind(Table, c(iter, f, log10(s), al))
        if (s < eps) 
            break
        al <- 2 * al
        for (i in 0:10) {
            X <- Tmat - al * Gp
            v <- 1/sqrt(c(rep(1, nrow(X)) %*% X^2))
            Tmatt <- X %*% diag(v)
            L <- A %*% t(solve(Tmatt))
            VgQt <- do.call(Method, append(list(L), methodArgs))
            VgQ.2 = vgQ.ESEM(Tmatt, BGWeight, BGTarget, PhiWeight, PhiTarget, wxt2)


            improvement <- f - ( VgQt$f + VgQ.2$f.ESEMT )
            if (improvement > 0.5 * s^2 * al) 
                break
            al <- al/2
        }
        Tmat <- Tmatt
        f1 <- VgQt$f
        G1 <- -t(t(L) %*% VgQt$Gq %*% solve(Tmatt))

    VgQ.2 = vgQ.ESEM(Tmatt, BGWeight, BGTarget, PhiWeight, PhiTarget, wxt2)
    f = f1 + VgQ.2$f.ESEMT
    G = G1 + VgQ.2$dQ2T                                     


    }
    convergence <- (s < eps)
    if ((iter == maxit) & !convergence) 
        warning("convergence not obtained in GPFoblq. ", maxit, 
            " iterations used.")
    if (normalize) 
        L <- L * W
    dimnames(L) <- dimnames(A)
    r <- list(loadings = L, Phi = t(Tmat) %*% Tmat, Th = Tmat, 
        Table = Table, method = VgQ$Method, orthogonal = FALSE, 
        convergence = convergence, Gq = VgQt$Gq)
    class(r) <- "GPArotation"
    r
} # ESEMpstQ

##########################################################################################
