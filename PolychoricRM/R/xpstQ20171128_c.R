

###########################################################################
### The function "GPFoblq.RCPhi" conduct oblique rotation with targets specified for factor correlations.
### The function is modified from "GPFoblq" in the R package "GPArotation",


xpstQ = function (A, Tmat = diag(ncol(A)), normalize = FALSE, eps = 1e-05, 
    maxit = 1000, method = "quartimin", methodArgs = NULL, PhiWeight = NULL, PhiTarget = NULL, wxt2 = 1e0) 
{

# A(p,m), the unrotated factor loading matrix, Input
# Tmat(m,m), the factor rotation matrix defined as in Jennrich2012Psychometrika, Input
# Normalize, logical, whether rows are standardized, Input
# esp, the user specified stopping criterion, Input
# maxit, the maximum number of iterations, Input
# method, characters, the rotation criterion, Input
# methodArgs, character, the target matrix and the weight matrix for factor loadings, Input
# PhiWeight(m,m), the weight matrix for factor correlations, Input
# PhiTarget(m,m), the target matrix for factor correlations, Input
# wxt2 = 1e0, the weight controling the contribution of phi toward to the whole target criterion function

# Loadings (p,m), the rotated factor loading matrix, Output
# Phi (m,m), the rotated factor correlation matrix, Output
# Th (m,m), the transformation matrix, Output
# Table, iteration details, Output

#---------------------------------------------------------------------------------
vgQ.pst <- function(L, W=NULL, Target=NULL){
   if(is.null(W))      stop("argument W must be specified.")
   if(is.null(Target)) stop("argument Target must be specified.")
   # Needs weight matrix W with 1's at specified values, 0 otherwise
   # e.g. W = matrix(c(rep(1,4),rep(0,8),rep(1,4)),8). 
   # When W has only 1's this is procrustes rotation
   # Needs a Target matrix Target with hypothesized factor loadings.
   # e.g. Target = matrix(0,8,2)
   Btilde <- W * Target
   list(Gq= 2*(W*L-Btilde), 
        f = sum((W*L-Btilde)^2),
        Method="Partially specified target")
}
#---------------------------------------------------------------------------------

### The function "vgQ.pstPhi" allows targets specified for factor correlations.
### The function is modified from vgQ.pst in the R package "GPArotation"

vgQ.pstPhi <- function(Transform, PhiW=NULL, PhiTarget=NULL, wxt2 = 1e0){

# Transform(m:m), the T matrix, defined as in Jennrich2002Psychometrika, Input
# PhiW(m,m), the weight matrix for Phi, input
# PhiTarget(m,m), the target matrix for Phi, Input
# wxt2 = 1e0, the weight controling the contribution of phi toward to the whole target criterion function, input

# dQ2T(m,m), the derivatives of the rotation crtieria with regard to T, Output
# f.Phi, the rotation criterion function value for Phi, only the UPPER triangular elements are considered, Output
# Method, character strings, output

   if(is.null(PhiW))      stop("argument W must be specified.")
   if(is.null(PhiTarget)) stop("argument Target must be specified.")
   if (max(abs(PhiTarget - t(PhiTarget)))>1.0e-10) stop(" PhiTarget must be symmetric.")
   if (max(abs(PhiW - t(PhiW)))>1.0e-10) stop(" PhiW must be symmetric.")

   # library(MASS)
   # Needs weight matrix W with 1's at specified values, 0 otherwise
   # e.g. W = matrix(c(rep(1,4),rep(0,8),rep(1,4)),8). 
   # When W has only 1's this is procrustes rotation
   # Needs a Target matrix Target with hypothesized factor loadings.
   # e.g. Target = matrix(0,8,2)
   Phi = t(Transform) %*% Transform
   Btilde <- PhiW * PhiTarget
   
   Gq2phi = 2*(PhiW * Phi - Btilde)
   f.Phi = sum((PhiW * Phi - Btilde)^2) / 2  
   Method="Partially specified target: Phi"   


   # Compute dQ2T
   m = dim(Transform)[1]   
   dQ2T = matrix(0,m,m)
  
   for (j in 2:m) {
      for (i in 1:j) {
        if (PhiW[i,j]==1) { 
         dQ2T[1:m,i] = dQ2T[1:m,i] + Transform[1:m,j] * Gq2phi[i,j] 
         dQ2T[1:m,j] = dQ2T[1:m,j] + Transform[1:m,i] * Gq2phi[i,j] 
        }
   }
  }



   list(dQ2T = dQ2T * wxt2, 
        f.Phi = f.Phi * wxt2,
        Method=Method)
} # vgQ.pstPhi
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
    VgQ.2 = vgQ.pstPhi(Tmat, PhiWeight, PhiTarget,wxt2)    
    f = f1 + VgQ.2$f.Phi                                   
    G = G1 + VgQ.2$dQ2T                                    
    Table <- NULL
    VgQt <- do.call(Method, append(list(L), methodArgs))   
    VgQ.2 = vgQ.pstPhi(Tmat, PhiWeight, PhiTarget,wxt2)    

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
           VgQ.2 = vgQ.pstPhi(Tmatt, PhiWeight, PhiTarget,wxt2)                  


            improvement <- f - ( VgQt$f + VgQ.2$f.Phi )     
            if (improvement > 0.5 * s^2 * al) 
                break
            al <- al/2
        }
        Tmat <- Tmatt
        f1 <- VgQt$f                                        
        G1 <- -t(t(L) %*% VgQt$Gq %*% solve(Tmatt))         

    VgQ.2 = vgQ.pstPhi(Tmatt, PhiWeight, PhiTarget,wxt2)     
    f = f1 + VgQ.2$f.Phi                                    
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
} # xpstQ

##########################################################################################
