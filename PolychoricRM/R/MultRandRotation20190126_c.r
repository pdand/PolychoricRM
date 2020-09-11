
MultRandRotation <- function(unrotated, epsilon = .01, nstart = 100, 
                             plot = TRUE, cex = .5,
                             rotation = 'geomin', rtype = 'oblique', 
                             normalize=F, MWeight=NULL, MTarget=NULL, 
                             PhiWeight = NULL, PhiTarget = NULL, wxt2 = 1){
  geomin.delta <- epsilon
  factors <- ncol(unrotated)
  p <- nrow(unrotated)
  

Make.Rot.Args <- function(rtype,rotation,normalize,p,m, geomin.delta, MTarget, MWeight,
                          PhiWeight, PhiTarget,wxt2=1e0,transformation=NULL) {
  maxit <- 100000 
  
  if ((rotation=='geomin') & (is.null(geomin.delta))) geomin.delta = 0.01
  if ((rotation=='target') & ((is.null(MWeight)) | (is.null(MTarget)))) stop ("MWeight or MTarget is not specified for target rotation")
  if ((rotation=='xtarget') & ((is.null(MWeight)) | (is.null(MTarget)) | (is.null(PhiWeight)) | (is.null(PhiTarget)) )) stop ("MWeight or MTarget is not specified for xtarget rotation") 
    
  if (rtype=='oblique') {
      
    if (rotation=='CF-varimax') {
      fnames = 'cfQ'
      rm("wxt2")
      Rot.Args <- list(Tmat=diag(m),kappa = 1/p, normalize=normalize,
                       eps=1e-6, maxit=maxit)   
        
      } else if (rotation=='CF-quartimax') {
        fnames = 'cfQ'
        rm("wxt2")
        Rot.Args <- list(Tmat=diag(m),kappa = 0, normalize=normalize, 
                         eps=1e-6, maxit=maxit)   
        
      } else if (rotation=='CF-facparsim') {
        fnames = 'cfQ'
        rm("wxt2")
        Rot.Args <- list(Tmat=diag(m),kappa = 1, normalize=normalize, 
                         eps=1e-6, maxit=maxit)   
        
      } else if (rotation=='CF-equamax') {
        fnames = 'cfQ'
        rm("wxt2")
        Rot.Args <- list(Tmat=diag(m),kappa = m/(2*p), normalize=normalize,
                         eps=1e-6, maxit=maxit)   
        
      } else if (rotation=='CF-parsimax') {
        fnames = 'cfQ'
        rm("wxt2")
        Rot.Args <- list(Tmat=diag(m),kappa = (m-1)/(p+m-2), normalize=normalize,
                         eps=1e-6, maxit=maxit)   
        
      } else if (rotation=='geomin') {
        fnames = 'geominQ'
        rm("wxt2")
        Rot.Args <- list(Tmat=diag(m),delta = geomin.delta, normalize=normalize, 
                         eps=1e-6, maxit=maxit)   
        
      } else if (rotation=='target') {
        fnames = 'pstQ'
        rm("wxt2")
        Rot.Args <- list(Tmat=diag(m), W = MWeight, Target=MTarget, 
                         normalize=normalize, eps=1e-6, maxit=maxit)   
        
      } else if (rotation=='xtarget') {
        fnames = 'xpstQ'
        Rot.Args <- list(Tmat=transformation, normalize=normalize, eps=1e-6, maxit=maxit,
                         method="pst", methodArgs = list(W = MWeight, Target = MTarget),
                         PhiWeight = PhiWeight, PhiTarget = PhiTarget, wxt2 = wxt2)   
      } else {
        stop (paste(rotation, ' has not been implemented yet.'))
      }
    } else {   
      
      if (rotation=='CF-varimax') {
        fnames = 'cfT'
        rm("wxt2")
        Rot.Args <- list(Tmat=diag(m),kappa = 1/p, normalize=normalize, 
                         eps=1e-6, maxit=maxit)   
        
      } else if (rotation=='CF-quartimax') {
        fnames = 'cfT'
        rm("wxt2")
        Rot.Args <- list(Tmat=diag(m),kappa = 0, normalize=normalize, 
                         eps=1e-6, maxit=maxit)   
        
      } else if (rotation=='CF-facparsim') {
        fnames = 'cfT'
        rm("wxt2")
        Rot.Args <- list(Tmat=diag(m),kappa = 1, normalize=normalize, 
                         eps=1e-6, maxit=maxit)   
        
      } else if (rotation=='CF-equamax') {
        fnames = 'cfT'
        rm("wxt2")
        Rot.Args <- list(Tmat=diag(m),kappa = m/(2*p), normalize=normalize, 
                         eps=1e-6, maxit=maxit)   
        
      } else if (rotation=='CF-parsimax') {
        fnames = 'cfT'
        rm("wxt2")
        Rot.Args <- list(Tmat=diag(m),kappa = (m-1)/(p+m-2), normalize=normalize, 
                         eps=1e-6, maxit=maxit)   
        
      } else if (rotation=='geomin') {
        fnames = 'geominT'
        rm("wxt2")
        Rot.Args <- list(Tmat=diag(m),delta = geomin.delta, normalize=normalize, 
                         eps=1e-6, maxit=maxit)   
        
      } else if (rotation=='target') {
        fnames = 'pstT'
        rm("wxt2")
        Rot.Args <- list(Tmat=diag(m), W = MWeight, Target=MTarget, 
                         normalize=normalize, eps=1e-6, maxit=maxit)   
        
      } else {
        stop (paste(rotation, ' has not been implemented yet.'))
      }
    } # End of the orthogonal rotations
    
    list(fnames = fnames, Rot.Args = Rot.Args)
  } # Make.Rot.Args
# Make.Rot.Args -----------------------------------
  
  Rot.Controls <- Make.Rot.Args(rtype, rotation, normalize, p, factors, geomin.delta,
                                MTarget, MWeight, PhiWeight, PhiTarget, wxt2)

  
  A <- unrotated               
  p <- dim(A)[1]
  m <- dim(A)[2] 
  K <- nstart                    # number of random starts
  Lam <- Phi <- ali <- list()   
  f.values <- rep(NA, K)        


  for(i in 1:K){

    RS <- qr.Q(qr(matrix(rnorm(m*m),m))) # random (m by m) orthogonal matrix
    A_start <- A %*% RS                  # random start (p by m) factor loadings

    Lambda.lst = do.call (Rot.Controls$fnames, append(list(A_start), Rot.Controls$Rot.Args))


    Lam[[i]] <- Lambda.lst$loadings 
    if(rtype=='oblique') Phi[[i]] <- Lambda.lst$Phi 
    if(rtype=='orthogonal') Phi[[i]] <- diag(m)
    
    f.values[i] <- min(Lambda.lst$Table[,2]) 

  }

  
  # Determine the number of solutions =================================
  deviation <- matrix(0, K, K)
  for(i in 1:K){
    if(i < K)
    for(j in (i+1):K){
      deviation[j,i] <- Align.Matrix(Order.Matrix = Lam[[i]], 
                                     Input.Matrix = rbind(Lam[[j]],diag(m)))[(p+m+1),1]
      deviation[i,j] <- deviation[j,i]
    }
  }
  diag(deviation) <- 1
  
  solution <- rep(0, K)
  for(j in 1:K){
    first.s <- which(solution==0)[1]
    solution[first.s] <- j
    for(i in 1:K){
      if(solution[i] == 0)
      if(round(deviation[first.s,i],3) == 0) solution[i] <- j
    }
  }
  
  # Check the number of solutions also using rotation criterion values -------
  solution2 <- rep(0, K)
  for(j in 1:K){
    first.s <- which(solution2==0)[1]
    solution2[first.s] <- j 
    for(i in 1:K){
      if(solution2[i] == 0)
      if(round(f.values[first.s] - f.values[i], 3) == 0) solution2[i] <- j
    }
  }

  
  # Basic descriptive outputs ==========================================
  solution <- solution2 
  nsol <- max(solution) 

  f.solutions <- rep(NA, nsol)
  for(j in 1:nsol) f.solutions[j] <- f.values[which(solution == j)[1]]
  count.sol <- count(solution)[,2]    # here we need library(plyr)
  count.sol <- count.sol[order(f.solutions)]
  f.solutions <- round(f.solutions[order(f.solutions)],3)
  
  # Align solutions =================================================
  LAMBDAS <- PHIS <- list()
  global.solution.number <- which(round(f.values,3) == f.solutions[1])[1]

  global.lambda <- Lam[[global.solution.number]]
  global.phi <- Phi[[global.solution.number]]
  LAMBDAS[[1]] <- round(global.lambda,3)
  PHIS[[1]] <- round(global.phi,3)

  if(nsol != 1){ 
    for(j in 2:nsol){
      solution.number <- which(round(f.values,3) == f.solutions[j])[1]
    
      lambdas <- Lam[[solution.number]]
      phis <- Phi[[solution.number]]
    
      aligned.lambdas2 <- Align.Matrix(Order.Matrix = global.lambda, 
                                       Input.Matrix = rbind(lambdas,phis))
      LAMBDAS[[j]] <- round(aligned.lambdas2[1:p,], 3)
      PHIS[[j]] <- round(aligned.lambdas2[-c(1:p,p+m+1),], 3)
    }
  } 
  
  # Make a bar plot =======================================================

  if(nsol == 1) {
    plot <- FALSE
    message('Plot is not created: only one solution was found!')
  }
  # ----------------------------
  
  if(plot == T){  
    TAB <- as.data.frame(cbind(f.solutions, count.sol))
    if(nsol > 8){    # Case A: 9 or more solutions --------------------
      global <- TAB[1,]; locals <- TAB[-1,]        
      n.locals <- dim(locals)[1]
      seven <- locals[order(locals[,2]),][-c(1:(n.locals-7)),]
      seven.ordered <- seven[order(seven[,1]),]
      other.locals <- c(NA, sum(locals[order(locals[,2]),][c(1:(n.locals-7)),2]))
      TAB2 <- rbind(global, seven.ordered, other.locals)
      p <- barplot(TAB2[,2]/K, ylab = 'Proportions', 
                   xlab = "Rotation Criterion Values",
                   names.arg = c(as.vector(TAB2[-9,1]), 'Others'), 
                   cex.names = cex, col = c('grey50',rep(NA,7),'gray85'))
      title(main=paste('Local Solutions (',K,' random starts, ',
                       nsol,' solutions)',sep=''), cex.main = .8)
    } else {    # Case B: 8 or less solutions: No 'others' ------------
      TAB2 <- TAB; n.bars <- nsol
      p <- barplot(TAB2[,2]/K, ylab = 'Proportions', 
                   xlab = "Rotation Criterion Values",
                   names.arg = c(as.vector(TAB2[,1])), 
                   cex.names = cex, col = c('grey50',rep(NA,n.bars-1)))
      title(main = paste('Local Solutions (',K,' random starts, ', 
                         nsol,' solutions)',sep=''), cex.main = .8)
    } 
  } # if(plot == T)
  
  # Output ==========================================================
  list(RotationCriterion = f.solutions, 
       N.Solutions = nsol,
       Frequencies = count.sol,
       loadings = LAMBDAS,
       Phi = PHIS)
} # MultRandRotation() -------------------------------------------------

