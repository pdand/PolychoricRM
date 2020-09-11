# efaMR: Exploratory factor analysis with multiple rotations


# ===================================================================================
efaMR <- function(x=NULL, factors=NULL, covmat=NULL, n.obs=NULL, 
                  dist='normal', fm='ols', rtype='oblique', rotation = 'CF-varimax', 
                  input.A=NULL, additionalRC = NULL, 
                  nstart = 100, compare = 'First', plot = T, cex = .5,
                  normalize=FALSE, geomin.delta=.01, 
                  MTarget=NULL, MWeight=NULL, PhiTarget = NULL, PhiWeight = NULL, 
                  useorder=FALSE, mnames=NULL, fnames=NULL, wxt2 = 1) {
  
  # 1. check input arguments -----------------------------------------------------------
  if ( (is.null(x)) & (is.null(covmat)) & (is.null(input.A)) ) stop ("Neither raw data, the correlation matrix, nor the unrotated factor loading matrix is provided!")
  if (is.null(input.A)){
    if (! (dist == 'normal') & (is.null(x)) ) stop ("Raw data are required for non-normal distributions!")
    if ( (is.null(x)) & (is.null(n.obs)) ) stop ("The sample size is not provided for the correlation matrix!")
  } else message ("Factor rotation is conducted with the input unrotated factor loading matrix; factor extraction is unnecessary.")
  
  if (!(is.null(x))) {
    p = ncol (x)
    n = nrow (x)
    if (n <p) stop ("The sample size is less than the number of manifest variables!")
    
    if (dist=='ordinal') {
      polychor= PolychoricRM(x)
      R0 = polychor$correlation
    } else {
      R0 = cor(x)
    }
    
    if (!(is.null(covmat))) {
      if(min(abs(R0 - covmat)) > 0.0001) message ('covmat is different from the one computed from the raw data! The one computed from raw data is used for EFA!') 
    } # if (!(is.null(covmat)))
  } else {            # (!(is.null(x)))
    
    if(is.null(input.A)){
      p = ncol(covmat)
      n = n.obs
      
      mvariance <- diag(covmat)
      
      if (all(mvariance==1)) {
        R0 = covmat
      } else {
        message('covmat is not a correlation matrix; EFA is conducted with the corresponding correlation matrix.')
        R0 = covmat
        msd = sqrt(mvariance)
        for (i in 1:p) {
          R0[i,1:p] = R0[i,1:p] / msd[i]
          R0[1:p,i] = R0[1:p,i] / msd[i] 
        }
      } # if the input matrix is a covariance matrix 
    }# if(is.null(input.A))
  } # is.null(x) 
  
  if(!is.null(input.A)){
    p <- dim(input.A)[1]
    factors <- dim(input.A)[2]
    n <- n.obs
  }
  
  if(is.null(input.A)){
    max.factors = floor(((2*p + 1) - sqrt(8*p+1))/2)
    
    ev = eigen(R0)$values
    
    if (is.null(factors)) {factors = length(which(ev > 1))
    } else {
      if (factors > max.factors) stop (paste(factors, "factors is too many for",p,"variables."))
    }
  }
  
  # 2. manifest variable names and factor names -------------------------------------
  if (is.null(mnames)) {
    mnames = rep(" ", p)  
    for(i in 1:p) {
      mnames[i] = paste("MV",i,sep="")
    }
  } else{
    if(length(mnames) != p) stop("The number of MV names is different from the number of MVs!")
  }
  
  if (is.null(fnames)) {
    fnames = rep(" ", factors)  
    for(i in 1:factors) {
      fnames[i] = paste("F",i,sep="")
    }
  } else{
    if(length(fnames) != factors) stop("The number of factor names is different from the number of factors!")
  }
  
  # 3. extract m factors -------------------------------------------------------
  if(is.null(input.A)){
    A.lst = fa.extract(R0, factors, extraction = fm)
    FE.Arg <- list(factors = factors, extraction = fm)
    unrotated <- A.lst$Unrotated
  } else { 
    # if unrotated factor loading matrix is provided instead of (in addition to) x
    unrotated <- input.A
  }
  
  
  # 4. factor rotation with multiple random starts ------------------------------
  
  
  
  
  
  # 4A. nstart!=1 --- mutiple random starts for one rotation criterion ----------
  
  if(factors == 1){
    MultipleSolutions <- 'No rotation has been done for one-factor models.'
    Comparisons <- NULL
  }  else  { # if(factors!=1)
    
    multiple <- MultRandRotation(unrotated, epsilon = geomin.delta, nstart = nstart, 
                                 plot = plot, cex = cex, rotation = rotation, rtype = rtype, 
                                 normalize = normalize, MWeight = MWeight, MTarget = MTarget,
                                 wxt2 = wxt2)
    
    NumberSolutions <- length(multiple$loadings)
    Solutions <- list()
    for(i in 1:NumberSolutions){
      Solutions[[i]] <- list(Lambda = multiple$loadings[[i]], Phi = multiple$Phi[[i]])
    }
    
    
    MultipleSolutions <- list(nstart = nstart,
                              RotationCriterion    = multiple$RotationCriterion, 
                              NumberSolutions      = multiple$N.Solutions,
                              FrequenciesSolutions = multiple$Frequencies,
                              Solutions            = Solutions)
    
    if(multiple$N.Solutions != 1){
      comparison <- CompareSolutions(multiple$loadings, compare)
      Comparisons <- list(MinimumCongruence = comparison$MinimumCongruence,
                          RawCongruences    = comparison$RawCongruences)
    } else { 
      Comparisons <- 'No comparison with only one solution!'  
    }
  } # if(factors!=1)
  
  
  
  # 5. organizing details -------------------------------------------------------
  details = list(manifest = p, factors = factors, n.obs = n, dist = dist, 
                 fm = fm, rtype = rtype, rotation = rotation, 
                 normalize = normalize, geomin.delta = geomin.delta, 
                 MTarget = MTarget, MWeight = MWeight)
  
  # 6. organizing outputs ------------------------------------------------------- 
  if(is.null(input.A)){
    fdiscrepancy = A.lst$f
    convergence = A.lst$convergence
  } else {
    fdiscrepancy <- NULL
    convergence <- NULL
  }
  
  dimnames(unrotated) = list(mnames,fnames)
  
  # 7. output ----------------------------------------------------------------  
  efaout = list(details = details, unrotated = round(unrotated,3), 
                fdiscrepancy = fdiscrepancy, convergence = convergence,
                MultipleSolutions = MultipleSolutions, Comparisons = Comparisons)
  
  class(efaout) <- "efaMRS" 
  
  return(efaout)
} # efaMRS ==================================================================
# ==========================================================================


