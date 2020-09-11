# TSCovCorr computes the asymptotic covariance matrix of lagged correlations
# X, input, n by p
# Internal functions: find.r is to find lagged correlations of log 0 and negative lages.
# Note that R does not allow 0 or negative indices in its arrays.

TSCovCorr <- function(X,nlagin=20) {


#--------------------------------------------------------------------

find.r <- function(dumlag, dumrow, dumcol) {

if (dumlag < 0) {
truelag = dumlag * (-1) + 1
dumtemp = dumrow
dumrow = dumcol
dumcol = dumtemp
} else {
truelag = dumlag + 1
}

Gamma[truelag,dumrow,dumcol]

} # find.r

#------------------------------------------------------------------


# declare variables

np = ncol(X)
nlagout = 0
nlagout1 = nlagout + 1


# Calculate lagged correlations
Gamma = acf(X,lag.max=nlagin,plot = FALSE)$acf


CovCorr = array( array(0, (nlagout1*np*np)**2), dim=c(nlagout1,np,np,nlagout1,np,np))

# Declare the logical array LogCovCorr and set its initial values
LogCovCorr <- CovCorr > 5
for (i in 1:np) {
LogCovCorr[1,i,i,,,] = TRUE
LogCovCorr[,,,1,i,i] = TRUE
}


nsumterm = nlagin - nlagout

for (m in 0:nlagout) {
 for (j in 1:np) {
    for (i in 1:np) {
      for (n in 0:nlagout) {
        for (l in 1:np) {
          for (k in 1:np) {
           m1 = m + 1
           n1 = n + 1

        if (! LogCovCorr[m1,i,j,n1,k,l]) {
          
          for (u in -nsumterm:nsumterm) {
          r.m.ij = Gamma[m1,i,j]
          r.n.kl = Gamma[n1,k,l]
                
          CovCorr[m1,i,j,n1,k,l] = CovCorr[m1,i,j,n1,k,l] + 
          0.5 * r.m.ij * r.n.kl * ((find.r(u,i,k))**2 + (find.r(u,j,k))**2 + (find.r(u,i,l))**2 + (find.r(u,j,l))**2) - # terms 1, 2, 3, 4
          r.n.kl * (find.r(u,j,k) * find.r(u+m,i,k)+find.r(u,j,l)*find.r(u+m,i,l)) - # terms 5 and 6
          r.m.ij * (find.r(u,i,l) * find.r(u-n,i,k)+find.r(u,j,l) * find.r(u-n,j,k)) + # terms 7 and 8
          find.r(u,j,l) * find.r(u-n+m,i,k) + find.r(u-n,j,k) * find.r(u+m,i,l)
          } # u

          CovCorr[n1,k,l,m1,i,j] = CovCorr[m1,i,j,n1,k,l]
         
          LogCovCorr[m1,i,j,n1,k,l] = TRUE
          LogCovCorr[n1,k,l,m1,i,j] = TRUE

        
        } # (! LogCovCorr[m,i,j,n,k,l])


      } # k
     } # l
    } # n
   } # i
  } # j
} # m


temp = array(CovCorr[1,,,1,,])
asc = array(temp,dim=c(np*np, np*np))


list(corr = Gamma[1,,], asc = asc)



} # TSCovCorr



