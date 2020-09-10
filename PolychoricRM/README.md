# EFAutilities

The goal of EFAutilities is to provide a number of utilities functions for exploratory factor analysis.  In particular, it computes standard errors for model parameters under a variety of conditions. 



This package can be installed directly from CRAN:

install.packages("EFAutilities")
library(EFAutilities)



## Example

Examples using the data sets included in the packages:

data("CPAI537")                   # Chinese personality assessment inventory (N = 537)

efa(x=CPAI537,factors=4, fm='ml') # normal, ml, oblique, CF-varimax, information, merror='NO'
```R
...
```
