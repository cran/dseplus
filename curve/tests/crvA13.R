# Tests of DSE curvature functions from dsecurvature.function.testsA
if(!require("dse2"))  stop("this test requires dse2.")
if(!require("curve"))stop("this test requires curve.")
 Sys.info()
 DSEversion()
 
fuzz.small <- 1e-12
fuzz.large <- 1e-6
fuzz.very.large <- 1e-2
digits <- 18
all.ok <- TRUE
test.rng <- list(kind="Wichmann-Hill",seed=c(979,1479,1542),normal.kind="Box-Muller")


# comparison values come only from a previous run of the 
#  code (theoretical values would be nice)...
# Test values have been changed with change to RNG when R 1.0.0 was released
#   (Feb. 29, 2000) and also previously.
  

# from user guide

  ARMAmodel1<-ARMA(A=array(c(1,.5,.3,0,.2,.1,0,.2,.05,1,.5,.3),c(3,2,2)),
             B=array(c(1,.2,0,.1,0,0,1,.3),c(2,2,2)), C=NULL) 

 ARMAmodel1<-l(ARMAmodel1,simulate(ARMAmodel1, rng=test.rng))
  SSmodel  <- l(toSS(ARMAmodel1),  ARMAmodel1$data)
 

cat("DSE curvature test A 13...\n")
  curvatureSS <- curvature(SSmodel, warn=FALSE)$stats  
  # neg sqrt in axis ratio produces warning if warn=T
#  good <- c(16, 200, 0.05, 1.492400232033412, 1.197222853812868,
#       1.945129668075371,   1.560408288784785,   1.000000000054351)
   good <- c(16, 200, 0.05, 1.42377065952031878,  1.11780582522466698,
        1.85568086289751188,  1.45689958170905332,  1.00000000346571949)
   tst  <- curvatureSS[-9]
   error <- max(abs(good-tst))
   cat("max. error ", max(error))
     
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.large < error) 
     {printTestValue(c(tst), digits=18); all.ok <- F }

  if (! all.ok) stop("some tests FAILED")
