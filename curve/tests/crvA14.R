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
  ARMAmodel<- l(toARMA(SSmodel), ARMAmodel1$data)
  ARMAmodel.fixed <- l(fixConstants(ARMAmodel), ARMAmodel1$data)
  

cat("DSE curvature test A 14...\n")
  curvatureARMA.fixed <- curvature(ARMAmodel.fixed, warn=FALSE)$stats
  # neg sqrt in axis ratio produces warning if warn=T
# R 0.64.1
#  good <- c(18, 200, 0.05, 2.381219927481035, 2.270893359206098,
#            3.06872275232684,  2.92654283591333,  1.001605964896068)
# Splus 3.3
# good <- c(18, 200, 0.05, 2.381219668018957, 2.270892355686405,
#           3.068722417953193, 2.926541542658688, 1.00160610746321)
# R 1.1.0 (devel)
  good <- c(18, 200, 0.05, 2.34112169445768892,  2.29032831656791913,
            3.0170474078587568,  2.95158902973962212,  1.00047805155781822)

   tst  <- curvatureARMA.fixed[-9]
   error <- max(abs(good-tst))
   cat("max. error ", max(error))
     
   if (any(is.na(error)) || any(is.nan(error)) || 10*fuzz.large < error) 
     {printTestValue(c(tst), digits=18); all.ok <- F }

  if (! all.ok) stop("some tests FAILED")
