# Tests of DSE curvature functions from dsecurvature.function.testsA
 require("dse2"); require("curve") #,  warn.conflicts=F)
 Sys.info()
 version.dse()
 
fuzz.small <- 1e-12
fuzz.large <- 1e-6
fuzz.very.large <- 1e-2
digits <- 18
all.ok <- T
test.rng <- list(kind="Wichmann-Hill",seed=c(979,1479,1542),normal.kind="Box-Muller")


# comparison values come only from a previous run of the 
#  code (theoretical values would be nice)...
# Test values have been changed with change to RNG when R 1.0.0 was released
#   (Feb. 29, 2000) and also previously.
  

# from user guide

  VARmodel<-ARMA(A=array(c(1,.5,.3,0,.2,.1,0,.2,.05,1,.5,.3),c(3,2,2)),
             B=array(c(1,.2,0,.1,0,0,1,.3),c(2,2,2)), C=NULL) 

# Note this gives a terrible fit.
  VARmodel<-l(VARmodel,simulate(VARmodel, rng=test.rng))
  SSmodel  <- l(to.SS(VARmodel),  VARmodel$data)
  ARMAmodel<- l(to.ARMA(SSmodel), VARmodel$data)

#  As of 2001.3 to.ARMA was changed to call fix.constants and eliminate some
#  near zero parameter. These were not only the cause of considerable numerical
#  instability in the curvature calculation, but also gave a different number of
#  parameters (20 vs 22) in R on Solaris and Linux. Some of the parameter value 
#  are about 1e-18 and were treated as zero in Linux but not Solaris. The result
#  now is 18 parameters, much improve stability, and more reasonable values.
 
good <- c(18, 200,  0.05,  2.34112169445768981,   2.29032831656791869,
  3.01704740785875769,  2.95158902973962167,  1.00047805155781822) #, NaN)

# above is R. S gives (within 1e-5 from R):
#        c(18, 200,   0.05, 2.3411215346688303285, 2.2903275022515097170,
#  3.0170472019358562932,  2.9515879803147995197,   1.0004779939595003224 )
  
 cat("DSE curvature test A 15...\n")
 
  curvatureARMA <- curvature(ARMAmodel, warn=F)$stats 
      
   tst  <- curvatureARMA [-9]
   error <- max(abs(good - tst))
   cat("max. error ", max(error))
     
   if (any(is.na(error)) || any(is.nan(error)) || 1e-5 < error) 
     {print.test.value(c(tst), digits=18); all.ok <- F }

  if (! all.ok) stop("some tests FAILED")
