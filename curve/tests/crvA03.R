 # Tests of DSE curvature functions from dsecurvature.function.testsA
if(!require("dse1"))  stop("this test requires dse1.")
if(!require("curve"))stop("this test requires curve.")
 Sys.info()
 DSEversion()
 
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
  

# simplified from user guide
  z<-ARMA(A=array(c(1,.5,.3),c(3,1,1)),
          B=array(1,c(1,1,1)),
          C=NULL, description="simplified guide example")
  VARmodel <-l(z,simulate(z, rng=test.rng))

  SSmodel  <- l(toSS(VARmodel),  VARmodel$data)
  ARMAmodel<- l(toARMA(SSmodel), VARmodel$data)

cat("DSE curvature test A 3a...")

  spanVAR <- span(VARmodel, compiled=.DSEflags()$COMPILED)
  spanARMA.f <- span(ARMAmodel, compiled=.DSEflags()$COMPILED)

   tst  <- spanARMA.f
   good <- spanVAR
   error <- max(abs(good-tst))
   cat("max. error ", max(error))
   
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {printTestValue(c(tst), digits=18)
      all.ok <- F  
     }

cat("DSE curvature test A 3b...")

  spanARMA <- span(ARMAmodel, compiled=F)

   tst  <- spanARMA
   good <- spanARMA.f
   error <- max(abs(good-tst))
   cat("max. error ", max(error))
   
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {printTestValue(c(tst), digits=18)
      all.ok <- F  
     }

  if (! all.ok) stop("some tests FAILED")
