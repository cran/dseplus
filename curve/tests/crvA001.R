 # Tests of DSE curvature functions 
if(!require("dse2"))  stop("this test requires dse2.")
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
  z  <-  ARMA(A=array(c(1,.5,.3),c(3,1,1)),
              B=array(1,c(1,1,1)),
              C=NULL, description="simplified guide example")

  VARmodel <-l(z,simulate(z, rng=test.rng))
  SSmodel  <- l(toSS(VARmodel),  VARmodel$data)

   spanSS.notcompiled <- span(SSmodel, compiled=F)
   cat("length spanSS.notcompiled ", length(spanSS.notcompiled), "\n")
   spanSS.compiled    <- span(SSmodel, compiled=T)
   cat("length spanSS.compiled ", length(spanSS.compiled), "\n")

   error <- max(abs(spanSS.compiled - spanSS.notcompiled))
   cat("max. error ", error, "\n")
   
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {printTestValue(c(tst), digits=18)
      all.ok <- F  
     }


  if (! all.ok) stop("some tests FAILED")
