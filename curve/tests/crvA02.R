 # Tests of DSE curvature functions from dsecurvature.function.testsA
 require("dse2"); require("curve") #,  warn.conflicts=F)
 Sys.info()
 version.dse()
 
fuzz.small <- 1e-11
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
  SSmodel  <- l(to.SS(VARmodel),  VARmodel$data)

cat("DSE curvature test A 2a...\n")

   good <- c( 11.5173146571970904,   9.8419970117338984,
               5.42943770603545595,  1.65151435950942771)
   spanSS <- span(SSmodel, compiled=F)
   error <- max(abs(good - spanSS))
   cat("max. error ", max(error), "\n")
   
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {print.test.value(c(spanSS), digits=18)
      all.ok <- F  
     }


cat("DSE curvature test A 2b...\n")

  spanSS.f <- span(SSmodel, compiled=.DSECOMPILED)
  error <- max(abs(good - spanSS.f))
  cat("max. error ", max(error), "\n")
   
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.large < error) 
     {print.test.value(c(spanSS.f), digits=18)
      all.ok <- F  
     }


  if (! all.ok) stop("some tests FAILED")

