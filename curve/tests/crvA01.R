# Tests of DSE curvature functions from dsecurvature.function.testsA
if(!require("dse2"))  stop("this test requires dse2.")
if(!require("curve"))stop("this test requires curve.")
 Sys.info()
 version.dse()

fuzz.small <- 1e-12
fuzz.large <- 1e-6
fuzz.very.large <- 1e-2
digits <- 18
all.ok <- T  

  random.number.test()

test.rng <- list(kind="Wichmann-Hill",seed=c(979,1479,1542),normal.kind="Box-Muller")


# comparison values come only from a previous run of the 
#  code (theoretical values would be nice)...
# Test values have been changed with change to RNG when R 1.0.0 was released
#   (Feb. 29, 2000) and also previously.



cat("DSE curvature test A 1 ...")


# simplified from user guide
  z <- ARMA(A=array(c(1,.5,.3),c(3,1,1)),
            B=array(1,c(1,1,1)),
            C=NULL, description="simplified guide example")
  VARmodel <-l(z,simulate(z, rng=test.rng))

# unstable model
  z<-ARMA(A=array(c(1,-0.5,-0.5),c(3,1,1)),
          B=array(1,c(1,1,1)),
          C=NULL, description="simplified guide example")
  VARmodel2 <-l(z,simulate(z, rng=test.rng))
# Mod(roots(VARmodel, by.poly=T))
# [1] 0.5477226 0.5477226
# Mod(roots(VARmodel2, by.poly=T))
# [1] 1.0 0.5

  if (fuzz.small < abs(VARmodel$estimates$like[1] - 143.78939695547763 ))
    {warning("VARmodel  likelihood  does not correspond to expected value.")
     cat("VARmodel  likelihood:")
     print(VARmodel$estimates$like[1], digits=digits)
    }
 
  if( fuzz.small < abs(sum(VARmodel$data$noise$w) - 9.384686064093962  ))
    {warning("Check sum of noise does not correspond to expected test value.")
     cat("Check sum of noise:")
     print(sum(VARmodel$data$noise$w), digits=digits)
    }

  SSmodel  <- l(to.SS(VARmodel),  VARmodel$data)
  ARMAmodel<- l(to.ARMA(SSmodel), VARmodel$data)
  if (fuzz.small < abs(ARMAmodel$estimates$like[1]- VARmodel$estimates$like[1]))
    {warning("ARMAmodel likelihood does not correspond to expected test value.")
     cat("ARMAmodel  likelihood:")
     print(ARMAmodel$estimates$like[1], digits=digits)
    }

   spanVAR <- span(VARmodel, compiled=F)

   good <- spanVAR
   tst  <- c(12.583380392358416,  9.21209438442244277)
   error <- max(abs(good-tst))
   cat("max. error ", error)
   
   if (any(is.na(error)) || any(is.nan(error)) || 10*fuzz.small < error) #10* for Linux vs Solaris R
     {printTestValue(c(tst), digits=18)
      all.ok <- F  
     }
 
cat("DSE curvature test A 1b...")

  spanVAR.f <- span(VARmodel, compiled=.DSECOMPILED)

   good <- spanVAR
   tst  <- spanVAR.f
   error <- max(abs(good-tst))
   cat("max. error ", error)

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {printTestValue(c(tst), digits=18)
      all.ok <- F  
     }


  if (! all.ok) stop("some tests FAILED")
