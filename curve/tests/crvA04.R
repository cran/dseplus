 # Tests of DSE curvature functions from dsecurvature.function.testsA
 require("dse2"); require("curve") #,  warn.conflicts=FALSE)
 Sys.info()
 version.dse()
 
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
  

# simplified from user guide
  z<-ARMA(A=array(c(1,.5,.3),c(3,1,1)),
          B=array(1,c(1,1,1)),
          C=NULL, description="simplified guide example")
  VARmodel <-l(z,simulate(z, rng=test.rng))

  SSmodel  <- l(to.SS(VARmodel),  VARmodel$data)
  ARMAmodel<- l(to.ARMA(SSmodel), VARmodel$data)

cat("DSE curvature test A 4a...")

# calculating from the score using genD might be quicker.
  hessianVAR <- hessian(VARmodel) 

   tst  <- hessianVAR
   good <- c( 9.83444900505639374,  -5.64196107342148512,
           -5.64196107342148512,  116.345053001385807)

   error <- max(abs(good-tst))
   cat("max. error ", max(error))
   
# relaxed from fuzz.large to 10*fuzz.large for R 1.3.0 in Linux
   if (any(is.na(error)) || any(is.nan(error)) || 10*fuzz.large < error) 
     {print.test.value(c(tst), digits=18)
      all.ok <- F  
     }

cat("DSE curvature test A 4b...")

  hessianSS <- hessian(SSmodel)
  tst  <- hessianSS

  good <-  c( 3.95368029144797717,  -0.418007621718194833,  -9.41644138321294832,
            36.9034481180689653,  -0.418007621718194833,  1.68828078196901132,  
	    39.9657123856431085,  -16.2517115170908539,  -9.41644138321294832,  
	    39.9657123856431085,  10.2620303429070088,  -65.7865038150653731,  
	    36.9034481180689653,  -16.2517115170908539,  -65.7865038150653731,  
	   107.778335544667868)

   error <- max(abs(good-tst))
   cat("max. error ", error)
   
   if (any(is.na(error)) || any(is.nan(error)) || 10*fuzz.large < error) 
     {print.test.value(c(tst), digits=18)
      all.ok <- F  
     }

cat("DSE curvature test A 4c...")

  hessianARMA <- hessian(ARMAmodel)

   tst  <- hessianARMA
   good <- hessianVAR
   error <- max(abs(good-tst))
   cat("max. error ", error)
   
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {print.test.value(c(tst), digits=18)
      all.ok <- F  
     }

  if (! all.ok) stop("some tests FAILED")
