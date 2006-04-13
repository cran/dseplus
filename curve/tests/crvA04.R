 # Tests of DSE curvature functions from dsecurvature.function.testsA
if(!require("dse1"))  stop("this test requires dse1.")
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
  

# simplified from user guide
  z<-ARMA(A=array(c(1,.5,.3),c(3,1,1)),
          B=array(1,c(1,1,1)),
          C=NULL, description="simplified guide example")
  VARmodel <-l(z,simulate(z, rng=test.rng))

  SSmodel  <- l(toSS(VARmodel),  VARmodel$data)
  ARMAmodel<- l(toARMA(SSmodel), VARmodel$data)

cat("DSE curvature test A 4a...")

# calculating from the score using genD might be quicker.
  printTestValue(c(hessianVAR <- hessian(VARmodel) ), digits=18)


#   good <- c( 9.83444900505639374,  -5.64196107342148512,
#           -5.64196107342148512,  116.345053001385807)
# DIFFERENCE FROM ABOVE FOUND MARCH 2006
#  WHEN SEPARATING OUT numDeriv. THERE WAS AN ERROR.
# HESSIAN USED FIRST COLUMNS IN genD RATHER THAN SKIPPING THEM. 
   good <- c(116.345053142090151,  -34.2425297946331142,  
             -34.2425297946331142, 115.240633888860629)
   
   
   error <- max(abs(good - hessianVAR))
   cat("max. error ", max(error))
   
# relaxed from fuzz.large to 10*fuzz.large for R 1.3.0 in Linux
   if (any(is.na(error)) || any(is.nan(error)) || 10*fuzz.large < error) 
     all.ok <- F  


cat("DSE curvature test A 4b...")

  printTestValue(c(hessianSS <- hessian(SSmodel)), digits=18)

#  good <-  c( 3.95368029144797717,  -0.418007621718194833,  -9.41644138321294832,
#            36.9034481180689653,  -0.418007621718194833,  1.68828078196901132,  
#	    39.9657123856431085,  -16.2517115170908539,  -9.41644138321294832,  
#	    39.9657123856431085,  10.2620303429070088,  -65.7865038150653731,  
#	    36.9034481180689653,  -16.2517115170908539,  -65.7865038150653731,  
#	   107.778335544667868)
# DIFFERENCE FROM ABOVE FOUND MARCH 2006 (see above).
 good <-  c( 
 39.9657210069456852,  10.2620249375422095,  -16.2517136853819615,  -14.065782999853182,  
 10.2620249375422095,  36.9034515031556225,  -65.7865066379252283,  -6.53008602892200951,
-16.2517136853819615, -65.7865066379252283,  107.778340252653535,   35.3477274147038472,  
-14.065782999853182,   -6.53008602892200951,  35.3477274147038472,  92.501773696778514)


 error <- max(abs(good - hessianSS))
 cat("max. error ", error, "\n")
   
 if (any(is.na(error))|| any(is.nan(error))||10*fuzz.large < error) all.ok <- F  

cat("DSE curvature test A 4c...")

 printTestValue(c(hessianARMA <- hessian(ARMAmodel)), digits=18)

 error <- max(abs(hessianVAR - hessianARMA))
 cat("max. error ", error, "\n")
   
 if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) all.ok <- F  

if (! all.ok) stop("some tests FAILED")
