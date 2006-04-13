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

 cat("DSE curvature test A 5a...\n")
   curvatureVAR <- curvature(VARmodel)$stats

   tst  <-  curvatureVAR
   good <-  c(2, 100, 0.05, 0, 0, 0, 0, 1, 1)
   error <- max(abs(good-tst))
   cat("max. error ", error)
   
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.large < error) 
     {printTestValue(c(tst), digits=18)
      all.ok <- F  
     }

 cat("DSE curvature test A 5b...\n")

# and for comparison
        
#  curvatureVAR.def <- curvature(function(coefficients, Shape, data)
#   {c(l(setArrays(Shape,coefficients=coefficients),data,result="pred")
#        - outputData(data))}, coef(VARmodel), 
#     method="Richardson", method.args=list(d=0.01, eps=1e-4, r=6, v=2),
#     show.details=FALSE,
#     Shape=TSmodel(VARmodel),data=TSdata(VARmodel))$stats

  curvatureVAR.def <- curvature(
       function(coefficients, Shape, data)
          {c(l(setArrays(Shape,coefficients=coefficients),
	       data,result="pred") - outputData(data))}, 
       coef(VARmodel), 
       method="Richardson",
       method.args=list(d=0.01, eps=1e-4,r=6, v=2, show.details=FALSE),
       Shape=TSmodel(VARmodel),
       data=TSdata(VARmodel)
       )$stats


   tst  <-  curvatureVAR.def
   good <-  curvatureVAR
   error <- max(abs(good-tst))
   cat("max. error ", error)
   
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.large < error) 
     {printTestValue(c(tst), digits=18)
      all.ok <- F  
     }

  if (! all.ok) stop("some tests FAILED")
