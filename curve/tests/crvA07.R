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

  curvatureVAR <- curvature(VARmodel)$stats

cat("DSE curvature test A 7a...\n")
  curvatureARMA <- curvature(ARMAmodel)$stats

   tst  <-  curvatureARMA
   good <-  curvatureVAR
   error <- max(abs(good-tst))
   cat("max. error ", max(error))
   
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {print.test.value(c(tst), digits=18)
      all.ok <- F  
     }

cat("DSE curvature test A 7b...\n")
 
  func.residual <- function(coefficients, Shape, data)
   {c(l(set.arrays(Shape,coefficients=coefficients), data, result="pred")
       - output.data(data))}

  curvatureARMA.def <- curvature.default(func.residual, coef(ARMAmodel), 
               func.args=list(Shape=TSmodel(ARMAmodel), data=TSdata(ARMAmodel)),
                     d=0.01, eps=1e-4,r=6, show.details=FALSE)$stats

   tst  <-  curvatureARMA.def
   good <-  curvatureARMA
   error <- max(abs(good-tst))
   cat("max. error ", max(error))
   
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.large < error) 
     {print.test.value(c(tst), digits=18)
      all.ok <- F  
     }

  if (! all.ok) stop("some tests FAILED")
