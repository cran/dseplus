# gctorture()

# Tests of DSE curvature functions from dsecurvature.function.testsA
if(!require("dse2"))  stop("this test requires dse2.")
if(!require("curve"))stop("this test requires curve.")
 Sys.info()
 DSEversion()
 
fuzz.large <- 1e-8
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


# previous rng
#  good <- c(4, 100, 0.05,  0.9992455609137233,  0.9470084226509906, 
#     1.56931709092204,  1.487278563994181,  1.000000000828572, 
#     1.321413857107741)
 
  good <- c(4, 100, 0.05, 1.02548364377805989,  0.979897817007443384,  
     1.61052405093533291,  1.53893142160253493,  1.00000000114714904)  #, NaN


cat("DSE curvature test A 6b...\n")
# gctorture()

  func.residual <- function(coefficients,Shape,data)
   {c(l(setArrays(Shape, coefficients=coefficients),data,result="pred")
      - outputData(data))} 

  curvatureSS.def <- curvature(func.residual, coef(SSmodel), 
               func.args=list(Shape=TSmodel(SSmodel), data=TSdata(SSmodel)),
                     d=0.01, eps=1e-4,r=6, show.details=FALSE, warn=FALSE)$stats
  # neg sqrt in axis ratio produces warning if warn=T

   tst  <-  curvatureSS.def[-9]
#   good <-  curvatureSS[-9]
   error <- max(abs(good - tst))
   cat("max. error ", max(error), "\n")
   
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.large < error) 
     {if (any(is.na(error)))  cat("contains na's: ",  is.na(error),  "\n")
      if (any(is.nan(error))) cat("contains nan's: ", is.nan(error), "\n")
      printTestValue(c(tst), digits=18)
      all.ok <- F  
     }

cat("DSE curvature test A 6a...\n")
  curvatureSS <- curvature(SSmodel, warn=FALSE)$stats
  # neg sqrt in axis ratio produces warning if warn=T
   tst  <-  curvatureSS[-9]

   error <- max(abs(good - tst))
   cat("max. error ", error, "\n")
   
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.large < error) 
     {if (any(is.na(error)))  cat("contains na's: ",  is.na(error),  "\n")
      if (any(is.nan(error))) cat("contains nan's: ", is.nan(error), "\n")
      printTestValue(c(tst), digits=18)
      all.ok <- F  
     }

  if (! all.ok) stop("some tests FAILED")
