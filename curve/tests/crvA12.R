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
  

# from user guide

  VARmodel<-ARMA(A=array(c(1,.5,.3,0,.2,.1,0,.2,.05,1,.5,.3),c(3,2,2)),
             B=array(c(1,.2,0,.1,0,0,1,.3),c(2,2,2)), C=NULL) 

# Note this gives a terrible fit.
  VARmodel<-l(VARmodel,simulate(VARmodel, rng=test.rng))
  SSmodel  <- l(to.SS(VARmodel),  VARmodel$data)
  ARMAmodel<- l(to.ARMA(SSmodel), VARmodel$data)


cat("DSE curvature test A 12a..\n")
  curvatureVAR <- curvature(VARmodel,warn=FALSE)$stats
#  good <- c(11, 200, 0.05, 0.8490449316698463, 0.712275843318316,
#        1.151573691713139,  0.9660715137831059, 1.000000000576390) #, NaN)


  good <- c(11, 200, 0.05, 0.807396116175452816,  0.681455079712046774,
         1.0950847140096498,  0.924268802049481475,  1.00000000053258109)

  good <- c(11, 200, 0.05, 0.807400747659583362,  0.681458988757835171,
         1.09509099576822266,  0.924274103952975268,  1.00000000053258109)

   tst  <- curvatureVAR[-9]
   error <- max(abs(good-tst))
   cat("max. error ", max(error))
     
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.large < error) 
     {print.test.value(c(tst), digits=18); all.ok <- F }

cat("DSE curvature test A 12b..\n")
  func.residual <- function(coefficients,Shape,data)
   {c(l(set.arrays(Shape,coefficients=coefficients),data,result="pred")
       - output.data(data))}
   
  curvatureVAR.def <- curvature.default(func.residual, coef(VARmodel), 
              func.args=list(Shape=TSmodel(VARmodel), data=TSdata(VARmodel)),
                     d=0.01, eps=1e-4,r=6, show.details=FALSE)$stats

   good <- curvatureVAR[-9]
   tst  <- curvatureVAR.def[-9]
   error <- max(abs(good-tst))
   cat("max. error ", max(error))
     
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {print.test.value(c(tst), digits=18); all.ok <- F }

  if (! all.ok) stop("some tests FAILED")
