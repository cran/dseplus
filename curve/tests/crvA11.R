# Tests of DSE curvature functions from dsecurvature.function.testsA
if(!require("dse2"))  stop("this test requires dse2.")
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
  

# from user guide

  ARMAmodel1<-ARMA(A=array(c(1,.5,.3,0,.2,.1,0,.2,.05,1,.5,.3),c(3,2,2)),
             B=array(c(1,.2,0,.1,0,0,1,.3),c(2,2,2)), C=NULL) 


  ARMAmodel1<-l(ARMAmodel1,simulate(ARMAmodel1, rng=test.rng))
  SSmodel  <- l(toSS(ARMAmodel1),  ARMAmodel1$data)
  ARMAmodel<- l(toARMA(SSmodel), ARMAmodel1$data)
  ARMAmodel.fixed <- l(fixConstants(ARMAmodel), ARMAmodel1$data)
  
cat("DSE curvature test A 11a..")
  #  following test values have all been set using 
  #    R0.63.3pre and gnu f77 on SunOS 5.6 (Solaris)
  #  and reset as of R 1.0.0 
  hess <- hessian(ARMAmodel1)  

# 1472.1174
   good <- 1397.19588043396584
   tst  <- sum(hess)
   error <- max(abs(good-tst))
   cat("max. error ", max(error))
     
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.very.large < error) 
     {printTestValue(c(tst), digits=18); all.ok <- F }

cat("DSE curvature test A 11b..")
  hess <- hessian(SSmodel)
  #1393.953399146646
   good <- 1351.52741934382425
   tst  <- sum(hess)
   error <- max(abs(good-tst))
   cat("max. error ", max(error))
     
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.very.large < error) 
     {printTestValue(c(tst), digits=18); all.ok <- F }
 
cat("DSE curvature test A 11c..")
  hess <- hessian(ARMAmodel) 
  # R Linux      -1409.96904582771185
  # R pre 1.0 Solaris   385220.9412876417
  # R 1.0.0   Solaris     -352.3883095147531
  # Splus          728.89769711477993
  #all.ok <- test(sum(diag(hess)), 385220.9412876417, all.ok, flag="11c", fuzz=fuzz.very.large,
  #               print.values=print.values)
  warning("Skipping 11c comparison. Problem is too ill-conditioned.")
  # print(sum(diag(hess)), digits=18)
 
cat("DSE curvature test A 11d..")
  hess <- hessian(ARMAmodel.fixed)
#  3039.23996170283772   using svd until R 1.6.1 (both Linux and Soaris, but large tolerance)
#  3039.22741633394708  # R 1.6.1 La.svd with Solaris
#  3039.24495763039704  # R 1.6.1 La.svd with Mandrake 9.0

good <- if (Sys.info()[["sysname"]] == "Linux")
               3039.24495763039704 else 3039.22741633394708

   tst  <- sum(hess)
   error <- max(abs(good-tst))
   cat("max. error ", max(error))
     
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.very.large < error) 
     {printTestValue(c(tst), digits=18); all.ok <- F }

  if (! all.ok) stop("some tests FAILED")
