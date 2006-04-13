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
  print(hess <- hessian(ARMAmodel1), digits=18)  

# 1472.1174
#   good <- 1397.19588043396584
# DIFFERENCE FROM ABOVE FOUND MARCH 2006
#  WHEN SEPARATING OUT numDeriv. THERE WAS AN ERROR.
# HESSIAN USED FIRST COLUMNS IN genD RATHER THAN SKIPPING THEM. 
   good <- 480.348751175030429
   
   printTestValue(tst  <- sum(hess), digits=18)
   error <- max(abs(good-tst))
   cat("max. error ", max(error), "\n")
     
   if (any(is.na(error)) || any(is.nan(error)) || 
                          fuzz.very.large < error)  all.ok <- FALSE 

cat("DSE curvature test A 11b..")
  print(hess <- hessian(SSmodel), digits=18)
  #1393.953399146646
  # good <- 1351.52741934382425
  # DIFFERENCE FROM ABOVE FOUND MARCH 2006 (see above).
   good <-  591.625298508878586  
   printTestValue(tst  <- sum(hess), digits=18)
   error <- max(abs(good-tst))
   cat("max. error ", max(error), "\n")
     
   if (any(is.na(error)) || any(is.nan(error)) || 
                          fuzz.very.large < error)  all.ok <- FALSE
 
cat("DSE curvature test A 11c..")
  print(hess <- hessian(ARMAmodel), digits=18) 
  # R Linux      -1409.96904582771185
  # R pre 1.0 Solaris   385220.9412876417
  # R 1.0.0   Solaris     -352.3883095147531
  # Splus          728.89769711477993
  #all.ok <- test(sum(diag(hess)), 385220.9412876417, all.ok, flag="11c", fuzz=fuzz.very.large,
  #               print.values=print.values)
  warning("Skipping 11c comparison. Problem is too ill-conditioned.")
  print(sum(diag(hess)), digits=18)
 
cat("DSE curvature test A 11d..")
  print(hess <- hessian(ARMAmodel.fixed), digits=18)
#  3039.23996170283772   using svd until R 1.6.1 (both Linux and Soaris, but large tolerance)
#  3039.22741633394708  # R 1.6.1 La.svd with Solaris
#  3039.24495763039704  # R 1.6.1 La.svd with Mandrake 9.0

good <- if (Sys.info()[["sysname"]] == "Linux")
               3039.24495763039704 else 3039.22741633394708
#  above works with R-2.2.0 on RH and Gentoo, but Gentoo with ACML BLAS
#    gives 3039.23450834219602 so error in test is /10 as of Dec 2005.
  # DIFFERENCE FROM ABOVE FOUND MARCH 2006 (see above).
   good <-  805.34949654355762
   printTestValue(tst  <- sum(hess), digits=18)
   error <- max(abs(good-tst))
   cat("max. error ", max(error), "\n")
     
   if (any(is.na(error)) || any(is.nan(error)) ||
                   fuzz.very.large < (error/10)) all.ok <- FALSE

if (! all.ok) stop("some tests FAILED")
