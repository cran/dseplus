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
  ARMAmodel.fixed <- l(fix.constants(ARMAmodel), VARmodel$data)
  
cat("DSE curvature test A 11a..")
  #  following test values have all been set using 
  #    R0.63.3pre and gnu f77 on SunOS 5.6 (Solaris)
  #  and reset as of R 1.0.0 
  hess <- hessian(VARmodel)  

# 1472.1174
   good <- 1397.19588043396584
   tst  <- sum(hess)
   error <- max(abs(good-tst))
   cat("max. error ", max(error))
     
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.very.large < error) 
     {print.test.value(c(tst), digits=18); all.ok <- F }

cat("DSE curvature test A 11b..")
  hess <- hessian(SSmodel)
  #1393.953399146646
   good <- 1351.52741934382425
   tst  <- sum(hess)
   error <- max(abs(good-tst))
   cat("max. error ", max(error))
     
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.very.large < error) 
     {print.test.value(c(tst), digits=18); all.ok <- F }
 
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
  # 3038.56

   good <- 3039.23996170283772
   tst  <- sum(hess)
   error <- max(abs(good-tst))
   cat("max. error ", max(error))
     
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.very.large < error) 
     {print.test.value(c(tst), digits=18); all.ok <- F }

  if (! all.ok) stop("some tests FAILED")
