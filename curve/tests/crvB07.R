# Tests of DSE curvature functions (previously dsecurvature.function.testsB )

# comparison values come only from a previous run of the 
#  code (theoretical values would be nice)...
# Test values have been changed with change to to.ARMA in 2001.2 which
# eliminates near zero parameter values using fix.constants. The result is
# much more stable and believable curvature results. The span results do not
# change much (as would be hoped) but do change more than the tolerance of 
# these tests. Old values in comments are  strictly for historical reference.

 require("dse2"); require("curve") #,  warn.conflicts=F)
 Sys.info()
 version.dse()

fuzz.small <- 1e-12
fuzz.large <- 1e-6
digits <- 18
all.ok <- T  

  if (is.R()) data("eg1.DSE.data.diff", package="dse1") else 
  if (is.S()) source(paste(DSE.HOME, "/data/eg1.DSE.data.diff.R", sep=""))

# data size affects memory constraints
  data <- eg1.DSE.data.diff
   input.data(data) <- NULL
  output.data(data) <- output.data(data)[1:50,1:2]

  VARmodel <- est.VARX.ls(data, re.add.means=F)
  SSmodel  <- l(to.SS(VARmodel),  data)
  ARMAmodel<- l(to.ARMA(SSmodel), data)


cat("DSE curvature test B 7 ...")

  hessianVAR <- hessian(VARmodel)
  
#  good <- if(is.Splus()) 7219.717083137912   else if(is.R()) 7219.19366223377

# values with R 1.2.2 and Splus 3.3 ( note Linux is different !)
   good <- if(is.Splus())                             7219.22183565129399 else 
           if(is.R()) {
	     if (Sys.info()[["sysname"]] == "Linux")  7219.22919493643258 else
	     if (Sys.info()[["sysname"]] == "SunOS" ) 7219.22210394543526 else
	                                              7219.22210394543526 #defaulat Solaris
	     } 
   tst  <- sum(hessianVAR)
   error <- max(abs(good - tst))
   cat("max. error ", max(error), "\n")
   
   print.test.value(c(tst), digits=18)

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {#print.test.value(c(tst), digits=18)
      all.ok <- F  
     }

cat("DSE curvature test B 8 ...")

  hessianSS <- hessian(SSmodel)

#  good <- if(is.Splus()) 7844.3395239153897  else if(is.R()) 7841.271986002

# values with R 1.2.2 and Splus 3.3
   good <- if(is.Splus())                             7840.99348875210035 else 
           if(is.R()) {
	     if (Sys.info()[["sysname"]] == "Linux")  7841.35186698713642 else
	     if (Sys.info()[["sysname"]] == "SunOS" ) 7841.24813340843411 else
	                                              7841.24813340843411 #defaulat Solaris
	     } 
   tst  <- sum(hessianSS)
   error <- max(abs(good - tst))
   cat("max. error ", max(error), "\n")
   
   print.test.value(c(tst), digits=18)

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {#print.test.value(c(tst), digits=18)
      all.ok <- F  
     }

cat("DSE curvature test B 9 ...")

  hessianARMA <- hessian(ARMAmodel)
  
#  good <- if(is.Splus()) 1789846677.6677122  else if(is.R()) 90636.84015934517
#  good <- 256440.198630697385

# values with R 1.2.2 and Splus 3.3
   good <- if(is.Splus())                             10711.2666306187384 else 
           if(is.R()) {
	     if (Sys.info()[["sysname"]] == "Linux")  10711.013899114736 else
	     if (Sys.info()[["sysname"]] == "SunOS" ) 10711.2557033145931 else
	                                              10711.2557033145931 #defaulat Solaris
	     } 

   tst  <- sum(hessianARMA)
   error <- max(abs(good - tst))
   cat("max. error ", max(error), "\n")
   
   print.test.value(c(tst), digits=18)

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {#print.test.value(c(tst), digits=18)
      all.ok <- F  
     }


  if (! all.ok) stop("some tests FAILED")
