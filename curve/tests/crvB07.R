# Tests of DSE curvature functions (previously dsecurvature.function.testsB )

# comparison values come only from a previous run of the 
#  code (theoretical values would be nice)...
# Test values have been changed with change to toARMA in 2001.2 which
# eliminates near zero parameter values using fixConstants. The result is
# much more stable and believable curvature results. The span results do not
# change much (as would be hoped) but do change more than the tolerance of 
# these tests. Old values in comments are  strictly for historical reference.

if(!require("dse1"))  stop("this test requires dse1.")
if(!require("curve"))stop("this test requires curve.")
 Sys.info()
 DSEversion()

# Linux versions vary a lot on these tests. see also test specific settings below
digits <- 18
all.ok <- TRUE 

 data("eg1.DSE.data.diff", package="dse1") 

# data size affects memory constraints
  data1 <- eg1.DSE.data.diff
   inputData(data1) <- NULL
  outputData(data1) <- outputData(data1)[1:50,1:2]

  VARmodel <- estVARXls(data1, re.add.means=FALSE)
  SSmodel  <- l(toSS(VARmodel),  data1)
  ARMAmodel<- l(toARMA(SSmodel), data1)


cat("DSE curvature test B 7 ...")

  print(hessianVAR <- hessian(VARmodel))
  

## R 1.2.3    if (Sys.info()[["sysname"]] == "Linux")  7219.22919493643258 else
## R 1.3.0    if (Sys.info()[["sysname"]] == "Linux")  7219.17128513884654 else
## R 1.4.1  and R 1.5.1 in Mandrake 8.1
## svd until R 1.6.1 in
##   Mandrake 9.0 gave 7219.17325840329431
##   Solaris	  gave 7219.22210394543526
## switching to La.svd 
##    atlas Linux goves 7219.2019407758653)

##  good <- if(!is.R()) 7219.717083137912   else  7219.19366223377

## values with R 1.2.2 and Splus 3.3 ( note Linux is different !)
#   good <- if(!is.R())                                7219.22183565129399 else 
#	     if (Sys.info()[["sysname"]] == "Linux")  7219.18588980847835 else
#	     if (Sys.info()[["sysname"]] == "Windows")7219.2222921110606  else
#	     if (Sys.info()[["sysname"]] == "SunOS" ) 7219.19791566255844 else
#	                                              7219.19791566255844 #default Solaris
#

# DIFFERENCE FROM ABOVE FOUND MARCH 2006
#  WHEN SEPARATING OUT numDeriv. THERE WAS AN ERROR.
# HESSIAN USED FIRST COLUMNS IN genD RATHER THAN SKIPPING THEM. 
   good <-  6293.68641833058791
	      
   printTestValue(tst  <- sum(hessianVAR), digits=18)
   error <- max(abs(good - tst))
   cat("max. error ", error, "\n")
   
# Linux seems to have large variance here
   if (any(is.na(error)) || any(is.nan(error)) ||  0.05 < error)
	{cat("test FAILED."); all.ok <- FALSE } 



cat("DSE curvature test B 8 ...")

  print(hessianSS <- hessian(SSmodel))

#  good <- if(!is.R()) 7844.3395239153897  else       7841.271986002


## R 1.2.3    if (Sys.info()[["sysname"]] == "Linux")  7841.35186698713642 else
## R 1.3.0    if (Sys.info()[["sysname"]] == "Linux")  7841.24925650512705 else
## R 1.4.1 in Mandrake 8.1
## svd until R 1.6.1 in
##	 Mandrake 9.0 gave 7841.22786316239581
##	 Solaris      gave 7841.24813340843411
## switching to La.svd  
##   Linux atlas gives 7841.58938345190381
# values with R 1.2.2 and Splus 3.3
#   good <- if(!is.R()) 			       7840.99348875210035 else 
#	      if (Sys.info()[["sysname"]] == "Linux")  7841.31080366921742 else
#	      if (Sys.info()[["sysname"]] == "Windows")7841.07715279903    else
#	      if (Sys.info()[["sysname"]] == "SunOS" ) 7841.20397396244061 else
#						       7841.20397396244061 #defaulat Solaris
#      
# DIFFERENCE FROM ABOVE FOUND MARCH 2006 (see above)
  good <- 6292.77943518256779 
  #This is now close to above test result !!!! check other matrix norms!!
   printTestValue(c(tst  <- sum(hessianSS)), digits=18)
   error <- max(abs(good - tst))
   cat("max. error ",error, "\n")
   
# Linux seems to have large variance here
   if (any(is.na(error)) || any(is.nan(error)) || 0.5 < error) 
	{cat("test FAILED."); all.ok <- FALSE } 


cat("DSE curvature test B 9 ...")

  print(hessianARMA <- hessian(ARMAmodel))
  
##  good <- if(is.Splus()) 1789846677.6677122  else if(is.R()) 90636.84015934517
##  good <- 256440.198630697385

##		 # 10711.013899114736	 R 1.2.3
##		 # 10711.8649535454842 # R 1.3.0
## values with R 1.2.2 and Splus 3.3
## R 1.4.1 and R 1.5.1 in Mandrake 8.1 and
## svd until R 1.6.1 in
##     Mandrake 9.0 gave 10710.1343804970293
##     Solaris      gave 10711.2557033145931
## switching to La.svd  
##     atlas Linux gives 10710.7368676075766
#   good <- if(!is.R()) 			       10711.2666306187384 else 
#	      if (Sys.info()[["sysname"]] == "Linux")  10711.8642678423475 else
#	      if (Sys.info()[["sysname"]] == "Windows")10711.271347343158  else
#	      if (Sys.info()[["sysname"]] == "SunOS" ) 10710.6870697151153 else
#						       10710.6870697151153 #defaulat Solaris
#	       
#	      #if (Sys.info()[["sysname"]] == "Linux")  10712.2643162879049 else
# DIFFERENCE FROM ABOVE FOUND MARCH 2006 (see above)
   good <- 4303.25552473151038
   printTestValue(tst  <- sum(hessianARMA), digits=18)
   error <- max(abs(good - tst))
   cat("max. error ", error, "\n")

# Linux seems to have large variance here
    if (any(is.na(error)) || any(is.nan(error)) || 0.5 < error)
	{cat("test FAILED."); all.ok <- FALSE } 


  if (! all.ok) stop("some tests FAILED")
