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


cat("DSE curvature test B 4 ...")

  curvatureVAR <- curvature(VARmodel)$stats

   good <- c(24, 100, 0.05, 0, 0, 0, 0, 1, 1 )
   tst  <- curvatureVAR
   error <- max(abs(good - tst))
   cat("max. error ", max(error))

   if (any(is.na(error)) || any(is.nan(error)) || 1e-5 < error) 
     {print.test.value(c(tst), digits=18)
      all.ok <- F  
     }

cat("DSE curvature test B 5 ...")

  curvatureSS <- curvature(SSmodel, warn=F)$stats

# if(is.R())     tst <- c(48, 100, 0.05, 323.99227471499745, 124.74717834454975,
#        409.3130651286905,  157.59835625487983, 1.0000000021695887, NaN)[-9]
# if(is.Splus()) tst <- c(48, 100, 0.05, 255.387386434704,    98.33696083718515,
#        322.6416248003061, 124.23321787877502, 1, 1.0000000004520695  ) [-9]

# values with R 1.2.2 and Splus 3.3 ( note Linux is different !)
   good <- if(is.Splus())   
               c(48, 100, 0.05, 323.992363637504695, 124.759997096571936,  
	        409.313177468233278, 157.614550723361958, 1.0000000014689332) else 
           if(is.R()) {
	     if (Sys.info()[["sysname"]] == "Linux") 
               c(48, 100, 0.05, 323.989506507851672, 124.76868587630355,  
	        409.309567936195094, 157.625527624176414, 1.00000000058753002) else 
	     if (Sys.info()[["sysname"]] == "SunOS" ) 
               c(48, 100, 0.05, 323.992363637504695, 124.759997096571936,  
	        409.313177468233278, 157.614550723361958, 1.0000000014689332) else 
               c(48, 100, 0.05, 323.992363637504695, 124.759997096571936,  
	        409.313177468233278, 157.614550723361958, 1.0000000014689332) #defaulat Solaris
	     } 

   tst  <- curvatureSS[-9]
   error <- max(abs(good - tst))
   cat("max. error ", max(error))

   print.test.value(c(tst), digits=18)

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {#print.test.value(c(tst), digits=18)
      all.ok <- F  
     }

cat("DSE curvature test B 6 ...")

  curvatureARMA <- curvature(ARMAmodel, warn=F)$stats

# previous test values were suspicious
#  if (print.values) print.test.value(curvatureARMA, digits=digits)
# if(is.R())     tst <- c(71,100,0.05, 8.083857768891898e+23, 354675866.7879350,
#          1.0653907038344889e+24, 467435699.82689363, 1, NaN )[-9]
# if(is.Splus()) tst <- c(72,100,0.05, 31947799733885313024, 1708051.5249938569,
#       42274456907383595008,  2260154.101077503, 1, 1.0000267463352852 )[-9] 
#   tst <- c(71, 100, 0.0500000000000000028, 1.48023299583791679e+23, 
#     7.16541930125915901e+21,  1.95083401806432887e+23, 9.443475294697275e+21,   1 )

# values with R 1.2.2 and Splus 3.3 ( note Linux is different !)
   good <- if(is.Splus())   
               c(48, 100, 0.05, 149.118676845175685, 73.5365097052148116,  
                188.387895177823538, 92.9017650583934795, 1.00005841302234355) else 
           if(is.R()) {
	     if (Sys.info()[["sysname"]] == "Linux") 
               c(48, 100, 0.05, 149.094838773755185, 73.5157177804486111,
	        188.357779539762589, 92.875497745497043, 1.00001977123294594) else 
	     if (Sys.info()[["sysname"]] == "SunOS" ) 
               c(48, 100, 0.05, 149.118676845175685, 73.5365097052148116,  
                188.387895177823538, 92.9017650583934795, 1.00005841302234355) else 
               c(48, 100, 0.05, 149.118676845175685, 73.5365097052148116,  
                188.387895177823538, 92.9017650583934795, 1.00005841302234355) #defaulat Solaris
	     } 

   tst  <- curvatureARMA[-9]
   error <- max(abs(good - tst))
   cat("max. error ", max(error))

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {print.test.value(c(tst), digits=18)
      all.ok <- F  
     }

  if (! all.ok) stop("some tests FAILED")
