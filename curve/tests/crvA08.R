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
  

cat("DSE curvature test A 8a...")
# from user guide

  VARmodel<-ARMA(A=array(c(1,.5,.3,0,.2,.1,0,.2,.05,1,.5,.3),c(3,2,2)),
             B=array(c(1,.2,0,.1,0,0,1,.3),c(2,2,2)), C=NULL) 

# Note this gives a terrible fit.
  VARmodel<-l(VARmodel,simulate(VARmodel, rng=test.rng))
  SSmodel  <- l(to.SS(VARmodel),  VARmodel$data)
  ARMAmodel<- l(to.ARMA(SSmodel), VARmodel$data)

  spanVAR.f <- span(VARmodel, compiled=.DSECOMPILED) 

#  if(is.Splus()) good <-  c(
#      18.28342261104111799, 16.86252848720547703, 15.47359559657559025,
#      12.60857465926911836, 10.85404071073360299, 10.80245325753396202,
#       9.79091689469922066,  9.21160820465702379,  2.32881291213965458,
#       1.97834187456894184,  1.50902538655984464 )
#  if(is.R())     good <-  c(
#        17.8366876297672725, 16.7671374021287249, 16.2282896291638892,
#        13.3626055311535605, 10.3267790597707432,  9.9307009769349008,
#         8.8813776725358053,  8.2133988282015906,  2.2202966179396091,
#         1.9141400784167535,  1.4838767428262707 )     #f77 on Sun5 R0.61.1
#     c( 1.783668762976735e+01, 1.676713740212845e+01, 1.622828962916468e+01,
#     1.336260553115620e+01, 1.032677905977513e+01, 9.930700976937008,
#     8.881377672535175, 8.213398828201516, 2.220296617939313,
#     1.914140078418723, 1.483876742824580)))  #R 0.49
#  good <-c(13.30038806708721,  12.23076792331929,  9.752849398484228, 
#         8.786241709293281,  8.583699014891037,  7.624249106943636,  
#         5.998914711634687,   5.61060655287252,  1.774590614465879,  
#         1.436530005339643,  0.9712521284803599)
   good <- c( 17.3672899077923439,  15.9424331056250885,  14.5699882827285894,
             11.6462671196918794,  10.7723733430260129,  10.5668728474602247, 
	      8.76191214115850769,  7.93300083097722553,  2.24655970809629668,  
	      1.88819440545349604,  1.47455500924911109)
 
   tst  <- spanVAR.f 
   error <- max(abs(good-tst))
   cat("max. error ", max(error))
   
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.large < error) 
     {print.test.value(c(tst), digits=18); all.ok <- F }

cat("DSE curvature test A 8b...")
  spanVAR <- span(VARmodel, compiled=FALSE) 

   tst  <-  spanVAR
   good <-  spanVAR.f
   error <- max(abs(good-tst))
   cat("max. error ", max(error))
   
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {print.test.value(c(tst), digits=18); all.ok <- F }

  if (! all.ok) stop("some tests FAILED")
