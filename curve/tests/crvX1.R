
if(!require("dse2"))  stop("this test requires dse2.")
if(!require("curve"))stop("this test requires curve.")
 Sys.info()
 version.dse()
 
fuzz.small <- 1e-9
fuzz.large <- 1e-6
fuzz.very.large <- 1e-2
digits <- 18
all.ok <- TRUE
test.rng <- list(kind="Wichmann-Hill",seed=c(979,1479,1542),normal.kind="Box-Muller")

  ARMAmodel1<-ARMA(A=array(c(1,.5,.3,  0,.2,.1,  0,.2,.05,  1,.5,.3),c(3,2,2)),
             B=array(c(1,.2,  0,.1,  0,0,  1,.3),c(2,2,2)), C=NULL) 

# sampleT=500

  ARMAmodel1<-l(ARMAmodel1,simulate(ARMAmodel1, rng=test.rng, sampleT=500))
  SSmodel  <- l(to.SS(ARMAmodel1),  ARMAmodel1$data)
  ARMAmodel2<- l(to.ARMA(SSmodel), ARMAmodel1$data)
  
  roots(ARMAmodel1)
  roots(SSmodel)
  roots(ARMAmodel2)
  
  summary(ARMAmodel1)
  summary(SSmodel)
  summary(ARMAmodel2)


  spanARMA1 <- span(ARMAmodel1, warn=FALSE)
  spanSS    <- span(SSmodel, warn=FALSE)
  spanARMA2 <- span(ARMAmodel2, warn=FALSE)

  print(spanARMA1, digits=18)
  print(spanSS, digits=18)
  print(spanARMA2, digits=18)
  if (fuzz.small < abs(sum(spanARMA1) - 235.93017349078193))
      {all.ok <- F; cat("Test FAILED\n"); print(sum(spanARMA1), digits=18)}
  if (fuzz.small < abs(sum(spanSS)    - 239.6366605681105))
      {all.ok <- F; cat("Test FAILED\n"); print(sum(spanSS), digits=18)}
  if (fuzz.small < abs(sum(spanARMA2) - 451.24472649037284))
      {all.ok <- F; cat("Test FAILED\n"); print(sum(spanARMA2), digits=18)}
  
  tfplot(spanARMA1)
  tfplot(spanSS)
  tfplot(spanARMA2)

  curvatureARMA1 <- curvature(ARMAmodel1, warn=FALSE)$stats
  curvatureARMA1def <- curvature(ARMAmodel1, warn=FALSE, compiled=FALSE)$stats
# defaults:      d=0.01, eps=1e-4,r=6, 
  
  if (! test.equal(curvatureARMA1, curvatureARMA1def, fuzz=1e-11))
    {print(curvatureARMA1, digits=18)
     print(curvatureARMA1def, digits=18)
     all.ok <- F
    }
    
  curvatureSS    <- curvature(SSmodel, warn=FALSE)$stats
  curvatureARMA2 <- curvature(ARMAmodel2, warn=FALSE)$stats
  tst  <- rbind(curvatureARMA1, curvatureSS, curvatureARMA2)
  print(tst, digits=18)

#Linux values work on Solaris with fuzz.large=1e-6
#Parms Sample Sign. level      RMS Parameter      RMS Intrinsic  c*sqrt(F) Parameter c*sqrt(F) Intrinsic     Min Axis Ratio Max Axis Ratio
  good <- t(matrix(c(  
11, 1000, 0.05, 0.3562155343985393, 0.3009020761040172, 0.47768990044815857, 0.4035137968406400, 1.0000000002443807, 1.32774403734037660,
16, 1000, 0.05, 0.6016763965085490, 0.4849528063707048, 0.77375522775520988, 0.6236488107582950, 1.0000000002181073,		     NaN,
18, 1000, 0.05, 0.9665129322089927, 0.9901210419889496, 1.22802952756715111, 1.2580254592652123, 1.0004326355017099,		     NaN),
9,3))

   error <- max(abs(good[, -9] - tst[, -9]))
   cat("max. error ", error)
     
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.large < error) 
     {printTestValue(c(t(tst)), digits=18); all.ok <- F ; cat("Test FAILED") }

# sampleT=100

  ARMAmodel1<-l(ARMAmodel1,simulate(ARMAmodel1, rng=test.rng, sampleT=100))
  SSmodel  <- l(to.SS(ARMAmodel1),  ARMAmodel1$data)
  ARMAmodel2<- l(to.ARMA(SSmodel), ARMAmodel1$data)
  curvatureARMA1 <- curvature(ARMAmodel1, warn=FALSE)$stats
  curvatureARMA1def <- curvature(ARMAmodel1, warn=FALSE, compiled=FALSE)$stats
  if (! test.equal(curvatureARMA1, curvatureARMA1def, fuzz=1e-11))
    {print(curvatureARMA1, digits=18)
     print(curvatureARMA1def, digits=18)
     all.ok <- F
    }
  curvatureSS    <- curvature(SSmodel, warn=FALSE)$stats
  curvatureARMA2 <- curvature(ARMAmodel2, warn=FALSE)$stats
  tst  <- rbind(curvatureARMA1, curvatureSS, curvatureARMA2)
  print(tst, digits=18)

#These values were on Linux Mandrake 9.0 (with gcc 2.9x? or at least < 3.2.2)
#They work on Solaris with fuzz 1e-6
#On Linux Mandrake 9.1 (with gcc 3.2.2) requires fuzz set to 1e-5
#Parms Sample Sign. level      RMS Parameter      RMS Intrinsic  c*sqrt(F) Parameter c*sqrt(F) Intrinsic     Min Axis Ratio Max Axis Ratio
  good <- t(matrix(c(
11, 200, 0.05, 0.8072684390607161, 0.68050736181498195, 1.09491154342599750, 0.9229833965818595, 1.00000000089819618, NaN,
16, 200, 0.05, 1.4200015989565475, 1.11328625807171511, 1.85076843299700444, 1.4510089741043570, 1.00000000417908952, NaN,
18, 200, 0.05, 2.3319228789506332, 2.28217641539719374, 3.00519272189918629, 2.9410835219078497, 1.00034641227181775, NaN),
9,3))

#Linux Mandrake 9.1 (with gcc 3.2.2) values are
#c(t(11,                   16,                   18,                  
#200,                  200,                  200,  
#0.0500000000000000028,  0.0500000000000000028,  0.0500000000000000028,  
#0.807268437301731923,  1.42000160081664561,  2.33192258659243423,  
#0.68050734798853596,  1.11328625926179581,  2.28217907095856853,  
#1.09491154104025812,  1.85076843542137537,  3.00519234513169531,  
#0.922983377828822005,  1.45100897565545672,  2.94108694417956729,  
#1.00000000063105121,   1.0000000020690627,  1.00034613602572842,                  
#NaN,                  NaN,                  NaN))

   error <- max(abs(good[, -9] - tst[, -9]))
   cat("max. error ", error)
     
   if (any(is.na(error)) || any(is.nan(error)) || 1e-5 < error) 
     {printTestValue(c(t(tst)), digits=18); all.ok <- F ; cat("Test FAILED") }


# sampleT=50

  ARMAmodel1<-l(ARMAmodel1,simulate(ARMAmodel1, rng=test.rng, sampleT=50))
  SSmodel  <- l(to.SS(ARMAmodel1),  ARMAmodel1$data)
  ARMAmodel2<- l(to.ARMA(SSmodel), ARMAmodel1$data)
  curvatureARMA1 <- curvature(ARMAmodel1, warn=FALSE)$stats
  curvatureARMA1def <- curvature(ARMAmodel1, warn=FALSE, compiled=FALSE)$stats
  if (! test.equal(curvatureARMA1, curvatureARMA1def, fuzz=1e-11))
    {print(curvatureARMA1, digits=18)
     print(curvatureARMA1def, digits=18)
     all.ok <- F
    }
  curvatureSS    <- curvature(SSmodel, warn=FALSE)$stats
  curvatureARMA2 <- curvature(ARMAmodel2, warn=FALSE)$stats
  tst  <- rbind(curvatureARMA1, curvatureSS, curvatureARMA2)
 print(tst, digits=18)

#Linux values work on Solaris with 1e-5
# Warning message: 
#acceration space dimension reduced for smaller sample space. in: curvature.Darra#Parms Sample Sign. level      RMS Parameter      RMS Intrinsic  c*sqrt(F) Parameter c*sqrt(F) Intrinsic     Min Axis Ratio Max Axis Ratio
  good <- t(matrix(c(
 11, 100, 0.05, 1.1902331226374676, 1.02413128081527893, 1.6397187390037589, 1.4108893631119579, 1.00000000108466902, NaN,
 16, 100, 0.05, 2.1698167661684429, 1.64496591418322957, 2.8829778015213545, 2.1856224400108535, 1.00000000056460969, NaN,
 18, 100, 0.05, 4.2278364145158402, 3.70068582243928734, 5.5624870633970582, 4.8689246684992131, 1.01584061323199637, NaN),
9,3))

   error <- max(abs(good[, -9] - tst[, -9]))
   cat("max. error ", max(error))
     
   if (any(is.na(error)) || any(is.nan(error)) || 1e-5 < error) 
     {printTestValue(c(tst), digits=18); all.ok <- F }

# larger model (sampleT=500)
  
  ARMAmodel1<-ARMA(A=array(c(1,.7,  0,.2,  0,.2,
                             0,.2,  1,.5,  0,.2,
			     0,.2,  0,.2,  1,.3),c(2,3,3)),
                   B=array(c(1,.2,  0,.1,  0,.1,  
	                     1,.2,  1,.1,  0,0,
			     0,0,   0,.1,  1,.3),c(2,3,3)), C=NULL) 

  ARMAmodel1<-l(ARMAmodel1,simulate(ARMAmodel1, rng=test.rng, sampleT=500))
  SSmodel  <- l(to.SS(ARMAmodel1),  ARMAmodel1$data)
  ARMAmodel2<- l(to.ARMA(SSmodel), ARMAmodel1$data)
  
  roots(ARMAmodel1)
  roots(SSmodel)
  roots(ARMAmodel2)
  
  summary(ARMAmodel1)
  summary(SSmodel)
  summary(ARMAmodel2)

  curvatureARMA1 <- curvature(ARMAmodel1, warn=FALSE)$stats
  curvatureSS    <- curvature(SSmodel, warn=FALSE)$stats
  curvatureARMA2 <- curvature(ARMAmodel2, warn=FALSE)$stats
  tst  <- rbind(curvatureARMA1, curvatureSS, curvatureARMA2)
 print(tst, digits=18)

#closer to unit root may be more interesting


data("egJofF.1dec93.data", package="dse1")

if(!exists("egJofF.1dec93.data")) stop("egJofF.1dec93.data does not exist")


  cat("Curvature calculations with bft estimated model \n")
   eg.data<- egJofF.1dec93.data
   output.data(eg.data) <- output.data(eg.data, series=c(1,2,6)) #CPI GDP employment
 # following is optional 
 # tframe(output.data(eg.data))<- tframe(output.data(egJofF.1dec93.data))

  bft.model <- bft(trim.na(eg.data), max.lag=3, verbose=FALSE) 

tfplot(bft.model)
summary(bft.model)
stability(bft.model)
span.bft <- span(bft.model)
span.bft
tfplot(span.bft)

z <- l(fix.constants(bft.model, fuzz=0.16), trim.na(eg.data))
summary(z)
span(z)
tfplot(z)

# comparison values come only from a previous run of the code.

  tst.span.bft <- span(bft.model)
  # good from R 1.6.2
#These values were on Linux Mandrake 9.0 (with gcc 2.9x? or at least < 3.2.2)
#They work on Solaris with fuzz 1e-7
#On Linux Mandrake 9.1 (with gcc 3.2.2) and R 1.7.0beta requires fuzz set to 1e-6
  good <- c(24.1219109405661953,  21.4893351371261616,  16.7891067422431135,
            13.4030221805392742,  11.5605466485710267,   9.4835864591854051,  
	     8.82882082257517808,  8.25373217057830288,   8.0290715128405985,  
	     7.57911845828720399,  7.23605710138632752,  6.90447759423024721,  
	     6.55052218588380164,  6.23909660487700801,  5.72720169335468654,  
	     5.45900757259505731,  5.12703117229429051,  4.92008418965641514,  
	     4.65487584871777216,  4.48083390079675414,  4.16509269292236795,  
	     4.05004413755998183,  3.90038798838493905,  3.61137294458389269,  
	     3.55746336585901934,  3.37159592090552485,  3.11904794121582984,  
	     2.99404224126606833,  2.65553768047008809,   2.4108538633271408,  
	     2.22737422060760393,  2.06171737746169814,  1.83057703088174772,  
	     1.57728269747513683,  1.40294643857774792,  4.34329853397508004e-07,  
	     1.73654008174314229e-10,  8.25307512598579223e-11,  7.31355766042255533e-11, 
	     6.84027327043274035e-11,  6.38814748388448426e-11,  5.52551798419003698e-11, 
	     5.238445497649314e-11,  5.09239749267603534e-11,  4.80811475415604721e-11,  
	     4.64794174624757404e-11,  4.33817394951145184e-11,  4.256295675280607e-11,  
	     3.7944872921781713e-11,  3.62910362674584484e-11,  3.57030600848564127e-11,  
	     3.2338477623622031e-11,  2.99157310546484913e-11,  2.67517904415068739e-11,  
	     2.6222061631566316e-11,  2.36791357839587041e-11,  2.21241760069945202e-11,  
	     2.10934060495354299e-11,  1.75140996365650319e-11,  1.52976389733687478e-11)

   error <- max(abs(good - tst.span.bft))
   cat("max. error ", error, "\n")
   
   if (any(is.na(error)) || any(is.nan(error)) || 1e-6 < error) 
     {if (any(is.na(error)))  cat("contains na's: ",  is.na(error),  "\n")
      if (any(is.nan(error))) cat("contains nan's: ", is.nan(error), "\n")
      printTestValue(c(tst.span.bft), digits=18)
      all.ok <- FALSE  
     }

  tst.bft  <- curvature(bft.model)$stats
# R 1.6.2
# Linux
  good <- c(60, 708,  0.0500000000000000028,  77360922883252.1562,  
  90560814571134.1719,  89527271605117.6875,  104803075515608.328,  
  1.00000000000055755, NaN)

# Solaris
#  good <- c(60, 708,  0.0500000000000000028,  83649469902256.2344,
#  86243510860756,  96804801861854.2812,  99806800814179.1562,  
#  1.00000000000008527, NaN)

# Warning messages: 
#1: eliminating degenerate subspace for R11. in: rel.curvature(s.sqr, R[1:p, 1:p], R[, (p + 1):m], show.extra.details = show.extra.details) 
#2: acceration space dimension reduced for smaller sample space. in: curvature.Darray(genD(func, x = x, func.args = func.args, d = d,  
#3: B is not symmetric. in: effective.curvature(cur, QRofD, residual, s.sqr, show.details = show.extra.details,  

   error <- max(abs(good - tst.bft)[-9])
   cat("max. error ", error, "\n")

# there is not much point in testing this yet   
#   if (any(is.na(error)) || any(is.nan(error)) || 1e-9 < error) 
#     {if (any(is.na(error)))  cat("contains na's: ",  is.na(error),  "\n")
#      if (any(is.nan(error))) cat("contains nan's: ", is.nan(error), "\n")
      printTestValue(c(tst.bft), digits=18)
#      all.ok <- FALSE  
#     }


   arma.model <- to.ARMA(bft.model)

   tst.span.arma <- span(arma.model)
   print(sum(tst.span.arma), digits=18)
   good <- 2219.292822557674

   error <- max(abs(good - sum(tst.span.arma)))
   cat("max. error ", max(error))

   if (any(is.na(error)) || any(is.nan(error)) || 1e-7 < error) 
     {if (any(is.na(error)))  cat("contains na's: ",  is.na(error),  "\n")
      if (any(is.nan(error))) cat("contains nan's: ", is.nan(error), "\n")
      printTestValue(c(sum(tst.span.arma)), digits=18)
      all.ok <- FALSE  
     }

   tst.arma  <- curvature(arma.model)$stats
#Warning messages: 
#1: acceration space dimension reduced for smaller sample space. in: curvature.Darray(genD(func, x = x, func.args = func.args, d = d,  
#2: N.B. (I-B) is not positive definite as it should be at a local min!  Axis ratio will evaluate to NA or NaN.  Eigenvalues of B  35.6610439830782 in: effective.curvature(cur, QRofD, residual, s.sqr, show.details = show.extra.details,  

   print(tst.arma, digits=18)
#         Parms Sample   Sign. level          RMS Parameter
#            RMS Intrinsic    c*sqrt(F) Parameter     c*sqrt(F) Intrinsic 
#	      Min Axis Ratio   Max Axis Ratio

# R 1.6.2 Linux Mandrake 9.0
   good <- c(75, 708, 0.0500000000000000028,  81.0320265417880279,
         32.0032463681021468, 92.6067497650618492, 36.574632951479721,
	  1.00012969286955733, NaN)
#  R 1.7.0beta Linux Mandrake 9.1
#   good <- c(75, 708, 0.0500000000000000028, 81.0346252569432721,
#        32.0087997790132945,  92.6097196841711536,  36.5809796190448751,  
#          1.00012337610921764, NaN)
 
# R 1.6.2 Solaris
#   good <- c(75, 708,  0.0500000000000000028,  81.0269058238850874,
#         31.9989906362728398,  92.6008975969538568,  36.5697693252148639,
#	 1.00013468994302435,  NaN)

   error <- max(abs(good - tst.arma)[-9])

   if (any(is.na(error)) || any(is.nan(error)) || 1e-1 < error) 
     {if (any(is.na(error)))  cat("contains na's: ",  is.na(error),  "\n")
      if (any(is.nan(error))) cat("contains nan's: ", is.nan(error), "\n")
      printTestValue(c(tst.arma), digits=18)
      all.ok <- FALSE  
     }


  if (! all.ok) stop("some tests FAILED")
