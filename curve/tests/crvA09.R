# Tests of DSE curvature functions from dsecurvature.function.testsA
if(!require("dse2"))  stop("this test requires dse2.")
if(!require("curve"))stop("this test requires curve.")
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

  ARMAmodel1<-ARMA(A=array(c(1,.5,.3,0,.2,.1,0,.2,.05,1,.5,.3),c(3,2,2)),
             B=array(c(1,.2,0,.1,0,0,1,.3),c(2,2,2)), C=NULL) 


  ARMAmodel1<-l(ARMAmodel1,simulate(ARMAmodel1, rng=test.rng))
  SSmodel  <- l(to.SS(ARMAmodel1),  ARMAmodel1$data)
  ARMAmodel<- l(to.ARMA(SSmodel), ARMAmodel1$data)


#  if(is.Splus()) good <-  c(
#         12.98465581800713764, 12.04846041797238776, 11.75207879184403659,
#         10.86851614803762800,  9.41015226452319986,  9.14102587417263912,
#         8.81980414072582874,  8.32921589456761780,   5.12046284222839354,
#         5.02836144306142163,  3.44776911209970471,   3.14630978631044567,
#         2.20492937840385128,  2.07975512652550254,   1.93518614772789666,
#         1.87713878204627571 )
#  if(is.R())     good <-  c(
#         14.3191707820078413, 12.7553932264152472, 12.0399559044862077,
#         11.1780750248459100,  8.7642374743163280,  8.4110650474333628,
#          8.1119238141018251,  7.7112054454257013,  4.5816849188761744,
#          4.4928124488652665,  3.0909505682400633,  2.7790878709102365,
#          2.2901166620448916,  2.1211736442161429,  1.8231113211019712,
#          1.7364603432999757 )         #f77 on Sun5 R0.61.1

#      c( 1.431917078200739e+01, 1.275539322641331e+01, 1.203995590448420e+01,
#       1.117807502484729e+01, 8.764237474316120, 8.411065047433565,
#       8.111923814106442,     7.711205445427577, 4.581684918879841,
#       4.492812448867578,     3.090950568240565, 2.779087870909759,
#       2.290116662046349,     2.121173644211920, 1.823111321101024,
#       1.736460343300366)))  # 0.49

#  good <-  c(9.470240142187398,  8.880447829525087,  7.552501916476404,
#            7.531576128069482,  7.169071366088333,  6.868663706072056,  
#            5.629811912135156,  5.343700909140345,  3.984898912643869,  
#            3.899294396457272,  2.444843033246666,  2.207333797038643,  
#            1.447912667824839,  1.363598776819928,  1.256253315103265,  
#            1.205277678709112)
  good <- c(12.6684354906294967,  11.4927966118834544,  11.0537104821112262,
           10.3371110331134979,  9.77605886766110288,  9.21103083988342064, 
	    8.16488126052858476,   7.7872986620231881,  4.79640068843280787, 
	    4.69148925882611323,  3.23108014824838419,  2.94195878161703162,  
	    2.11511171224067507,  2.01921578796399359,  1.85314036275455907,  
	    1.79333516556799233)

#gctorture()  # to help localize segfault

cat("DSE curvature test A 9a...\n")

   spanSS <- span(SSmodel, compiled= F)  # should be .DSECOMPILED)
   error <- max(abs(good - spanSS))
   cat("max. error ", max(error))
    
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.large < error) 
     {printTestValue(c(spanSS), digits=18); all.ok <- F }

cat("DSE curvature test A 9b...\n")

   spanSS <- span(SSmodel, compiled=FALSE) 
   error <- max(abs(good - spanSS))
   cat("max. error ", max(error))
    
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.large < error) 
     {printTestValue(c(spanSS), digits=18); all.ok <- F }

  if (! all.ok) stop("some tests FAILED")
