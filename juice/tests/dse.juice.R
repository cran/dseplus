  Sys.info()

  require("mva"); require("ts")
  require("dse2") # adds dse1, tframe, and syskern
  version.dse()

  require("juice") 

  if (is.R()) data("egJofF.1dec93.data", package="dse1")
  if (is.S()) source(paste(DSE.HOME, "/data/egJofF.1dec93.data.R", sep=""))

  if (is.S()) test.rng <- list(kind="default", normal.kind="default", 
                       seed=c(13,44,1,25,56,0,6,33,22,13,13,0))
  if (is.R()) test.rng <- list(kind="default", normal.kind="default",
                     seed=c( 979)) #, 1479, 1542))

 
fuzz.small <- 1e-12
digits <- 18
all.ok <- TRUE  

ets <- FALSE

  set.RNG(test.rng)
  data00 <- matrix(rnorm(300), 100,3)
  data0 <- TSdata(output=data00)
  data1 <- TSdata(output=data0$output %*% matrix(rnorm(9),3,3))


cat("dse juice test 0 ... ")
  z <- concentrate(data00, conc=est.projection(data00, n=3))
  # next is true because all PCs are used
  ok <-      test.equal(data00, reconstitute(z), fuzz=fuzz.small)
# T but does this have meaning  ok <- ok & test.equal(data00,reconstitute(z), fuzz=fuzz.small)
  ok <- ok & all(seriesNames(reconstitute(z))
                == paste("recon.", seriesNames(data00)))
  ok <- ok & all(concentrated.seriesNames(z)
                == paste("concentrate", seq(length=concentrated.dimension(z))))
  ok <- ok & all(seriesNames(z) == seriesNames(data00))		   		   

  z <- concentrate(output.data(egJofF.1dec93.data)) # p=1
  z <-tfwindow(z, start=c(1981,1))
  ok <- all(c(1981,1) == start(z)) & all(c(1981,1) == start(reconstitute(z)))
  all.ok <- all.ok & ok
  if (ok) cat("ok\n") else  cat("failed!\n") 


  cat("dse juice test 1 ... ")
  z <- concentrate(data0, conc=est.projection(data0, p=3))
  # next is true because all PCs are used
  ok <-      test.equal(data0, reconstitute(z), fuzz=fuzz.small)
# T but does this have meaning  ok <- ok & test.equal(data0,TSdata(z), fuzz=fuzz.small)
  ok <- ok & all(seriesNamesOutput(reconstitute(z))
                == paste("recon.", seriesNamesOutput(data0)))
  ok <- ok & all(concentrated.seriesNamesOutput(z)
                == paste("concentrate", seq(length=concentrated.nseriesOutput(z))))
  ok <- ok & all(seriesNamesOutput(z) == seriesNamesOutput(data0))		   		   
  all.ok <- all.ok & ok 
  if (ok) cat("ok\n") else  cat("failed!\n") 


  cat("dse juice test 2 ... ")
  z <- concentrate(data0, conc=est.projection(data0, p=1))
  ok <-       all(seriesNamesOutput(data0) == seriesNamesOutput(z))
  ok <- ok & "concentrate 1" == concentrated.seriesNamesOutput(z)
  ok <- ok &  ! test.equal(data0, reconstitute(z), fuzz=fuzz.small)
  all.ok <- all.ok & ok 
  if (ok) cat("ok\n") else  cat("failed!\n") 


  cat("dse juice test 3 ... ")
  z <- concentrate(data0, conc=est.projection(data0, p=2))
  z <- concentrate(data0, conc=est.projection(data0)) # p=1

  ok <- test.equal(z,   concentrate(data0,conc=z), fuzz=fuzz.small)
  ok <- ok & (3 == nseriesOutput(z)) & (3 == nseriesOutput(TSdata(z))) &
             (1 == concentrated.nseriesOutput(z))               &
             (100 == periods(z)) & (100 == periods(TSdata(z)))
  all.ok <- all.ok & ok
  if (ok) cat("ok\n") else  cat("failed!\n") 


  cat("dse juice test 4 ... ")
  z <- concentrate(egJofF.1dec93.data) # p=1
  ok <- all(c(1974,2) == start(z)) & all(c(1974,2) == start(reconstitute(z)))
  all.ok <- all.ok & ok
   {if (ok) cat("ok\n") else  cat("failed!\n") }


  cat("dse juice test 5 ... ")
  z <-tfwindow(z, start=c(1981,1))
  ok <- all(c(1981,1) == start(z)) & all(c(1981,1) == start(reconstitute(z)))
  all.ok <- all.ok & ok
  if (ok) cat("ok\n") else  cat("failed!\n") 


  cat("dse juice test 6 ... ")
#  zm <- est.concentrated.model(z, scale=TRUE, center=TRUE, estimation="bft",
#  z is already a concentrated object so center and scale are not used
  zm <- est.concentrated.model(z, estimation="bft",
                                estimation.args=list(max.lag=2, verbose=FALSE))
  #tfplot(zm)
  zmr <- check.residuals(zm, ac=FALSE, pac=FALSE, plot.=FALSE)
# S values
tst  <- zmr$skewness 
good <- c(-1.979814544376786,  0.4241402252105713, -0.2758381743875303,
         -0.5754191393269468, 1.983488003995619,   0.1565333254836311,
         -0.116925097663112,  0.1326786266039132, -0.2649736832254932,
         -0.3710824450042611) 

   error <- max(abs(good - tst))
   cat("max. error ", max(error), "\n")
   
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {if (any(is.na(error)))  cat("na's: ",  is.na(error), "\n")
      if (any(is.nan(error))) cat("nan's: ", is.nan(error), "\n")
      if (fuzz.small < error) cat("error: ", error, "\n")
      print.test.value(c(tst), digits=18)
      all.ok <- FALSE 
     }

  
tst  <- zmr$kurtosis
good <- c(11.01131667332406,   3.503769190958172, 5.024697760083122,
          4.972991895880567, 15.5906633209087,   2.3188280244819,
          2.935474704110169,  4.131056248843817, 3.823872187369519,
          4.900819828268634)

   error <- max(abs(good - tst))
   cat("max. error ", max(error), "\n")
   
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {if (any(is.na(error)))  cat("na's: ",  is.na(error), "\n")
      if (any(is.nan(error))) cat("nan's: ", is.nan(error), "\n")
      if (fuzz.small < error) cat("error: ", error, "\n")
      print.test.value(c(tst), digits=18)
      all.ok <- FALSE
     }

  zmf <- featherForecasts(zm)
#  ok <- ???

  if (ets)
    {cat("  juice test using ets ...\n")
     JofF.VAR.data.names <- TSPADIdata(
	input = "B14017", input.transform="diff", input.names="R90",
	output =   c("B820600", "I37026", "b1627", "b14013", "b4237", "D767608",
                      "b3400", "M.JQIND", "M.CUSA0"), # P484549=CPI discontinued
	output.names = c( "CPI", "GDP", "M1", "RL", "TSE300", "employment", 
                       "PFX", "US ind. prod.", "US CPI"),

	output.transforms = c("percent.change", "percent.change", 
                  "percent.change", "diff", "diff", "percent.change", 
                     "percent.change", "percent.change",  "percent.change"),
        server="ets", start.server=TRUE, server.process="fame.server", 
        cleanup.script="cleanup.fame.server", stop.on.error=TRUE, warn=TRUE)


#         c("ets","", "M.BCPI",  "percent.change", "com. price ind."),

#   availability(JofF.VAR.data.names)
      JofF.VAR.data <- freeze(JofF.VAR.data.names)

      z <- concentrate(JofF.VAR.data, p=2)
      z <- concentrate(JofF.VAR.data, p=3)
      z <- concentrate(JofF.VAR.data, p=2, scale=FALSE)
      ok <- test.equal(seriesNames(reconstitute(z)), 
                       seriesNames(JofF.VAR.data)) 

}



cat("juice graphics tests ...\n")
    # If no device is active then write to postscript file 
    if ( dev.cur() == 1 ) {
                postscript(file="zot.postscript.test.ps", width=6, 
                        height=6, pointsize=10, onefile=FALSE, 
                        print.it=FALSE, append=FALSE)
                on.exit((function() {
                        dev.off()
                        synchronize(1)
                        rm("zot.postscript.test.ps")
                }))
        }


   set.RNG(test.rng)
    data0 <- TSdata(output=matrix(rnorm(300), 100,3))
    seriesNames(data0)<- 
                  list(output=paste("data0", seq(length=nseriesOutput(data0))))
    data1 <- TSdata(output=data0$output %*% matrix(rnorm(9),3,3))
    seriesNames(data1)<- 
                  list(output=paste("data1", seq(length=nseriesOutput(data1))))

cat("  juice graphics test 1 ...\n")
    z <- concentrate(data0) # p=1
    ok <-  all( seriesNamesOutput(z) 
                    == paste("data0", seq(nseriesOutput(z))))

    concentrated.tfplot(z)
    tfplot(z)
    all.ok <- all.ok  & ok
    if (ok) cat("ok\n") else  cat("failed!\n") 

cat("  juice graphics test 2 ...\n")
    tfplot(reconstitute(z))
    ok <-  all( seriesNamesOutput(reconstitute(z))
                      == paste("recon. data0", seq(nseriesOutput(data0))))
    all.ok <- all.ok  & ok
    if (ok) cat("ok\n") else  cat("failed!\n") 

cat("  juice graphics test 3 ...\n")
    z <- concentrate(data0, p=2)
    tfplot(z)
    all.ok <- all.ok  & ok
    if (ok) cat("ok\n") else  cat("failed!\n") 

cat("  juice graphics test 4 ...\n")
    tfplot(reconstitute(z))
    all.ok <- all.ok  & ok
    if (ok) cat("ok\n") else  cat("failed!\n") 

cat("  juice graphics test 5 ...\n")
    z <- concentrate(data0, p=3)
    tfplot(z)
    all.ok <- all.ok  & ok
    if (ok) cat("ok\n") else  cat("failed!\n") 
    cat("   (note that actual and reconstructed coincide.)\n")

cat("  juice graphics test 6 ...\n")
    tfplot(reconstitute(z))
    all.ok <- all.ok  & ok
    if (ok) cat("ok\n") else  cat("failed!\n") 

cat("  juice graphics test 7 ...\n")
    z <- concentrate(egJofF.1dec93.data, p=2)
    tfplot(z)
    z <- tfwindow(z, start=c(1992,1))
    tfplot(z)
    tfplot(reconstitute(z))
    all.ok <- all.ok  & ok
    if (ok) cat("ok\n") else  cat("failed!\n") 

  if (ets)
    {cat("  juice graphics test 7 using ets ...\n")
     JofF.VAR.data.names <- TSPADIdata(
	input = "B14017", input.transform="diff", input.names="R90",
	output =   c("B820600", "I37026", "b1627", "b14013", "b4237", "D767608",
                      "b3400", "M.JQIND", "M.CUSA0"), # P484549=CPI discontinued
	output.names = c( "CPI", "GDP", "M1", "RL", "TSE300", "employment", 
                       "PFX", "US ind. prod.", "US CPI"),

	output.transforms = c("percent.change", "percent.change", 
                  "percent.change", "diff", "diff", "percent.change", 
                     "percent.change", "percent.change",  "percent.change"),
        server="ets", start.server=TRUE, server.process="fame.server", 
        cleanup.script="cleanup.fame.server", stop.on.error=TRUE, warn=TRUE)

#         c("ets","", "M.BCPI",  "percent.change", "com. price ind."),
#   availability(JofF.VAR.data.names)
      JofF.VAR.data <- freeze(JofF.VAR.data.names)

      z <- concentrate(JofF.VAR.data, p=2)
      ok <- all(seriesNamesOutput(reconstitute(z)) 
                   == paste("recon.", seriesNamesOutput(JofF.VAR.data))) 
      tfplot(z)
      tfplot(reconstitute(z))

      all.ok <- all.ok  & ok
      if (ok) cat("ok\n") else  cat("failed!\n") 

}


  if (! all.ok) stop("some tests FAILED")


