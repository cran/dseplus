# Expecting results from  tests

#R --nsize=1M --vsize=6M 
  require("mva"); require("ts"); require("dse2"); require("juice") # adds dse, tframe, and syskern
 #x11()
 # postscript(file="juice.out.ps",  paper="letter", horizontal=F, onefile=T)
 #           # width=6, height=8, pointsize=10,





######################################################################

#  test function

######################################################################


juice.function.tests <- function( verbose=TRUE, synopsis=TRUE, fuzz.small=1e-12,
                                  graphics=TRUE, ets=F)
{
  all.ok <-  T
  if (synopsis & !verbose) cat("All dse juice tests ...")
  if      (is.R()) data("egJofF.1dec93.data", package="dse1")
  else if (is.S()) source(paste(DSE.HOME, "/data/egJofF.1dec93.data.R", sep=""))

  if (is.S()) test.rng <- list(kind="default", normal.kind="default", 
                       seed=c(13,44,1,25,56,0,6,33,22,13,13,0))
  if (is.R())
    test.rng <- list(kind="default", normal.kind="default",
                     seed=c( 979)) #, 1479, 1542))

  set.RNG(test.rng)
  data00 <- matrix(rnorm(300), 100,3)
  data0 <- TSdata(output=data00)
  data1 <- TSdata(output=data0$output %*% matrix(rnorm(9),3,3))


  if (verbose) cat("dse juice test 0 ... ")
  z <- concentrate(data00, conc=est.projection(data00, n=3))
  # next is true because all PCs are used
  ok <-      test.equal(data00, reconstitute(z), fuzz=fuzz.small)
# T but does this have meaning  ok <- ok & test.equal(data00,reconstitute(z), fuzz=fuzz.small)
  ok <- ok & all(series.names(reconstitute(z))
                == paste("recon.", series.names(data00)))
  ok <- ok & all(concentrated.series.names(z)
                == paste("concentrate", seq(length=concentrated.dimension(z))))
  ok <- ok & all(series.names(z) == series.names(data00))		   		   

  z <- concentrate(output.data(egJofF.1dec93.data)) # p=1
  z <-tfwindow(z, start=c(1981,1))
  ok <- all(c(1981,1) == start(z)) & all(c(1981,1) == start(reconstitute(z)))
  all.ok <- all.ok & ok
  if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }

  if (verbose) cat("dse juice test 1 ... ")
  z <- concentrate(data0, conc=est.projection(data0, p=3))
  # next is true because all PCs are used
  ok <-      test.equal(data0, reconstitute(z), fuzz=fuzz.small)
# T but does this have meaning  ok <- ok & test.equal(data0,TSdata(z), fuzz=fuzz.small)
  ok <- ok & all(output.series.names(reconstitute(z))
                == paste("recon.", output.series.names(data0)))
  ok <- ok & all(concentrated.output.series.names(z)
                == paste("concentrate", seq(length=concentrated.output.dimension(z))))
  ok <- ok & all(output.series.names(z) == output.series.names(data0))		   		   
  all.ok <- all.ok & ok 
  if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }

  if (verbose) cat("dse juice test 2 ... ")
  z <- concentrate(data0, conc=est.projection(data0, p=1))
  ok <-       all(output.series.names(data0) == output.series.names(z))
  ok <- ok & "concentrate 1" == concentrated.output.series.names(z)
  ok <- ok &  ! test.equal(data0, reconstitute(z), fuzz=fuzz.small)
  all.ok <- all.ok & ok 
  if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }

  if (verbose) cat("dse juice test 3 ... ")
  z <- concentrate(data0, conc=est.projection(data0, p=2))
  z <- concentrate(data0, conc=est.projection(data0)) # p=1

  ok <- test.equal(z,   concentrate(data0,conc=z), fuzz=fuzz.small)
  ok <- ok & (3 == output.dimension(z)) & (3 == output.dimension(TSdata(z))) &
             (1 == concentrated.output.dimension(z))               &
             (100 == periods(z)) & (100 == periods(TSdata(z)))
  all.ok <- ok
  if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }

  if (verbose) cat("dse juice test 4 ... ")
  z <- concentrate(egJofF.1dec93.data) # p=1
  ok <- all(c(1974,2) == start(z)) & all(c(1974,2) == start(reconstitute(z)))
  all.ok <- ok
  if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }

  if (verbose) cat("dse juice test 5 ... ")
  z <-tfwindow(z, start=c(1981,1))
  ok <- all(c(1981,1) == start(z)) & all(c(1981,1) == start(reconstitute(z)))
  all.ok <- ok
  if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }

  if (verbose) cat("dse juice test 6 ... ")
#  zm <- est.concentrated.model(z, scale=T, center=T, estimation="bft",
#  z is already a concentrated object so center and scale are not used
  zm <- est.concentrated.model(z, estimation="bft",
                                estimation.args=list(max.lag=2, verbose=F))
  #tfplot(zm)
  zmr <- check.residuals(zm, ac=F, pac=F, plot.=F)
# S values
  ok <-     all(zmr$skewness - 
       c(-1.979814544376786,  0.4241402252105713, -0.2758381743875303,
         -0.5754191393269468, 1.983488003995619,   0.1565333254836311,
         -0.116925097663112,  0.1326786266039132, -0.2649736832254932,
         -0.3710824450042611) < fuzz.small)
  ok <- ok & all(zmr$kurtosis -
       c(11.01131667332406,   3.503769190958172, 5.024697760083122,
          4.972991895880567, 15.5906633209087,   2.3188280244819,
          2.935474704110169,  4.131056248843817, 3.823872187369519,
          4.900819828268634) < fuzz.small)
# print(zmr$skewness, digits=16)
# print(zmr$kurtosis, digits=16)
  zmf <- featherForecasts(zm)
#  ok <- ???
  all.ok <- ok
  if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }

  if (ets)
    {if (verbose) cat("  juice test using ets ...")
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
        cleanup.script="cleanup.fame.server", stop.on.error=TRUE, warn=T )


#         c("ets","", "M.BCPI",  "percent.change", "com. price ind."),

#   availability(JofF.VAR.data.names)
      JofF.VAR.data <- freeze(JofF.VAR.data.names)

      z <- concentrate(JofF.VAR.data, p=2)
      z <- concentrate(JofF.VAR.data, p=3)
      z <- concentrate(JofF.VAR.data, p=2, scale=F)
      ok <- test.equal(series.names(reconstitute(z)), 
                       series.names(JofF.VAR.data)) 

}

  if (synopsis) 
    {if (verbose) cat("All dse juice tests completed")
     if (all.ok) cat(" OK\n") else    cat(", some FAILED!\n")
    }
    
  if (graphics)
    all.ok <- all.ok & juice.graphics.tests(verbose=verbose, synopsis=synopsis,
     				 ets=ets)

  if (all.ok) invisible(T)  else stop("FAILED")
}

juice.graphics.tests <- function( verbose=TRUE, synopsis=TRUE, ets=F)
  {# graphics tests do not do any value comparisons

# R
# library("tframe"); library("syskern"); library("padi"); 
# library("dse"); library("dsex1")

# S
# attach(paste(Sys.getenv("PADI"),".Data", sep="/"), first=TRUE)
# attach("/home/res/gilp/dse/pub/Sdse/.Data", first=TRUE)

#   source("/home/res/gilp/dse/my/src/personal.utils.s")
#  hsource("/home/res/gilp/dse/my/src/juice.hs")

    if (synopsis & !verbose)  cat("juice graphics tests ...")
    # If no device is active then write to postscript file 
    if ( dev.cur() == 1 ) {
                postscript(file = "zot.postscript.test.ps", width = 6, 
                        height = 6, pointsize = 10, onefile = F, 
                        print.it = F, append = F)
                on.exit((function() {
                        dev.off()
                        synchronize(1)
                        rm("zot.postscript.test.ps")
                }))
        }

    all.ok <-  T
   if (is.S()) test.rng <- list(kind="default", normal.kind="default", 
                       seed=c(13,44,1,25,56,0,6,33,22,13,13,0))
   if (is.R()) test.rng <- list(kind="default", normal.kind="default",
                       seed=c( 979)) #, 1479, 1542))
    set.RNG(test.rng)
    data0 <- TSdata(output=matrix(rnorm(300), 100,3))
    series.names(data0)<- 
                  list(output=paste("data0", seq(length=output.dimension(data0))))
    data1 <- TSdata(output=data0$output %*% matrix(rnorm(9),3,3))
    series.names(data1)<- 
                  list(output=paste("data1", seq(length=output.dimension(data1))))

    if (verbose) cat("  juice graphics test 1 ...")
    z <- concentrate(data0) # p=1
    ok <-  all( output.series.names(z) 
                    == paste("data0", seq(output.dimension(z))))

    concentrated.tfplot(z)
    tfplot(z)
    all.ok <- all.ok  & ok
    if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }

    if (verbose) cat("  juice graphics test 2 ...")
    tfplot(reconstitute(z))
    ok <-  all( output.series.names(reconstitute(z))
                      == paste("recon. data0", seq(output.dimension(data0))))
    all.ok <- all.ok  & ok
    if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }

    if (verbose) cat("  juice graphics test 3 ...")
    z <- concentrate(data0, p=2)
    tfplot(z)
    all.ok <- all.ok  & ok
    if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }

    if (verbose) cat("  juice graphics test 4 ...")
    tfplot(reconstitute(z))
    all.ok <- all.ok  & ok
    if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }

    if (verbose) cat("  juice graphics test 5 ...")
    z <- concentrate(data0, p=3)
    tfplot(z)
    all.ok <- all.ok  & ok
    if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }
    if (verbose) cat("   (note that actual and reconstructed coincide.)\n")

    if (verbose) cat("  juice graphics test 6 ...")
    tfplot(reconstitute(z))
    all.ok <- all.ok  & ok
    if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }

    if (verbose) cat("  juice graphics test 7 ...")
    z <- concentrate(egJofF.1dec93.data, p=2)
    tfplot(z)
    z <- tfwindow(z, start=c(1992,1))
    tfplot(z)
    tfplot(reconstitute(z))
    all.ok <- all.ok  & ok
    if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }

  if (ets)
    {if (verbose) cat("  juice graphics test 7 using ets ...")
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
        cleanup.script="cleanup.fame.server", stop.on.error=TRUE, warn=T )

#         c("ets","", "M.BCPI",  "percent.change", "com. price ind."),
#   availability(JofF.VAR.data.names)
      JofF.VAR.data <- freeze(JofF.VAR.data.names)

      z <- concentrate(JofF.VAR.data, p=2)
      ok <- all(output.series.names(reconstitute(z)) 
                   == paste("recon.", output.series.names(JofF.VAR.data))) 
      tfplot(z)
      tfplot(reconstitute(z))

      all.ok <- all.ok  & ok
      if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }

}


  if (synopsis) 
    {if (verbose) cat("All juice graphics tests completed")
     if (all.ok) cat(" OK\n") else    cat(", some FAILED!\n")
    }
  if (all.ok) invisible(T)  else stop("FAILED")
}



   juice.function.tests( verbose=F, graphics=T)

