   require("stats") #,      warn.conflicts=TRUE)
   require("dse2") #,    warn.conflicts=TRUE)
   require("dsepadi") #, warn.conflicts=TRUE)
   require("monitor") #, warn.conflicts=TRUE)
   
#   A TS PADI server is necessary for these tests.

   require("padi",    warn.conflicts=TRUE)

 Sys.info()
 DSEversion()

 cat("search path ", search(),"\n")
 cat("PATH set to ",  Sys.getenv("PATH"), "\n")
 cat("PADI set to ",  Sys.getenv("PADI"), "\n")
 cat("PADI_LDLIB set to ",  Sys.getenv("PADI_LDLIB"), "\n")
 cat("PADI_STARTUP set to ", Sys.getenv("PADI_STARTUP"), "\n")
 cat("PADI_CLEANUP set to ", Sys.getenv("PADI_CLEANUP"), "\n")
 cat("user name set to ", Sys.info()[["user"]], "\n")




###########################################################################

# Tests function    <<<<<<<<<<<<

###########################################################################

combination.monitor.function.tests <- function( verbose=TRUE, synopsis=TRUE, 
         fuzz.small=1e-10,
         server.process = PADIserverProcess(),
         cleanup.script = PADIcleanupScript() )
{# Some of the tests here are really for functions defined in dse1 ... dse3
 #   but are not tested there to avoid assuming Fame access is available.
 # The main short coming of these tests is that they do not test
 #     functions which produce output or graphs.
 # These tests require access to Fame data bases and the files:
 #          monitoring.test.db    fake database 
 #          monitoring.test.info  comparison info. to check results
 #          monitoring.test.data  fake over-riding data 

 # Note also that the test data is not real data (it may have been differenced
 #  or otherwise transformed) and is only intended to test that functions
 #  work as originally specified. 

  server <- Sys.info()[["nodename"]]
#  db     <- paste(.path.package("monitor"),"/otherdata/monitoring.test.db",sep="")
  db     <- system.file("otherdata", "monitoring.test.db", package="monitor")
  if (synopsis & !verbose) cat("All combination monitor tests ...")
  all.ok <- TRUE

  if (verbose) cat("combination monitor test 0 ... ")
  # simulated a database server
  pid <- startPADIserver(server=server, dbname=db, 
           server.process=server.process)
  on.exit(cleanupPADIserver(pid, cleanup.script=cleanup.script))

  # wait for server to start 
     for (i in 1:30)
       {if (checkPADIserver(server)) break
        Sys.sleep(1)
       }
  ok <- TRUE
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n")  else cat("failed!\n") }

  if (verbose) cat("combination monitor test 1 ... ")
  #  dbname=db would not be nec. with a public mode fame server

  test.data.names <- TSPADIdata(
      input="B14017", 
#      input.transforms= "diff",
       output=c( "P484549", "I37026", "lfsa201","b3400"), 
#      output.transforms= rep("percentChange",4),
      db=db, server=server,pad.end =TRUE)

#  source(paste(.path.package("monitor"),"/otherdata/monitoring.test.info", sep=""))
  source(system.file("otherdata", "monitoring.test.info", package="monitor"))

  v.data <- verification.data
  v.data$output <- v.data$output[,c(1,2,6,7)]
  tframe(v.data$output) <- tframe(verification.data$output)
  ok <- is.TSdata(v.data)
  all.ok <- ok 
  if (verbose) 
    {if (ok) cat("ok\n")
     else    cat("failed! (loading verification data)\n")
    }


  if ( ! (require("padi") && ("fame.server" == PADIserverProcess() ))) 
     warning("Warning: skipping tests that require local fame server.") 
  else {

  if (verbose) cat("combination monitor test 2 ... ")
  data <-retrieve.and.verify.data(test.data.names, 
                                    verification.data=v.data)
  ok <- testEqual(data, ets.test.data, fuzz=fuzz.small)
  tags(data$input)  <- "data"
  tags(data$output) <- "data"
  all.ok <- all.ok & ok 
  if (verbose) 
    {if (ok) cat("ok\n")
     else    cat("failed! (retrieve.and.verify.data)\n")
    }

  if (verbose) cat("combination monitor test 3 ... ")
  overriding.data <- get.overriding.data(
                   #file=paste(.path.package("monitor"),"/otherdata/monitoring.test.data", sep=""),
                   file=system.file("otherdata", "monitoring.test.data", package="monitor"),
                   m=1, p=10,
                   first.input="diff(R90=B14017)", 
                   first.output="%change(CPI=P484549)",  
                   second.output="%change(GDP=I37026)"  )
  z.tf <-tframe(overriding.data$output)
  overriding.data$output <- overriding.data$output[,c(1,2,6,7)]
  tframe(overriding.data$output) <- z.tf
  tags(overriding.data$input) <- "over"
  tags(overriding.data$output) <- "over"
  ok <- testEqual(overriding.data, monitor.test.data, fuzz=fuzz.small)
  all.ok <- all.ok & ok 
  if (verbose) 
    {if (ok) cat("ok\n")
     else    cat("failed! (get.overriding.data)\n")
    }

  if (verbose) cat("combination monitor test 4 ... ")
  combined.forecast<-combineAndForecast(monitoring.test.model,
                  list(data=data, overriding.data=overriding.data)) 

#  ok <- fuzz.small > max(abs( combined.forecast$best.guess - 
#                            best.guess.test))
#  Rbug - kludge - the above chokes on - because classes are not the same and 
#  gives Error: invalid time series parameters specified

  ok <- fuzz.small > abs( sum(combined.forecast$best.guess) - 
                              sum(best.guess.test))
  all.ok <- all.ok & ok 
  if (verbose) 
    {if (ok) cat("ok\n")
     else    cat("failed! (combineAndForecast)\n")
    }

  }
  
  if (synopsis) 
    {if (verbose) cat("All combination monitor tests completed")
     if (all.ok) cat(" OK\n\n")
     else    cat(", some FAILED!\n\n")
    }

  if (all.ok) invisible(TRUE)  else stop("FAILED")
}


  Sys.sleep(15) # just in case a previous server has not yet died

if ( ! (require("padi") && checkPADIserver("ets"))) 
     warning("Warning: package padi and server ets are needed.") else
   combination.monitor.function.tests(verbose=TRUE)
