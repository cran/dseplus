# These  tests are skipped unless the PADI interface is
#   also installed (available at http://www.bank-banque-canada.ca/pgilbert)
#   and they fail it it is installed but not working.

   require("stats",      warn.conflicts=TRUE)
   require("dse2",    warn.conflicts=TRUE)
   require("dsepadi", warn.conflicts=TRUE)
   require("padi",    warn.conflicts=TRUE)

 Sys.info()
 DSEversion()

#   if (is.S()) {
#	# the next 2 lines remove old versions of PADI in the search path
# 	invisible(if(0!=length(grep("b*/PADI/.Data",search())))
#                         detach(grep("b*/PADI/.Data",search()))  )
#	load.padi()           # this gets the version indicated by PADI_LDLIB
#   }

#    a TS PADI server is necessary for the following

 cat("search path ", search(),"\n")
 cat("PATH set to ",  Sys.getenv("PATH"), "\n")
 cat("PADI set to ",  Sys.getenv("PADI"), "\n")
 cat("PADI_LDLIB set to ",  Sys.getenv("PADI_LDLIB"), "\n")
 cat("PADI_STARTUP set to ", Sys.getenv("PADI_STARTUP"), "\n")
 cat("PADI_CLEANUP set to ", Sys.getenv("PADI_CLEANUP"), "\n")
 cat("user name set to ", Sys.info()[["user"]], "\n")
 

#######################################################################

#    tfPADI interface tests (from Brief User's Guide)   <<<<<<<<<<

#######################################################################


tfPADI.function.tests <- function( verbose=TRUE, synopsis=TRUE,
      fuzz.small=1e-14, fuzz.large=1e-6)
{# test for TSPADI access using simple.server

 # These tests only check that the tfPADI structures work. For a more
 #   complete set of PADI tests see the file padi.s distributed 
 #   with the TS PADI software.


  if (synopsis & !verbose) cat("tfPADI tests ...")

  scratch.db <-"zot123456.db"
  unlink(scratch.db, recursive = TRUE)
  server <- Sys.info()[["nodename"]]

 if (verbose) cat("tfPADI test 0 ... ")
  if (checkPADIserver(server))
     stop("A server is already running. Testing stopped. Use cleanupPADIserver() or killPADIserver() to terminate it.")

  pid <- startPADIserver(server=server, dbname="", 
                 server.process=paste("simple.server ", scratch.db))
  on.exit(cleanupPADIserver(pid, cleanup.script="cleanup.simple.server"))

  # wait to ensure padi server is started
     for (i in 1:30)
       {if (checkPADIserver(server)) break
        Sys.sleep(1)
       }
  all.ok <- ok <- TRUE
  if (verbose) 
    {if (ok) cat("ok\n")
     else  cat("failed! starting server\n")
    }


  if (verbose) cat("tfPADI test 1 ... ")

  eg.put.data <- tframed(matrix(c(1*exp(1:20),2*exp(1:20)),20,2), 
                         list(start=c(1950,1),freq=1))
  seriesNames(eg.put.data) <- c("exp1", "exp2")

  if (any(seriesNames(eg.put.data) != c("exp1", "exp2")))
    stop("series.name setting is not working properly. Other tests will fail.")

  eg.names <- tfputpadi(eg.put.data,
                      dbname=scratch.db, server=server,
                      start.server=TRUE, server.process="simple.server", 
                      cleanup.script="cleanup.simple.server",
                      stop.on.error=TRUE, warn=TRUE)
  ok<-is.tfPADIdata(eg.names) 
  all.ok <- ok
  if (verbose) 
    {if (ok) cat("ok\n")
     else  cat("failed! tfputpadi\n")
    }

  if (verbose) cat("tfPADI test 2 ... ")
  eg.data <- freeze(eg.names)
  ok <- is.tfPADIdata(eg.names) &
            testEqual(eg.data, eg.put.data, fuzz=fuzz.large)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n")  else cat("failed!\n") }


  on.exit()
  cleanupPADIserver(pid, cleanup.script="cleanup.simple.server")

  if (synopsis) 
    {if (verbose) cat("All tfPADI tests completed")
     if (all.ok) cat(" OK\n") else cat(", some FAILED!\n")
    }

  if (!all.ok) stop("FAILED")
  # database is left in place to help debug if there is an error
  unlink(scratch.db, recursive = TRUE)
  invisible(TRUE)
}
 
 
#    TS PADI interface tests 

TSPADI.function.tests <- function( verbose=TRUE, synopsis=TRUE,
      fuzz.small=1e-14, fuzz.large=1e-6)
{# test for TSPADI access using simple.server

 # These tests only check that the DSE structures work with PADI. For a more
 #   complete set of PADI tests see the file padi.s distributed 
 #   with the TS PADI software.


  if (synopsis & !verbose) cat("DSE TSPADI tests ...")

  scratch.db <-"zot123456.db"
  unlink(scratch.db, recursive = TRUE)
  server <- Sys.info()[["nodename"]]

 if (verbose) cat("DSE TSPADI test 0 ... ")
  if (checkPADIserver(server))
     stop("A server is already running. Testing stopped. Use cleanupPADIserver() or killPADIserver() to terminate it.")

  pid <- startPADIserver(server=server, dbname="", 
                 server.process=paste("simple.server ", scratch.db))
  on.exit(cleanupPADIserver(pid, cleanup.script="cleanup.simple.server"))
  
  # wait to ensure padi server is started
     for (i in 1:30)
       {if (checkPADIserver(server)) break
        Sys.sleep(1)
       }

  exp1 <- tframed(matrix(1*exp(1:20),20,1), list(start=c(1950,1),freq=1))
#  exp1 <- tframed(1*exp(1:20), list(start=c(1950,1),freq=1))
#  tframe(exp1) <- tframe(exp1)
  eg.put.data <- TSdata(input= exp1, 
                       output= tframed(tbind(2*exp1, 3*exp1),tframe(exp1)))
  seriesNames(eg.put.data) <- list(input="exp1", output=c("exp2","exp3"))

  if (any(seriesNamesInput(eg.put.data) != "exp1"))
    stop("series.name setting is not working properly. Other tests will fail.")

  if (any(seriesNamesOutput(eg.put.data) != c("exp2","exp3")))
    stop("series.name setting is not working properly. Other tests will fail.")

#  exp1 <- tframed(1*exp(1:20), list(start=c(1950,1),freq=1))
#  eg.put.data <- list(input= tsmatrix(exp1), 
#                      input.names="exp1",
#                      output= tsmatrix(2*exp1, 3*exp1), 
#                      output.names=c("exp2","exp3"))
  eg.names <- putpadi.TSdata(eg.put.data,
                      dbname=scratch.db, server=server,
                      start.server=TRUE, server.process="simple.server", 
                      cleanup.script="cleanup.simple.server",
                      stop.on.error=TRUE, warn=TRUE)
  ok<-is.TSPADIdata(eg.names) 
  all.ok <- ok
  if (verbose) 
    {if (ok) cat("ok\n")
     else  cat("failed! putpadi server started\n")
    }

  if (verbose) cat("DSE TSPADI test 1 ... ")
  eg.data <- freeze(eg.names)
  ok <- is.TSdata(eg.data ) & testEqual(eg.data, eg.put.data, fuzz=fuzz.large)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n")  else cat("failed!\n") }

  if (verbose) cat("DSE TSPADI test 2 ... ")

#If server= is supplied in the next, it should be "" and not NULL as previously 
eg.names <- TSPADIdata(input=c( "exp1","exp2"), output=c( "exp1","exp2","exp3"),
              frequency=1,
              db=scratch.db, stop.on.error=TRUE, warn=TRUE)

# z <- freeze(eg.names$input)
  eg.data <- freeze(eg.names)
  ok <- is.TSdata(eg.data ) 
warning("skipping something broken")
#&
#    (max(abs(outputData(eg.data) - 
#              cbind(exp(1:20),2*exp(1:20),3*exp(1:20)) ))<fuzz.large)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n")  else cat("failed!\n") }

  if (verbose) cat("DSE TSPADI test 3 ... ")
  avail <- availability(eg.names, verbose=FALSE)
  ok <- all(c(avail$start ==  t(matrix(c(1950,1),2,5)),
              avail$end   ==  t(matrix(c(1969,1),2,5)),
              avail$frequency ==  rep(1,5)))

  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n")  else cat("failed!\n") }

  on.exit()
  cleanupPADIserver(pid, cleanup.script="cleanup.simple.server")

  if (synopsis) 
    {if (verbose) cat("All DSE TSPADI tests completed")
     if (all.ok) cat(" OK\n") else cat(", some FAILED!\n")
    }

  if (!all.ok) stop("FAILED")
  # database is left in place to help debug if there is an error
  unlink(scratch.db, recursive = TRUE)
  invisible(TRUE)
}


if ( ! require("padi") ) warning("Warning: package padi is needed.") else {
   Sys.sleep(5)
   tfPADI.function.tests()
   Sys.sleep(5)
   TSPADI.function.tests(verbose=TRUE)   # all ok, 
   }

