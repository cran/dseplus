# These set of tests will not work unless the PADI interface is
#   also installed (available at http://www.bank-banque-canada.ca/pgilbert)

   require("ts",      warn.conflicts=T)
   require("dse2",    warn.conflicts=T)
   require("dsepadi", warn.conflicts=T)
   require("padi",    warn.conflicts=T)

 Sys.info()
 version.dse()

   if (is.S()) {
	# the next 2 lines remove old versions of PADI in the search path
 	invisible(if(0!=length(grep("b*/PADI/.Data",search())))
                         detach(grep("b*/PADI/.Data",search()))  )
	attach(paste(getenv("PADI"),"/.Data", sep=""), pos=3)
	#load.padi(from=".")    # this gets the version in pwd
	load.padi()           # this gets the version indicated by PADI
	# load.padi does the following two dynamic loads 
	#dyn.load.shared("/usr/lib/libnsl.so")     # splus 3.3 on SunOS5
	# If the shared library is not loaded then the next has missing symbols
	#dyn.load(paste(getenv("PADI"),"/lib/splusclnt.o", sep=""))
   }

#    a TS PADI server is necessary for the following

 cat("search path ", search(),"\n")
 cat("PATH set to ",  Sys.getenv("PATH"), "\n")
 cat("PADI set to ",  Sys.getenv("PADI"), "\n")
 cat("PADI_STARTUP set to ", Sys.getenv("PADI_STARTUP"), "\n")
 cat("PADI_CLEANUP set to ", Sys.getenv("PADI_CLEANUP"), "\n")
 cat("user name set to ", Sys.info()[["user"]], "\n")
 

#######################################################################

#    tfPADI interface tests (from Brief User's Guide)   <<<<<<<<<<

#######################################################################


tfPADI.function.tests <- function( verbose=T, synopsis=T,
      fuzz.small=1e-14, fuzz.large=1e-6)
{# test for TSPADI access using simple.server

 # These tests only check that the tfPADI structures work. For a more
 #   complete set of PADI tests see the file padi.s distributed 
 #   with the TS PADI software.


  if (synopsis & !verbose) cat("tfPADI tests ...")

  scratch.db <-"zot123456.db"
  syskern.rm(scratch.db)
  server <- Sys.info()[["nodename"]]

 if (verbose) cat("tfPADI test 0 ... ")
  if (check.padi.server(server))
     stop("A server is already running. Testing stopped. Use cleanup.padi.server() or kill.padi.server() to terminate it.")

  pid <- start.padi.server(server=server, dbname="", 
                 server.process=paste("simple.server ", scratch.db))
  on.exit(cleanup.padi.server(pid, cleanup.script="cleanup.simple.server"))

  # wait to ensure padi server is started
     for (i in 1:30)
       {if (check.padi.server(server)) break
        Sys.sleep(1)
       }
  all.ok <- ok <- T
  if (verbose) 
    {if (ok) cat("ok\n")
     else  cat("failed! starting server\n")
    }


  if (verbose) cat("tfPADI test 1 ... ")

  eg.put.data <- tframed(matrix(c(1*exp(1:20),2*exp(1:20)),20,2), 
                         list(start=c(1950,1),freq=1))
  series.names(eg.put.data) <- c("exp1", "exp2")

  if (any(series.names(eg.put.data) != c("exp1", "exp2")))
    stop("series.name setting is not working properly. Other tests will fail.")

  eg.names <- tfputpadi(eg.put.data,
                      dbname=scratch.db, server=server,
                      start.server=T, server.process="simple.server", 
                      cleanup.script="cleanup.simple.server",
                      stop.on.error=T, warn=T )
  ok<-is.tfPADIdata(eg.names) 
  all.ok <- ok
  if (verbose) 
    {if (ok) cat("ok\n")
     else  cat("failed! tfputpadi\n")
    }

  if (verbose) cat("tfPADI test 2 ... ")
  eg.data <- freeze(eg.names)
  ok <- is.tfPADIdata(eg.names) &
            test.equal(eg.data, eg.put.data, fuzz=fuzz.large)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n")  else cat("failed!\n") }


  on.exit()
  cleanup.padi.server(pid, cleanup.script="cleanup.simple.server")

  if (synopsis) 
    {if (verbose) cat("All tfPADI tests completed")
     if (all.ok) cat(" OK\n") else cat(", some FAILED!\n")
    }

  if (all.ok) invisible(T)  else stop("FAILED")
}
 
 
#    TS PADI interface tests 

TSPADI.function.tests <- function( verbose=T, synopsis=T,
      fuzz.small=1e-14, fuzz.large=1e-6)
{# test for TSPADI access using simple.server

 # These tests only check that the DSE structures work with PADI. For a more
 #   complete set of PADI tests see the file padi.s distributed 
 #   with the TS PADI software.


  if (synopsis & !verbose) cat("DSE TSPADI tests ...")

  scratch.db <-"zot123456.db"
  syskern.rm(scratch.db)
  server <- Sys.info()[["nodename"]]

 if (verbose) cat("DSE TSPADI test 0 ... ")
  if (check.padi.server(server))
     stop("A server is already running. Testing stopped. Use cleanup.padi.server() or kill.padi.server() to terminate it.")

  pid <- start.padi.server(server=server, dbname="", 
                 server.process=paste("simple.server ", scratch.db))
  on.exit(cleanup.padi.server(pid, cleanup.script="cleanup.simple.server"))

  # wait to ensure padi server is started
     for (i in 1:30)
       {if (check.padi.server(server)) break
        Sys.sleep(1)
       }

  exp1 <- tframed(matrix(1*exp(1:20),20,1), list(start=c(1950,1),freq=1))
#  exp1 <- tframed(1*exp(1:20), list(start=c(1950,1),freq=1))
#  tframe(exp1) <- tframe(exp1)
  eg.put.data <- TSdata(input= exp1, 
                       output= tframed(tbind(2*exp1, 3*exp1),tframe(exp1)))
  series.names(eg.put.data) <- list(input="exp1", output=c("exp2","exp3"))

  if (any(input.series.names(eg.put.data) != "exp1"))
    stop("series.name setting is not working properly. Other tests will fail.")

  if (any(output.series.names(eg.put.data) != c("exp2","exp3")))
    stop("series.name setting is not working properly. Other tests will fail.")

#  exp1 <- tframed(1*exp(1:20), list(start=c(1950,1),freq=1))
#  eg.put.data <- list(input= tsmatrix(exp1), 
#                      input.names="exp1",
#                      output= tsmatrix(2*exp1, 3*exp1), 
#                      output.names=c("exp2","exp3"))
  eg.names <- putpadi.TSdata(eg.put.data,
                      dbname=scratch.db, server=server,
                      start.server=T, server.process="simple.server", 
                      cleanup.script="cleanup.simple.server",
                      stop.on.error=T, warn=T )
  ok<-is.TSPADIdata(eg.names) 
  all.ok <- ok
  if (verbose) 
    {if (ok) cat("ok\n")
     else  cat("failed! putpadi server started\n")
    }

  if (verbose) cat("DSE TSPADI test 1 ... ")
  eg.data <- freeze(eg.names)
  ok <- is.TSdata(eg.data ) & test.equal(eg.data, eg.put.data, fuzz=fuzz.large)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n")  else cat("failed!\n") }

  if (verbose) cat("DSE TSPADI test 2 ... ")

#If server= is supplied in the next, it should be "" and not NULL as previously 
eg.names <- TSPADIdata(input=c( "exp1","exp2"), output=c( "exp1","exp2","exp3"),
              frequency=1,
              db=scratch.db, stop.on.error=T, warn=T)

# z <- freeze(eg.names$input)
  eg.data <- freeze(eg.names)
  ok <- is.TSdata(eg.data ) 
warning("skipping something broken")
#&
#    (max(abs(output.data(eg.data) - 
#              cbind(exp(1:20),2*exp(1:20),3*exp(1:20)) ))<fuzz.large)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n")  else cat("failed!\n") }

  if (verbose) cat("DSE TSPADI test 3 ... ")
  avail <- availability(eg.names, verbose=F)
  ok <- all(c(avail$start ==  t(matrix(c(1950,1),2,5)),
              avail$end   ==  t(matrix(c(1969,1),2,5)),
              avail$frequency ==  rep(1,5)))

  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n")  else cat("failed!\n") }

  on.exit()
  cleanup.padi.server(pid, cleanup.script="cleanup.simple.server")

  if (synopsis) 
    {if (verbose) cat("All DSE TSPADI tests completed")
     if (all.ok) cat(" OK\n") else cat(", some FAILED!\n")
    }

  if (all.ok) invisible(T)  else stop("FAILED")
}


if ( ! require("padi") ) warning("Warning: package padi is needed.") else {
   Sys.sleep(5)
   tfPADI.function.tests()
   Sys.sleep(5)
   padi.function.tests.simple(verbose=T)     # all ok
   Sys.sleep(5)
   TSPADI.function.tests(verbose=T)   # all ok, 
   }

