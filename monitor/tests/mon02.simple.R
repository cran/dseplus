   require("ts")
   require("dse2")
   require("dsepadi")
   require("monitor")
#   A TS PADI server is necessary for these tests.
#   The next line is only necessary to remove this in an old version which set
#     home in frame 0. (I'll never do that again.)
#   if (is.S()) remove("DSE.HOME", where=0) 


 Sys.info()
 version.dse()

 cat("search path ", search(),"\n")
 cat("PATH set to ",  Sys.getenv("PATH"), "\n")
 cat("PADI set to ",  Sys.getenv("PADI"), "\n")
 cat("PADI_STARTUP set to ", Sys.getenv("PADI_STARTUP"), "\n")
 cat("PADI_CLEANUP set to ", Sys.getenv("PADI_CLEANUP"), "\n")
 cat("user name set to ", Sys.info()[["user"]], "\n")

   require("padi")
   if (is.S()) {
	# the next 2 lines remove old versions of PADI in the search path
 	invisible(if(0!=length(grep("b*/PADI/.Data",search())))
                         detach(grep("b*/PADI/.Data",search()))  )
	attach(paste(Sys.getenv("PADI"),"/.Data", sep=""), pos=2)
	#load.padi(from=".")    # this gets the version in pwd
	load.padi()           # this gets the version indicated by PADI
	# load.padi does the following two dynamic loads 
	#dyn.load.shared("/usr/lib/libnsl.so")     # splus 3.3 on SunOS5
	# If the shared library is not loaded then the next has missing symbols
	#dyn.load(paste(Sys.getenv("PADI"),"/lib/splusclnt.o", sep=""))
	search()
   }



###########################################################################

# Tests function for data retrieval for simple monitoring    <<<<<<<<<<<<

###########################################################################


simple.monitor.function.tests <- function( verbose=T, synopsis=T, 
         fuzz.small=1e-14, fuzz.large=1e-8,
         server.process = padi.server.process(),
         cleanup.script = padi.cleanup.script() )
{# Some of the tests here are really for functions defined in dse1 ... dse3
 #   but are not tested there to avoid assuming TSPADI (or Fame) access is
 # available. The main short coming of these tests is that they do not test
 #     functions which produce output or graphs.
 # These tests require access to Fame data bases and the files:
 #          monitoring.test.db    fake database 
 #          monitoring.test.info  comparison info. to check results

 # Note also that the test data is not real data (it may have been differenced
 #  or otherwise transformed) and is only intended to test that functions
 #  work as originally specified. 

  server <- Sys.info()[["nodename"]]
  db     <- paste(DSE.HOME,"/data/monitoring.test.db",sep="")

  if (synopsis & !verbose) cat("All simple monitor tests ...")
  all.ok <- T

  if (verbose) cat("simple monitor test 0 ... ")
  # simulate a database server
  pid <- start.padi.server(server=server,
           dbname=db, 
           server.process=server.process)
  on.exit(cleanup.padi.server(pid, cleanup.script=cleanup.script))

  # wait for server to start 
     for (i in 1:30)
       {if (check.padi.server(server)) break
        Sys.sleep(1)
       }
  ok <- T
  all.ok <- all.ok & ok 
  if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }


  if (verbose) cat("simple monitor test 1 ... ")
  #  db=db would not be nec. with a public mode fame server   
  test.data.names <- TSPADIdata(
      input  ="B14017", 
      output = c( "P484549", "I37026", "lfsa201","b3400"), 
      server=server, db=db, pad.end =T)
   
  z <-availability(test.data.names, verbose=F) 
  ok <- all(c(z$start == t(matrix(c(1974,2),2,5)), 
              z$end   == t(matrix(c(1993,9),2,5)), 
              z$freq==rep(12,5) ))
  all.ok <- all.ok & ok 
  if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }


# the following sets ets.test.data, monitor.test.data, verification.data
#      and  monitoring.test
  source(paste(DSE.HOME,"/data/monitoring.test.info", sep=""))

  if (verbose) cat("simple monitor test 2 ... ") 
  v.data <- verification.data
  output.data(v.data) <- output.data(v.data)[,c(1,2,6,7)]
  tframe(output.data(v.data)) <- tframe(output.data(verification.data))
  ok <- is.TSdata(v.data)
  all.ok <- all.ok & ok 
  if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }

  if (verbose) cat("simple monitor test 3 ... ")
  hist.data <-retrieve.and.verify.data(test.data.names, 
                                    verification.data=v.data)
  ok <- test.equal(hist.data, ets.test.data, fuzz=fuzz.small)
  all.ok <- all.ok & ok 
  if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }


  if (verbose) cat("simple monitor test 4 ... ")
  monitoring<-simple.monitoring (monitoring.test.model, test.data.names, 
        previous.data=NULL, mail.list=Sys.info()[["user"]], error.mail.list=Sys.info()[["user"]]) 
  ok <-  monitoring$status == "Simple monitoring initialized."   
  if (verbose) cat("\n This test produces a warning: Input is not longer than output data. No forecasts produced...")
  # note that the following does not result in forecasts (and the forecast
  #   function produces a warning) because the input data does not extend
  #   beyond the output data.
  monitoring<-simple.monitoring (monitoring.test.model, test.data.names, 
           previous.data=monitoring$data, 
	   mail.list=Sys.info()[["user"]], error.mail.list=Sys.info()[["user"]]) 
  ok <- ok & (monitoring$status == "Simple monitoring updates not necessary.")
  monitoring<-simple.monitoring (monitoring.test.model, test.data.names, 
               previous.data=monitoring$data, 
               mail.list=Sys.info()[["user"]], error.mail.list=Sys.info()[["user"]], run.again=T) 
  ok <- ok & (monitoring$status == "Simple monitoring re-run.")
  ok <- ok & monitoring$message[7] == 
          "1993 Sep   0.110000   0.383440   0.397520   0.355500   0.947460 "
  ok <- ok & sum(output.data(monitoring$data))==235.64806565791809589
  output.data(monitoring$data) <- 
               tfwindow(output.data(monitoring$data), end=c(1993,8))
  monitoring<-simple.monitoring (monitoring.test.model, test.data.names, 
          previous.data=monitoring$data, 
	  mail.list=Sys.info()[["user"]], error.mail.list=Sys.info()[["user"]]) 
  ok <- ok & (monitoring$status == "Simple monitoring updated.") &
      sum(output.data(monitoring$data)) == 235.64806565791809589
  all.ok <- all.ok & ok 
  if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }


  if (verbose) cat("simple monitor test 5 ... ")

  watch <- watch.data(test.data.names, previous.data=NULL, mail.list=Sys.info()[["user"]])
  ok <- (watch$status == "System watch.data initialized.") & 
         sum(output.data(watch$data))== 235.64806565791809589
  watch <- watch.data(test.data.names, previous.data=watch, mail.list=Sys.info()[["user"]])
  ok <- ok & (watch$status == "No data updates.") & 
           sum(input.data(watch$data))== -4.1300000572204575988
  watch$data <- tfwindow(watch$data, end=c(1993, 8))
  watch <- watch.data(test.data.names, previous.data=watch, mail.list=Sys.info()[["user"]])
  ok <- ok & (watch$status == "Data has been updated.") & 
          sum(output.data(watch$data))== 235.64806565791809589

  all.ok <- all.ok & ok 
  if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }

  if (synopsis) 
    {if (verbose) cat("All simple monitor tests completed")
     if (all.ok) cat(" OK\n\n") else    cat(", some FAILED!\n\n")
    }

  if (all.ok) invisible(T)  else stop("FAILED")
}


  Sys.sleep(15) # just in case a previous server has not yet died

if ( ! require("padi") ) warning("Warning: package padi is needed.") else
   simple.monitor.function.tests(verbose=T) 

