 require("padi")
 Sys.info()
 Sys.getenv()["PATH"]
 Sys.getenv()["PADI"]
 Sys.getenv()["PADI_LDLIB"]
 PADIserverProcess()
 
padi.function.tests.simple <- function( verbose=TRUE, synopsis=TRUE, fuzz=1e-6, 
          scratch.db="zot123456.db")
{# A short set of tests for the S TS PADI client interface
  cat("Assuming load.padi has been executed (eg. in .First.lib)\n")
# simple.server checks that the basedir (given at startup) is contained in
#   the dbname (passed by requests) so a fully qualified path works.

#  Note: Sys.getenv()["PWD"] does not work if the run.tests script is run from a sh
#   script in another directory if cd is aliased, so the following is better:
  pwd <- getwd() # previously present.working.directory()
  fqpscratch.db <- paste(pwd,"/",scratch.db, sep="")
  unlink(scratch.db, recursive = TRUE)

  server <- Sys.info()[["nodename"]]
  if (checkPADIserver(server))
     stop("A server is already running. Testing stopped. Use cleanupPADIserver() or killPADIserver() to terminate it.")

  wait.for.server.to.terminate <- function(server)
    {# wait to ensure padi server is terminated
     for (i in 1:30)
       {if (!checkPADIserver(server)) break
        Sys.sleep(1)
       }
    }

  wait.for.server.to.start <- function(server)
    {# wait to ensure padi server is started
     for (i in 1:30)
       {if (checkPADIserver(server)) break
        Sys.sleep(1)
       }
    }


#    simple.server tests

  if (synopsis & !verbose) cat("All S to simple.server padi tests ...")
  all.ok <-TRUE
  max.error <- NA
  if (verbose) cat("S to simple.server padi test 1a... ")
  ok <- putpadi(ts(exp(1:20), start=c(1950,1),freq=1), series="exp", 
            server=server, server.process="simple.server",
            cleanup.script="cleanup.simple.server",
            dbname=scratch.db, start.server=TRUE)
  all.ok <- all.ok & ok
  if (verbose) 
    {if (ok) cat("ok\n")
     else  cat("failed! putpadi and starting server\n")
    }

  wait.for.server.to.terminate(server)

  if (verbose) cat("S to simple.server padi test 1b... ")
  ok <- putpadi(ts(exp(1:20), start=c(1950,1),freq=1), series="exp", 
            server=server,  server.process="simple.server",
            cleanup.script="cleanup.simple.server",
            dbname=fqpscratch.db, start.server=TRUE)
  all.ok <- all.ok & ok
  if (verbose) 
    {if (ok) cat("ok\n")
     else  cat("failed! putpadi and starting server\n")
    }

  wait.for.server.to.terminate(server)

  if (verbose) cat("S to simple.server padi test 2 ... ")
  z <- getpadi("exp",server=server,  server.process="simple.server",
                cleanup.script="cleanup.simple.server",
                dbname=scratch.db, start.server=TRUE)
  error <- max(abs((z - exp(1:20))))
  ok <-  fuzz > error
  ok <- ok & (frequency(z) == 1) & all(start(z) == c(1950,1))
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok
  if (verbose) 
    {if (ok) cat("ok\n")
     else  cat("failed! getpadi and starting server\n")
    }

  wait.for.server.to.terminate(server)

  #       simulated shared database ie. a public mode server
  # (this is not the way one would normally set up a public mode server)

 if (verbose) cat("S to simple.server padi test 3 ... ")
  pid <- startPADIserver(server=server, dbname="", 
                 server.process=paste("simple.server ", scratch.db))
  on.exit(cleanupPADIserver(pid, cleanup.script="cleanup.simple.server"))

  wait.for.server.to.start(server)

  ok <- putpadi(ts(2*exp(1:20), start=c(1950,1),freq=1), series="exp2", 
            server=server, dbname=scratch.db, start.server=FALSE)
  all.ok <- all.ok & ok
  if (verbose) 
    {if (ok) cat("ok\n")
     else  cat("failed! putpadi server started\n")
    }

  if (verbose) cat("S to simple.server padi test 4 ... ")
  z <- getpadi("exp2",server=server,  server.process="simple.server",
                  cleanup.script="cleanup.simple.server",
                  dbname=scratch.db, start.server=FALSE)
  error <- max(abs((z - 2*exp(1:20))))
  ok <-  fuzz > error
  ok <- ok & (frequency(z) == 1) & all(start(z) == c(1950,1))
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok
  if (verbose) 
    {if (ok) cat("ok\n")
     else  cat("failed! getpadi server started\n")
    }

  cleanupPADIserver(pid, cleanup.script="cleanup.simple.server")
  on.exit()
  unlink(scratch.db, recursive = TRUE) # this is a directory
#  system(paste("rm -r ", scratch.db))

  if (synopsis) 
    {if (verbose) cat("All S to simple.server padi tests completed")
     if (all.ok) cat(" ok\n\n")
     else    cat(", some failed! (max. error magnitude= ", max.error,")\n")
    }

  wait.for.server.to.terminate(server)
  invisible(all.ok)
}


   Sys.sleep(5)
   padi.function.tests.simple(verbose=TRUE)     # all ok
