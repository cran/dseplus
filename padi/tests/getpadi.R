#  This is really a test of the getpadi standalone command, not the R
#   interface, but this is a convenient way to test. 

  require("dsepadi") # only to write data
  require("padi") # only to write data
  Sys.info()
  Sys.getenv()["PATH"]
  Sys.getenv()["PADI"]

  pwd <- getwd() # previously present.working.directory()
  scratch.db <- "zot123456b.db"
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

  ok <- putpadi(ts(exp(1:20), start=c(1950,1),freq=1), series="exp", 
            server=server, server.process="simple.server",
            cleanup.script="cleanup.simple.server",
            dbname=scratch.db, start.server=TRUE)

  if (ok) cat("ok\n") else  cat("failed! putpadi and starting server\n")

  wait.for.server.to.terminate(server)

  cat("getpadi  test simple ... ")
  
  pid <- system(paste("simple.server ", scratch.db), intern = TRUE)
  pid
  z <- system(paste("getpadi ", server," exp", scratch.db), intern = TRUE)
  z  
  z <- system(paste("cleanup.simple.server ", pid), intern = TRUE)
  z


if (checkPADIserver("ets"))
 {cat("getpadi  test ets ... ")
  z <- system(paste("getpadi ets B2001"), intern = TRUE)
  z  
 }
