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

  wait.for.server.to.terminate <- function(server)
    {# wait to ensure padi server is terminated
     for (i in 1:60)
       {if (!checkPADIserver(server)) break
        Sys.sleep(2)
       }
    }

  wait.for.server.to.start <- function(server)
    {# wait to ensure padi server is started
     for (i in 1:30)
       {if (checkPADIserver(server)) break
        Sys.sleep(1)
       }
    }
  

  server <- Sys.info()[["nodename"]]
  wait.for.server.to.terminate(server)
  if (checkPADIserver(server))
     stop("A server is already running. Testing stopped. Use cleanupPADIserver() or killPADIserver() to terminate it.")


  #  this would be better with freq=12, but simple server does not support it.
  ok <- putpadi(ts(1000* 1:50, start=c(1950,3), freq=1), series="tst1", 
            server=server, server.process="simple.server",
            cleanup.script="cleanup.simple.server",
            dbname=scratch.db, start.server=TRUE)

  if (ok) cat("ok\n") else  cat("failed! putpadi and starting server\n")

  wait.for.server.to.terminate(server)

  cat("getpadi  test  ... ")
  
  pid <- system(paste("simple.server ", scratch.db), intern = TRUE)
  pid
  z <- system(paste("x11arima ", server," tst1", scratch.db), intern = TRUE)
  z  
  z <- system(paste("cleanup.simple.server ", pid), intern = TRUE)
  z


if (checkPADIserver("ets"))
 {cat("x11arima  test ets ... ")
  z <- system(paste("x11arima ets B2001"), intern = TRUE)
  z  
 }
