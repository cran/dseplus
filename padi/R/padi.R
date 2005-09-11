.onLoad  <- function(library, section) {
   # next require is necessary for bundle check to run examples,  
   # but does not seem to be necessary when packages are not bundled
   ok <- require("syskern", warn.conflicts=TRUE)
   ok <- ok & require("tframe", warn.conflicts=TRUE) 
   if(!ok) warning("This package requires the syskern and tframe packages.")
   invisible(TRUE)
   }

#.First.lib <- function(library, section)
#  {from <- paste(library, "/", section, sep="")
#   ok <- require("syskern", warn.conflicts=TRUE)
#   ok <- ok & require("tframe", warn.conflicts=TRUE) 
#   if(!ok) warning("This package requires the syskern and tframe packages.")
#   invisible(load.padi(from=from))
#  }


#load.padi <- function(from = Sys.getenv("PADI_LDLIB"))
#{# load C routines for PADI calls. 
# # NB. In S the default assumes the object code is in a subdirectory
# #	named $PADIabandonning/lib. In R this is ignored.
#
#   if("/" != substring(from, first=nchar(from))) from<-paste(from,"/", sep="")
#
#   if (is.R())  r <- library.dynam("padi")  # does not use from
#   else if(is.S())
#     {if ( Sys.info()$OSversion == "SunOS5"  & !is.Splus.pre3.3() )
#	   r <- dyn.load.shared("/usr/lib/libnsl.so")	  # splus 3.3 on SunOS5
#      ld <- paste(from,"splusclnt.o",sep="")
#      if (file.exists(ld)) r <- dyn.load(ld)
#      else {warning(paste(ld, " file to be loaded does not exist")); r <- FALSE}
#     }
#   invisible(r)
#}

##################################################################

PADIserverProcess <- function()
 {foo <- Sys.getenv("PADI_STARTUP")
  if(""== foo) foo <- "simple.server"
  paste(Sys.getenv("PADI"), "/", foo, sep = "")
 }
PADIcleanupScript <- function()
 {foo <- Sys.getenv("PADI_CLEANUP")
  if(""== foo) foo <- "cleanup.simple.server"
  paste(Sys.getenv("PADI"), "/", foo, sep = "")
 }

# the function user.name is removed. Use  Sys.info()[["user"]]

#   There is also a small dependance on tframe.s (seriesNames in putpadi)
#   and a more important dependency in getpadi when use.tframe=T.


##################################################################

startPADIserver <-function(server=Sys.info()[["nodename"]], 
			      server.process=PADIserverProcess(), dbname=NULL)
{# This function executes the server.process on server to start a PADI server.
 # The script or executable server.process should be on the Unix path. It
 #    would typically be found in $PADI/
 # The default is a "local mode" Fame server.
 # The function returns the machine and pid which can be used to
 #   kill (terminate) the process.
 # Note that starting the server takes a few seconds. If the interface is 
 #   being used a lot, it may be faster to start the server once and keep
 #   it running. Then specify different databases as necessary to getpadi and 
 #    putpadi. 
 # If dbname is not specified then a tempfile will be used as the initial
 #  database.
 # If dbname is specified, then the file will be created (by the 
 #   fame.server script) if it does not exist. Also, the fame.server
 #   script will append the fame convention ".db" so that should not be 
 #   specified.
 if (is.null(dbname)) dbname <- tempfile()
 # pid <-unix(paste("rsh",server,server.process,dbname,"&")) 
 # rsh above seems to cause problems terminating.
 # unix(paste("which ", server.process))
 pid <-system(paste(server.process, dbname), intern=TRUE)
 if(!is.na(charmatch("Error", pid)))
    stop(paste(server.process," did not initialize correctly.\n",
       " If the process started it may not have been terminated.\n",
       " Error messages follow:\n",
                paste(pid, collapse="\n")))
 list(server=server, pid=pid)
}

cleanupPADIserver <-function(process, cleanup.script=PADIcleanupScript())
{#This function terminates a server process specified by process, which should
 # be  a list as return by  startPADIserver. (process is just the argument
 # to cleanup.script and not actually the process ID - for the script  
 # cleanup.fame.server it is actually the pid of the process which 
 # started the server.)
#  r <- unix(paste("rsh", process$server,cleanup.script, process$pid))
  r <- system(paste(cleanup.script, process$pid), intern=TRUE)
  invisible(r)
}

killPADIserver <-function(kill.script="killserver", server.process=PADIserverProcess())
{#This function looks for a server with the name server.process and terminates 
 #  it. The function cleanupPADIserver is preferable as it cleans up
 #  scratch files, etc. (but it requires additional info.)
  r <- system(paste(kill.script, server.process), intern=TRUE)
  invisible(r)
}

checkPADIserver <-function(server=Sys.info()[["nodename"]],
        user=Sys.info()[["user"]], timeout=50)
{# This function attempts to get a non existent object and examines the error 
 #    msg to determine if a PADI server is running. It returns T or F.
 # Note that the error message string may not be standard on all operating
 #   systems, so the result may be ambiguous.

 msg <- .C("getpadi",user, "passwd", as.character(server), "","",
#	       working.freq=as.integer(0),  v2003.10 needs this, but...
	       as.integer(1),as.integer(1),as.integer(1),
               as.integer(0),as.integer(0),as.integer(0),as.integer(0),
               as.integer(10),double(1), as.integer(1),
               msg=paste(rep(" ",120), collapse=""), # error message buffer
               buffsz=as.integer(120),      # size of error message buffer
               as.integer(timeout),  # wait before generating error
	       PACKAGE="padi")$msg

 if (!is.na(charmatch("missing object name", msg))) return(TRUE)
 cat(msg)
 if (!is.na(charmatch(": RPC: Remote system error - Connection refused\n", msg)))
    return(FALSE)
 if (!is.na(charmatch(": RPC: Remote system error - Connection timed out\n", msg)))
    return(FALSE)
 if (!is.na(charmatch(": RPC: Program not registered", msg)))
    return (FALSE)
 if (!is.na(charmatch(": RPC: Unknown host", msg)))
    return (FALSE)
# else if(!is.na(charmatch("RPC: Server can't decode arguments",
#      msg)))return(TRUE)
# else  warning(paste(
#   "Ambiguous result from checkPADIserver. Assuming server is running.", 
#    msg,sep="\n"))
# TRUE
  FALSE
}



getpadi.default <-function(series,server=Sys.info()[["nodename"]], dbname="",
        start.server=TRUE, server.process=PADIserverProcess(),
        cleanup.script=PADIcleanupScript(),
        starty=0,startm=0,startd=1, endy=0,endm=0,endd=1, 
        nobs=0,max.obs=2000, transformations=NULL, pad=FALSE,
        user=Sys.info()[["user"]], passwd="",
        stop.on.error=TRUE, use.tframe=FALSE, warn=TRUE, timeout=60)
{   if (is.null(user)) stop("user cannot be NULL.")
    if (is.null(passwd)) stop("passwd cannot be NULL.")
    for (i in 1:length(server))
       if (is.null(server[i])) stop("server elements cannot be NULL.")
    for (i in 1:length(series))
       if (is.null(series[i])) stop("series elements cannot be NULL.")
    for (i in 1:length(dbname))
       if (is.null(dbname[i])) stop("dbname elements cannot be NULL.")

    if (all(server == server[1])) N <- 1
    else N <- length(server)

    for (i in 1:N)
      {server.running <- checkPADIserver(server[i])
       #  The above can falsely indicate T when a server is just terminating.
       if (!server.running) 
         {if (start.server && (server[i]==Sys.info()[["nodename"]]))
            {pid<-startPADIserver(server=server[i], dbname[i],
                                    server.process=server.process)
             on.exit(cleanupPADIserver(pid,               
                     cleanup.script=cleanup.script))
            }
          else stop(paste("Server",server[i], "is not running."))
         }
       }
     N <- length(series)
     if(length(server)==1) server <-rep(server,N)
     if(length(server)!=N) 
           stop("The length of server must correspond to the length of series")
     if(length(dbname)==1) dbname <-rep(dbname,N)
     if(length(dbname)!=N) 
           stop("The length of dbname must correspond to the length of series")
     if(!is.null(transformations))
        {if(length(transformations)==1) transformations <-rep(transformations,N)
         if(length(transformations)!=N) 
           stop("The length of transformations must correspond to the length of series.")
        } 

 
     ok <- NULL
     mat <-NULL
     buffsz <-120
     arsize <- max(nobs,max.obs)  # max. number of data points. Crude but...
     for (i in 1:length(series))
       {ns<-.C("getpadi",
               as.character(user),
               as.character(passwd),
               as.character(server[i]),
               as.character(series[i]),
               as.character(dbname[i]),
#	       working.freq=as.integer(0),  v2003.10 needs this, but...
               starty=as.integer(starty),
               startm=as.integer(startm),
               startd=as.integer(startd),
               endy=as.integer(endy),
               endm=as.integer(endm),
               endd=as.integer(endd),
               freq=as.integer(0),
               nobs=as.integer(nobs), # no. of observations
               data=double(arsize),   # returned data in double precision
               as.integer(arsize),    # max. no. of observations (size of data)
               msg=paste(rep(" ",buffsz), collapse=""), # error message buffer
               buffsz=as.integer(buffsz),    # size of error message buffer
               as.integer(timeout),  # wait before generating error
	       PACKAGE="padi")

        if (ns$buffsz > 0) 
          {if (stop.on.error) stop(paste("Error getting series ",series[i],
                                " on server ",server[i],
                                " data base ", dbname[i],"!\n", ns$msg))
           else return(ns$msg)
          }
        else 
          {data<- ns$data[1:ns$nobs]
           # replace NA indicator. 
           #  Note: the PADI protocol can control indicator values
           data[data ==
            170141507979022890158414086871904681984.000000]<-NA #Not Available
           data[data ==
            170141670238299719371777478449914970112.000000]<-NA #No Data
           data[data ==
            170141345719746060945050695293894393856.000000]<-Inf #Not Computable
           freq <- ns$freq
           if (freq == 52)
             {if (warn) warning("The period conversion for weekly data may not work properly in all cases.")
              m <- c(31,28,31,30,31,30,31,31,30,31,30,31)
              period<- 1+((sum(m[seq(12) < ns$startm])+ns$startd-1)*52) %/% 365
             }
           else period <- 1+ (((ns$startm-1)*freq) %/%12)
           if (use.tframe)
              z<-tframed(data, list(start=c(ns$starty,period),frequency=freq))
           else
              z<-ts(data, start=c(ns$starty,period),frequency=freq)
	   if(!is.null(transformations))
             {if(0!=nchar(transformations[i]))
                {if (!use.tframe)
                   warning("transformations may not work with data which is not tfamed. Set use.tframe=T in getpadi.")
                 z<-eval(call(transformations[i], z))
                }
             }
           ok <- c(ok,i)
           mat <- append(mat,list(z))
          }
       }

   data <-mat[[1]]
   if(1==length(mat)) data <-tbind(data,NULL) # for vector case
   for (i in ok[-1])  data <-tbind(data,mat[[i]], pad.start=pad, pad.end=pad)
   seriesNames(data) <- series[ok]
   data
}


putpadi.default <-function(data,  server=Sys.info()[["nodename"]], dbname="",
        series=seriesNames(data),
        start.server=TRUE, server.process=PADIserverProcess(),
        cleanup.script=PADIcleanupScript(),
        user=Sys.info()[["user"]], passwd="",
        stop.on.error=TRUE, warn=TRUE, timeout=60)
{   N <- nseries(data)
    if (is.null(series)) stop("series names must be supplied.")
    if (N != length(series) )
       stop("The length of series must correspond to the number of series.")
    if (1 == length(server) ) server <- rep(server, N)
    if (N != length(server) )
       stop("The length of server must correspond to the number of series.")
    if (is.null(dbname)) stop("dbname must be supplied.")
    if (1 == length(dbname) ) dbname <- rep(dbname, N)
    if (N != length(dbname) )
       stop("The length of dbname must correspond to the number of series.")

#    if (server.process != substring(cleanup.script,9))
#      warning(paste("cleanup.script(",cleanup.script,
#                    ") is not named cleanup.",server.process, sep="") )

   # next 3 lines are to look after older style name forms at the BOC
    ets <- "ets" == substring(dbname,1,3)
    server[ets] <-"ets"
    dbname[ets] <- ""

    if (all(server == server[1])) Ns <- 1
    else Ns <- length(server)
    for (i in 1:Ns)
      {server.running <- checkPADIserver(server[i])
       if (!server.running) 
         {if (start.server && (server[i]==Sys.info()[["nodename"]]))
            {pid<-startPADIserver(server=server[i], dbname[i],
                                    server.process=server.process)
             on.exit(cleanupPADIserver(pid,
                     cleanup.script=cleanup.script))
            }
          else stop(paste("Server",server[i], "is not running."))
         }
       }

    nobs <-periods(data)
    starty<- tfstart(data)[1]
    freq  <- tffrequency(data)
    startd <- 1
    if (freq == 1) startm<- 1
    else if (freq == 4)  startm<- tfstart(data)[2] * 3
    else if (freq == 12) startm<- tfstart(data)[2]
    else if (freq == 52) 
      {startm<- 1 + ( ((tfstart(data)[2]-1) *12) %/% 52)
       m <- sum(c(31,28,31,30,31,30,31,31,30,31,30,31)[seq(12) < startm])
       startd <- 2+  7*(tfstart(data)[2]-1) - m  # near beginning of week
       if (warn)
          warning("The period conversion for weekly data may not work properly in all cases.")
      }
    else 
      {freq <- 1
       starty <- 1
       startm <- 1
       startd <- 1
      }
    if((server.process=="fame.server") && (1==starty))
      {warning("Fame does not like starting year = 1 for time series.")
       #freq<-??? In Fame starty=1 is may be valid for vectors
      }
    buffsz <- 120  # warning message buffer
    if (!is.matrix(data)) data <- matrix(data, periods(data),nseries(data))
    ok <- rep(NA,N)
    for (i in 1:N)
       {ns<-.C("putpadi",
               as.character(user),
               as.character(passwd),
               as.character(server[i]),
               as.character(series[i]),
               as.character(dbname[i]),
               freq=as.integer(freq),
               starty=as.integer(starty),
               startm=as.integer(startm),
               startd=as.integer(startd),
               nobs=as.integer(nobs),
               data=as.double(data[,i]), # data in double precision
               msg=paste(rep(" ",buffsz), collapse=""),
               buffsz=as.integer(buffsz),
               as.integer(timeout),  # wait before generating error
	       PACKAGE="padi")
         if (ns$buffsz > 0) 
           {if (stop.on.error)
               stop(paste("Error writing series: ",series[i],
                     " on server ",server[i]," data base ", dbname[i],
                     "\n",ns$msg, sep=""))
            ok[i] <- FALSE
           }
         else ok[i] <- TRUE
        }
   invisible(ok)
}    


