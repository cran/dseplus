
##############################################################################

#  The first section of this file contains generic definitions for objects 
#   which describe a database source for time series data. 

#  Following that are methods for the object tfPADIdata.

############################################################################

# Note: the constructors (e.g. tfPADIdata, TSPADIdata) cannot be generic.


# freeze and freeze.default moved to tframe (so everything else can be
#  included with dsepadi).


availability <- function(obj, ...) UseMethod("availability")


availability.default <- function(obj, names=NULL, server="ets", dbname="",
                verbose=TRUE, timeout=60, stop.on.error=TRUE, warn=TRUE, ...)  
{#  (... further arguments, currently disregarded)
 # Indicate  dates for which data is available. 
 # obj should be a character vector of data identifiers  
 if(!require("padi")) stop("This function requires the padi package.")
 if (1== length(server)) server  <- rep(server, length(obj))
 if (1== length(dbname)) dbname  <- rep(dbname, length(obj))

 # next 3 lines are to look after older style name forms at the BOC
 ets <- "ets" == substring(dbname,1,3)
 server[ets] <-"ets"
 dbname[ets] <- ""
 
 # server[server=="ets"] <- "padi"   # temporary kludge at the BOC

 s <- e <- f <- NULL
 for (i in 1:length(obj))
      {data <- getpadi(obj[i], dbname=dbname[i], server=server[i],
                stop.on.error=stop.on.error, use.tframe=TRUE, warn=warn, 
                pad=FALSE, timeout=timeout)
       s <- rbind(s, tfstart(data))
       e <- rbind(e, tfend(data))
       f <- c(f,tffrequency(data))
       if (verbose)
         {cat(obj[i]," from: ",tfstart(data))
          cat("  to: ",tfend(data))
          cat("   frequency ", tffrequency(data))
          if (!is.null(names)) cat("  ",names[i])
          cat("\n")
      }  }
  invisible(list(start=s, end=e, frequency=f, series=obj))
}



refresh <- function(data)
{src <- sourceInfo(data)
 if (is.null(src)) stop("data must include source information to use refresh.")
 freeze(src)
}


# extract series sourceInfo  (this is used by refresh so the result should
#   be the correct class, etc.
sourceInfo <- function(obj)UseMethod("sourceInfo")

sourceInfo.default <- function(obj){
    if (is.null(attr(obj, "source"))) stop("object does not have source information.")
    attr(obj, "source")
   }

# extract series identifiers
identifiers <- function(obj)UseMethod("identifiers")

identifiers.default <- function(obj){
    if (is.null(attr(obj, "source"))) stop("object does not have source information.")
    attr(obj, "source")[1,]
   } 
   

# extract series sourcedb
sourcedb <- function(obj)UseMethod("sourcedb")

sourcedb.default <- function(obj){
    if (is.null(attr(obj, "source"))) stop("object does not have source information.")
    attr(obj, "source")[3,]
   } 
   

# extract series sourceserver
sourceserver <- function(obj)UseMethod("sourceserver")

sourceserver.default <- function(obj){
    if (is.null(attr(obj, "source"))) stop("object does not have source information.")
    attr(obj, "source")[2,]
   }
   
############################################################################

#    functions define a time series matrix class "tfPADIdata"    <<<<<<<<<<
#      which uses the TSPADI  data interface. See also dsepadi   <<<<<<<<<<
#      file uses the methods here to define TSdata.              <<<<<<<<<<

############################################################################

# The PADI interface uses some calls to operating system specific functions:
#    -the function Sys.sleep is used in TSPADI.function.tests  
#    -the function local.host.netname defined in syskern.s
#    -previously the function user.name defined in the PADI interface 
#        software called a
#        program (getpwuid) in the $PADI/bin. This was previously done with
#        whoami() in syskern.s, which uses /usr/ucb/whoami (not system V unix).
#        It is important that Sys.info()[["user"]] return the same result as the C
#        function getpwuid in order for the padi interface to work properly.


############################################################

#   Definition of class c("tfPADIdata") <<<<<<<<<<

############################################################



tfPADIdata <- function(series,  server = "", db= "", transforms= "",  
           start=NA, end=NA, frequency=NA, names=NULL, 
           pad=FALSE, pad.start=pad, pad.end=pad,
           use.tframe=TRUE,
           start.server=FALSE, 
	   server.process=PADIserverProcess(), 
	   cleanup.script=PADIcleanupScript(),
           stop.on.error=TRUE, warn=TRUE)
  {# This is the constructor (but see settfPADIdata for a prompter).
   if (is.null(series)) return(NULL)
   if (is.null(names))   names <- series
   if(length(series) != length(names) )
           stop("number of names does not match number of series.")
   r <- rbind(series, server, db, transforms)
   dimnames(r) <- list(c("series", "server", "db", "transforms") ,names)
   attr(r, "start")     <- start
   attr(r, "end")       <- end
   attr(r, "frequency") <- frequency
   attr(r, "pad.start") <- pad.start
   attr(r, "pad.end")   <- pad.end
   attr(r,"use.tframe") <- use.tframe 

   attr(r, "start.server") <- start.server
   attr(r, "server.process") <- server.process
   attr(r, "cleanup.script") <- cleanup.script
   attr(r, "stop.on.error") <- stop.on.error
   attr(r, "warn")      <- warn
   class(r) <- "tfPADIdata"
   r
   }

settfPADIdata <- function(preamble=TRUE)
 {# prompt for series identifiers, set class, etc.
  if (preamble) 
    {cat("This function prompts for the names and database locations for\n")
     cat("series, until an empty line is entered.\n\n")
     cat("Variables...\n")
    }
  series <- server <- db <- transforms <- NULL
  repeat
    {key <- readline("  series..")
# cat(":",key,":") there seems to be a bug here. readline is not flushing
     if (""== key) break  else series <-c(series, key)
     server     <- c(server,     readline("  server.."))
     db         <- c(db,         readline("  database.."))
     transforms <-c(transforms,  readline("  transformation.."))
    } 
  if (is.null(series)) return(NULL) 
  cat("  starting year..");key <- readline()
     if (!(""== key)) 
       {start <- as.integer(key)
        cat("  starting period..");key <- readline()
        start <- c(start, as.integer(key))
        if(any(is.na(start)))
            cat("Warning: start improperly specified. NOT set!")
        }
     else start <- NA
  cat("  ending year..");key <- readline()
     if (!(""== key)) 
       {end. <- as.integer(key)
        cat("  ending period..");key <- readline()
        end. <- c(end., as.integer(key))
        if(any(is.na(end.))) cat("Warning: end improperly specified. NOT set!")
        }
     else end. <- NA

  data <- tfPADIdata(series, server=server, db=db, transforms=transforms,
                     start=start, end=end)
  if (preamble) 
    {cat("The series may now be retrieved, in which case the data is\n")
     cat("  fixed as currently available, or they may be left `dynamic',\n")
     cat("  in which case they are retrieved using freeze.\n")
     cat("Retrieve data y/n:");key <- readline()
     if ((key =="y") | (key=="Y")) data <- freeze(data)
    }
  data
}


modify <- function(obj, start=NA, end=NA, frequency=NA, 
                         pad=NA, pad.start=pad, pad.end=pad,
                         server=NA, db=NA, start.server=NA, 
                         server.process=NA, cleanup.script=NA,
                         stop.on.error=NA, warn=NA,
			 append=NA, use.tframe=NA, ...) UseMethod("modify") 

modify.tfPADIdata <- function(obj,
                         start=NA, end=NA, frequency=NA, 
                         pad=NA, pad.start=pad, pad.end=pad,
                         server=NA, db=NA, start.server=NA, 
                         server.process=NA, cleanup.script=NA,
                         stop.on.error=NA, warn=NA,
			 append=NA, use.tframe=NA,
                     series=NA, transforms=NA, names=NA, ...)
  {#  (... further arguments, currently disregarded)
   if (!is.na(series))     obj[1,] <- series
   if (!is.na(server))     obj[2,] <- server
   if (!is.na(db))         obj[3,] <- db
   if (!is.na(transforms)) obj[4,] <- transforms
   
   if (!is.na(names))
       dimnames(obj) <- list(c("series", "server", "db", "transforms") ,names)

  if (!all(is.na(append))) 
    {if (is.null(append$series))     append$db <- ""
     if (is.null(append$db))         append$db <- ""
     if (is.null(append$transforms)) append$transforms <- ""
     if (is.null(append$names))      append$names <- append$series
     if(length(append$series) != length(append$names) )
           stop("number of new names does not match number of new series.")
     newr <- rbind(append$series, append$server, append$db, append$transforms)
     newr <- cbind(obj,newr)
     dimnames(newr) <- list(c("series", "server", "db", "transforms") ,
                            c(dimnames(obj)[[2]], append$names))
     attr(newr, "start")     <- attr(obj, "start") 
     attr(newr, "end")       <- attr(obj, "end")
     attr(newr, "frequency") <- attr(obj, "frequency")
     attr(newr, "pad.start") <- attr(obj, "pad.start")
     attr(newr, "pad.end")   <- attr(obj, "pad.end")
     attr(newr,"use.tframe") <- attr(obj,"use.tframe")  

     attr(newr, "start.server")   <- attr(obj, "start.server")
     attr(newr, "server.process") <- attr(obj, "server.process")
     attr(newr, "cleanup.script") <- attr(obj, "cleanup.script")
     attr(newr, "stop.on.error")  <- attr(obj, "stop.on.error") 
     attr(newr, "warn")           <- attr(obj, "warn")  
     class(newr) <- "tfPADIdata"
     obj <- newr
    }
   
   if (!any(is.na(start)))          attr(obj, "start")     <- start
   if (!any(is.na(end)) )           attr(obj, "end")       <- end
   if (!    is.na(frequency))       attr(obj, "frequency") <- frequency

   if (!is.na(pad))            pad.start<- pad.end<- pad
   if (!is.na(pad.start))      attr(obj, "pad.start") <- pad.start
   if (!is.na(pad.end))        attr(obj, "pad.end")   <- pad.end
   if (!is.na(use.tframe))     attr(obj,"use.tframe") <- use.tframe 

   if (!is.na(start.server))   attr(obj, "start.server") <- start.server
   if (!is.na(server.process)) attr(obj, "server.process") <- server.process
   if (!is.na(cleanup.script)) attr(obj, "cleanup.script") <- cleanup.script
   if (!is.na(stop.on.error))  attr(obj, "stop.on.error") <- stop.on.error
   if (!is.na(warn))           attr(obj, "warn")      <- warn
   obj
   }




############################################################

#     methods for tfPADIdata class objects <<<<<<<<<<

# See also freeze.tfPADIdata and availability.tfPADIdata further below

############################################################


is.tfPADIdata <- function(obj) {inherits(obj, "tfPADIdata") }

print.tfPADIdata <- function(x, ...)
  {print.default(x)
   invisible(x)
  }

tframe.tfPADIdata <- function(x) 
   if(is.null(attr(x, "tframe"))) NA else attr(x, "tframe")

tfstart.tfPADIdata <- function(x, ...)
     {if(is.null(attr(x, "start"))) NA else attr(x, "start")}
tfend.tfPADIdata <- function(x, ...)
     {if(is.null(attr(x, "end")))   NA else attr(x, "end")}
tffrequency.tfPADIdata <- function(x, ...)
     {if(is.null(attr(x, "frequency")))   NA else attr(x, "frequency")}
tfperiods.tfPADIdata <- function(x) NA  # could be better
seriesNames.tfPADIdata <- function(x) {dimnames(x)[[2]]}
# nseries default should work


identifiers.tfPADIdata <- function(obj)  {obj[1,]}
sourceserver.tfPADIdata <- function(obj)  {obj[2,]}
sourcedb.tfPADIdata <- function(obj)  {obj[3,]}
sourceInfo.tfPADIdata <- function(obj)  {attr(obj,"source")} #used by refresh


"[.tfPADIdata" <- function(x, i, j, drop = FALSE) #N.B. FALSE
   {a <- attributes(x)
    y <- NextMethod("[")
    a$dim      <- dim(y)
    a$dimnames <- dimnames(y)
    attributes(y) <- a
    y
   }

tsp.tfPADIdata <- function(x)
  {start <-tfstart(x)
   end   <-  tfend(x)
   f <- tffrequency(x)
   if (length(start)==2) start <- start[1] + (start[2]-1)/f
   if (length(end)==2)   end   <- end[1]   + (end[2]-1)/f
   c(start, end, f)
  }

 

############################################################

#      Database interface for tfPADIdata  <<<<<<<<<<

############################################################



freeze.tfPADIdata <- function(data, timeout=60)
{ # This function retreives data from a PADI server using getpadi
  # A server specified as NULL or as "" is expanded to the localhost.
 if(!require("padi")) stop("This function requires the padi package.")

   # next 3 lines are to look after older style name forms at the BOC
   ets <- "ets" == substring(data["db",],1,3)
   data["server", ets] <- "ets"
   data["db",     ets] <- ""

   data["server", data["server",] ==""] <- Sys.info()[["nodename"]] 

   # missing attr is NULL but should be translated to getpadi defaults:
   IfNull <- function(a,b) {c(a,b)[1]}

   r  <- getpadi( data["series",], server=data["server",], dbname=data["db",],
     start.server=   IfNull(attr(data,"start.server"), TRUE),
     server.process= IfNull(attr(data,"server.process"), padi.server.process()),
     cleanup.script= IfNull(attr(data,"cleanup.script"), padi.cleanup.script()),
     starty= if(any(is.na(tfstart(data)))) 0 else tfstart(data)[1],
     startm= if(any(is.na(tfstart(data)))) 0 else tfstart(data)[2],
     endy=   if(any(is.na(tfend(data))))   0 else tfend(data)[1],
     endm=   if(any(is.na(tfend(data))))   0 else tfend(data)[2],
     transformations = data["transforms",],
     pad  = (attr(data,"pad.start") | attr(data,"pad.end") ),
     user =          IfNull(attr(data,"user"), Sys.info()[["user"]] ),
     passwd=         IfNull(attr(data,"passwd"),       ""  ),
     stop.on.error = IfNull(attr(data,"stop.on.error"), TRUE  ),
     use.tframe=     IfNull(attr(data,"use.tframe"),    FALSE  ), 
     warn=           IfNull(attr(data,"warn"),          TRUE  ),
     timeout= timeout)

 if (is.character(r)) stop(r)
 if (!attr(data,"pad.start")) r <- trimNA(r, startNAs=TRUE,  endNAs=FALSE)
 if (!attr(data,"pad.end") )  r <- trimNA(r, startNAs=FALSE, endNAs=TRUE)
 if (dim(r)[2] != dim(data)[2]) stop("Error retrieving data.")
 if ( !is.na(tffrequency(data)) && (tffrequency(data)) != tffrequency(r))
       warning("returned data frequency differs from request.")
 seriesNames(r) <- seriesNames(data)
 attr(r, "source") <- data 
 attr(r, "retrieval.date") <- dateParsed() 
 r
}

availability.tfPADIdata <- function(obj, verbose=TRUE, timeout=60, ...)  
{#  (... further arguments, currently disregarded)
# Indicate  dates for which data is available.
 # This requires retrieving series individually so they are not truncated.
 if(!require("padi")) stop("This function requires the padi package.")

   # next 3 lines are to look after older style name forms at the BOC
   ets <- "ets" == substring(obj["db",],1,3)
   obj["server", ets] <- "ets"
   obj["db",     ets] <- ""

   obj["server", obj["server",] ==""] <- Sys.info()[["nodename"]] 
   series <- obj["series", ]
   s <- e <- f <- NULL
   for (i in 1:length(series))
     {r <- getpadi( series[i], server=obj["server",i], dbname=obj["db",i], 
        start.server   = attr(obj,"start.server"), 
        server.process = attr(obj,"server.process"),
        cleanup.script = attr(obj,"cleanup.script"),
        starty=if(any(is.na(tfstart(obj)))) 0 else tfstart(obj)[1],
        startm=if(any(is.na(tfstart(obj)))) 0 else tfstart(obj)[2],
        endy=if(any(is.na(tfend(obj))))  0 else tfend(obj)[1],
        endm=if(any(is.na(tfend(obj))))  0 else tfend(obj)[2],
        transformations = obj["transforms",i],
        pad  = (attr(obj,"pad.start") | attr(obj,"pad.end")) ,
        user =if(is.null(attr(obj,"user"))) Sys.info()[["user"]] else attr(obj,"user"),
        passwd=if(is.null(attr(obj,"passwd")))  ""    else attr(obj,"passwd"),
        stop.on.error = attr(obj,"stop.on.error"),
        use.tframe=attr(obj,"use.tframe"), 
        warn=attr(obj,"warn"), timeout=timeout)

       s <- rbind(s, tfstart(r))
       e <- rbind(e, tfend(r))
       f <- c(f,tffrequency(r))
       if (verbose)
         {cat(series[i]," from: ",tfstart(r))
          cat("  to: ",tfend(r))
          cat("   frequency ", tffrequency(r))
          cat("  ",seriesNames(obj)[i])
          cat("\n")
      }  }
  invisible(list(start=s, end=e, frequency=f, series=series))
}



tfputpadi <- function(data,  
         server = Sys.info()[["nodename"]],
         dbname = "", 
         series = seriesNames(data),
         start.server = TRUE,
         server.process = padi.server.process(), 
         cleanup.script = padi.cleanup.script(),
         user = Sys.info()[["user"]], passwd= "",
         stop.on.error = TRUE, warn = TRUE)   
  {# This is just putpadi with a tfPADIdata object returned suitable for 
   #   retrieving the data.
 if(!require("padi")) stop("This function requires the padi package.")

   ok <- putpadi(data, server=server, dbname=dbname, series=series,
         start.server = start.server, server.process=server.process, 
         cleanup.script=cleanup.script,
         user=user, passwd=passwd,
         stop.on.error=stop.on.error, warn=warn ) 

   if (!all(ok)) stop("error putting data on database.")
  
   tfPADIdata( series, server=server, db=dbname, transforms="",  
           start=tfstart(data), end=tfend(data), frequency=tffrequency(data), 
           names=series, pad=FALSE, 
           use.tframe=TRUE, stop.on.error=stop.on.error, warn=warn)
  }


#   The following function is supplied separately (with PADI ). The 
#   documentation is included here so it will integrate with DSE.

#######################################################################

#     functions for converting defunct format FAMEdata structure
#         (these are primarily for use at the BOC)

#######################################################################

freeze.FAMEdata <- function(data)
  {stop("FAMEdata is defunct. Use FAMEdata.to.tfPADIdata to convert the structure")}


#######################################################################

