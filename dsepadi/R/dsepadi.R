
##############################################################################

############################################################################

#    functions for DSE interface to Time Series Protocol for      <<<<<<<<<<
#      Application Database Interface (TSPADI) data interface    <<<<<<<<<<

############################################################################

# Functions in this file now handle only the TSdata aspects. The database
#  interface has been pushed down into a time series matrix class "tfPADIdata"
#  defined in the file tfpadi.

# It is not clear that the c("TSPADIdata", "TSdata") class is necessary,
#  as the tfPADIdata should be just another time series matrix class,
#  but in an attempt to smooth the transition that is not being attempted 
#  in one step. Also, the name tfPADIdata might better be TSPADIdata,
#  but that would be confusing for the time being.

# Also, freeze is used to do both freeze and TSdata since these have been 
#  done together previously. This should be cleaned up.

############################################################

#   Definition of class c("TSPADIdata", "TSdata") <<<<<<<<<<

############################################################




# This TSPADIdata constructor was once called make.TSnames

TSPADIdata <- function( output=NULL,           input=NULL,
                        output.server=server,  input.server=server,
                        output.db=db,          input.db=db,
                        output.transforms="",  input.transforms="", 
                        output.names=NULL,     input.names=NULL,
                         start=NA, end=NA, frequency=NA, 
                         pad=FALSE, pad.start=pad, pad.end=pad,
                         server="", db="", start.server=NULL, 
                         server.process=NULL, cleanup.script=NULL,
                         stop.on.error=TRUE, warn=TRUE)
  {i <- if (is.null(input)) NULL else tfPADIdata(input, 
      transforms=input.transforms, names=input.names, 
      start=start, end=end, frequency=frequency,
      pad.start=pad.start, pad.end=pad.end, 
      server=input.server, db=input.db, start.server=start.server, 
      server.process=server.process,  cleanup.script=cleanup.script,
      stop.on.error=stop.on.error, warn=warn)

   o <- if (is.null(output)) NULL else tfPADIdata(output, 
      transforms=output.transforms, names=output.names, 
      start=start, end=end, frequency=frequency,
      pad.start=pad.start, pad.end=pad.end,
      server=output.server, db=output.db, start.server=start.server, 
      server.process=server.process, cleanup.script=cleanup.script,
      stop.on.error=stop.on.error, warn=warn)

    classed(list(input=i, output=o), c("TSPADIdata", "TSdata")) # constructor 
   }





TSPADIdata2 <- function(input=NULL, output=NULL,
    start = NA, end = NA, frequency = NA, pad.start = FALSE, 
    pad.end = FALSE,  start.server = NULL, 
    server.process = NULL, cleanup.script = NULL, stop.on.error = TRUE, 
    warn = TRUE)
  {i <- o <- NULL
   for (j in seq(length=length( input))) i <- cbind(i,  input[[j]])
   for (j in seq(length=length(output))) o <- cbind(o, output[[j]])
   TSPADIdata(input =     i[3,], output=            o[3,],
        input.server=     i[1,], output.server=     o[1,],
        input.db=         i[2,], output.db=         o[2,],
        input.transforms= i[4,], output.transforms= o[4,],
        input.names=      i[5,], output.names=      o[5,],
      start = start, end = end, frequency = frequency, pad.start = pad.start, 
      pad.end = pad.end, start.server = start.server, 
      server.process = server.process, cleanup.script = cleanup.script,
      stop.on.error = stop.on.error, warn = warn)
   }
  



modify.TSPADIdata <- function(obj,
                         start=NA, end=NA, frequency=NA, 
                         pad=NA, pad.start=pad, pad.end=pad,
                         server=NA, db=NA, start.server=NA, 
                         server.process=NA, cleanup.script=NA,
                         stop.on.error=NA, warn=NA,
			 append=NA, use.tframe=NA,
                        output=NA,             input=NA,
                        output.server=NA,  input.server=NA,
                        output.db=NA,          input.db=NA,
                        output.transforms=NA,  input.transforms=NA, 
                        output.names=NA,       input.names=NA,
			...)
  {#  (... further arguments, currently disregarded)
   if( (!all(is.na(c(input, input.server, input.db, input.transforms))))  |
        !all(is.na(c(start, end, frequency, pad.start, pad.end, server, db, 
           start.server, server.process, cleanup.script, stop.on.error, warn))))
    inputData(obj) <-  modify.tfPADIdata(inputData(obj),
      series=input, transforms=input.transforms, names=input.names, 
      start=start, end=end, frequency=frequency,
      pad.start=pad.start, pad.end=pad.end, 
      server=input.server, db=input.db, start.server=start.server, 
      server.process=server.process,  cleanup.script=cleanup.script,
      stop.on.error=stop.on.error, warn=warn, 
      append=append, use.tframe=use.tframe)

   if( (!all(is.na(c(output, output.server, output.db, output.transforms))))  |
        !all(is.na(c(start, end, frequency, pad.start, pad.end, server, db, 
           start.server, server.process, cleanup.script, stop.on.error, warn))))
    outputData(obj) <-  modify.tfPADIdata(outputData(obj),
      series=output, transforms=output.transforms, names=output.names, 
      start=start, end=end, frequency=frequency,
      pad.start=pad.start, pad.end=pad.end, 
      server=input.server, db=input.db, start.server=start.server, 
      server.process=server.process,  cleanup.script=cleanup.script,
      stop.on.error=stop.on.error, warn=warn,
      append=append, use.tframe=use.tframe)
    obj
   }



############################################################

#     methods for TSPADIdata class objects <<<<<<<<<<

############################################################

print.TSPADIdata <- function(x, ...) {print.default(x) }

is.TSPADIdata <- function(obj) {inherits(obj, "TSPADIdata") }

# TSdata methods should work for start, end, frequency


tsp.TSPADIdata <- function(x)
  {i <- tsp( inputData(x))
   o <- tsp(outputData(x))
   if (is.null(o)) return(i)
   if (is.null(i)) return(o)
   if (!all(i == o)) 
      warning("tsp results differ for input and output data. Using output")
   o
}

periodsInput.TSPADIdata <- function(x) periods( inputData(x))  
periodsOutput.TSPADIdata <- function(x) periods(outputData(x))  
periods.TSPADIdata <- function(x) periods(outputData(x))

 
inputData.TSPADIdata <- function(x, series=seq(length=nseriesInput(x)))
{if(is.null(x$input))  NULL else  x$input[ , series, drop=FALSE]}

outputData.TSPADIdata <- function(x,series=seq(length=nseriesOutput(x)))
{if(is.null(x$output)) NULL else  x$output[ , series, drop=FALSE]}


#  default should work
# nseriesInput.TSPADIdata <- function(x) {nseries( inputData(x))}
#nseriesOutput.TSPADIdata <- function(x) {nseries(outputData(x))}

# seriesNamesInput, seriesNamesOutput default should work

identifiers.TSPADIdata <- function(obj) 
	{list(input=identifiers(obj$input), output=identifiers(obj$output))}
sourcedb.TSPADIdata <- function(obj) 
	{list(input=sourcedb(obj$input), output=sourcedb(obj$output))}
sourceserver.TSPADIdata <- function(obj) 
	{list(input=sourceserver(obj$input), output=sourceserver(obj$output))}
sourceInfo.TSPADIdata <- function(obj) 
	{list(input=sourceInfo(obj$input), output=sourceInfo(obj$output))}


############################################################

#      Database interface for TSPADIdata  <<<<<<<<<<

############################################################



freeze.TSPADIdata <- function(data, timeout=60)
{ # This function retreives data from a PADI server using getpadi
  # See freeze.
  x <- data # temp fix so arg is called data as in generic
  if (is.null(x$input))
    {z <- TSdata(output=freeze(x$output))
     z$source <- x
     return(z)
    }
  if (is.null(x$output)) 
    {z <-TSdata(input=freeze(x$input))
     z$source <- x
     return(z)
    }
  # now so that input and output are aligned ...

   if (! testEqual(attr(x$input,  "start") , attr(x$output, "start")))
         warning("input and output start values do no match. Using outputs.")
   if (! testEqual(attr(x$input,  "end") , attr(x$output, "end")))
         warning("input and output end values do no match. Using outputs.")
   if (! testEqual(attr(x$input,  "frequency"), attr(x$output, "frequency")))
        warning("input and output frequency values do no match. Using outputs.")

   if(attr(x$input, "pad.start") != attr(x$output, "pad.start") |
      attr(x$input, "pad.end")   != attr(x$output, "pad.end")   )
      warning ("input and output padding attibutes do not match. Using outputs")
      
   if(!(testEqual(attr(x$input, "use.tframe"),    attr(x$output, "use.tframe"))  &
        testEqual(attr(x$input, "start.server"),  attr(x$output, "start.server"))&
        testEqual(attr(x$input, "server.process"),attr(x$output, "server.process"))&
        testEqual(attr(x$input, "cleanup.script"),attr(x$output,"cleanup.script"))&
        testEqual(attr(x$input, "stop.on.error"), attr(x$output,"stop.on.error"))&
        testEqual(attr(x$input, "warn"),          attr(x$output, "warn") )))
      warning ("input and output server attibutes do not match. Using outputs")

  r <- freeze(modify(x$output, # output first so attributes are used
         append=list(series=x$input[1,],server=x$input[2,],db=x$input[3,],
	             transforms=x$input[4,],names=seriesNames(x$input))))
  r <- TSdata(input=r, output=r)
  r$source <- x
  inputData(r)  <-  inputData(r, series=ncol(x$output)+seq(length=ncol(x$input)))
  outputData(r) <- outputData(r, series=seq(length=ncol(x$output)))
  r
}


availability.TSPADIdata <- function(obj, verbose=TRUE, timeout=60, ...)  
{#  (... further arguments, currently disregarded)
 # Indicate  dates for which data is available. 

 i <- if (0 ==  nseriesInput(obj)) NULL
      else availability( inputData(obj), verbose=verbose)
 o <- if (0 == nseriesOutput(obj)) NULL
      else availability(outputData(obj), verbose=verbose)
 if (is.null(i) & is.null(o)) stop("No data.")
 invisible(list(start = rbind(i$start, o$start),
                end   = rbind(i$end, o$end),
                frequency=c(i$frequency, o$frequency),
                series=c(i$series, o$series)))
}


putpadi.TSdata <- function(data, dbname, server=Sys.info()[["nodename"]], 
                   start.server=TRUE, server.process=padi.server.process(), 
                   cleanup.script=padi.cleanup.script(),
                   series=seriesNames(data),
                   user=Sys.info()[["user"]], passwd="",
                   stop.on.error=TRUE, warn=TRUE){   
   #dbname and server can be a single string in which case it is applied to
   # all series. Otherwise it should be a structure like series: a list with
   # elements input and output, each vectors with a string for each series.

   # This function uses tfputpadi and returns an TSPADIdata object which can 
   #  be used to fetch the data from the database. tfputpadi in turn 
   #  uses putpadi.default.

   m <-nseriesInput(data)
   p <-nseriesOutput(data)

   if(!is.list(dbname)) 
     {z <-dbname[1]
      dbname  <- list(input  = if (m==0) NULL else rep(z,m),
                      output = if (p==0) NULL else rep(z,p) )
     }

   if(!is.list(server)) 
     {z <-server[1]
      server <- list(input  = if (m==0) NULL else rep(z,m),
                     output = if (p==0) NULL else rep(z,p) )
     }

   if (m == 0) i <- NULL else
     {if(all (1 == start(inputData(data))))
         warning("Fame may choke on a start date of 1,1")
      mat <- tframed(inputData(data), list(start = start(inputData(data)), 
                frequency=frequency(inputData(data))))

      i <- tfputpadi(mat, server=server$input, dbname=dbname$input, 
         series=series$input,
         start.server = start.server, server.process = server.process, 
         cleanup.script = cleanup.script,
         user=user, passwd=passwd, stop.on.error=stop.on.error, warn=warn)   
     }
   if (p == 0) o <- NULL else
     {if(all (1 == start(outputData(data))))
         warning("Fame may choke on a start date of 1,1")
      mat <- tframed(outputData(data), list(start = start(outputData(data)),
                 frequency=frequency(data)))
      o <- tfputpadi(mat,  server=server$output, dbname=dbname$output, 
         series = series$output,
         start.server = start.server, server.process = server.process, 
         cleanup.script = cleanup.script,
         user=user, passwd=passwd, stop.on.error=stop.on.error, warn=warn)   

     }
   #This bypasses the constructor (structures are already built by tfputpadi):
   invisible(classed(list(input=i, output=o), c("TSPADIdata", "TSdata"))) # bypass constructor 
  }


#   The following function is supplied separately (with PADI ). The 
#   documentation is included here so it will integrate with DSE.



setTSPADIdata <- function()
 {# prompt for input and output series identifiers, sets class, etc.
  cat("This function prompts for the names and database locations for\n")
  cat("input and output series, until an empty line is entered.\n")
  cat("If your model has no input or no output then return an empty line.\n\n")
  cat("Input (exogenous) variables...\n")
  i <- settfPADIdata(preamble=FALSE)
  cat("Output (endogenous) variables...\n")
  o <- settfPADIdata(preamble=FALSE)
  data <- classed(list(input=i, output=o),   # bypass constructor (setTSPADIdata)
                  c("TSPADIdata", "TSdata"))  
  cat("The series may now be retrieved, in which case the data is\n")
  cat("  fixed as currently available, or they may be left `dynamic',\n")
  cat("  in which case they are retrieved using freeze.\n")
  cat("Retrieve data y/n:");key <- readline()
  if ((key =="y") | (key=="Y")) data <- freeze(data)
  data
}




retrieve.and.verify.data <- function(data.names,
             verification.data=verification.data, fuzz=1e-10)  
{# retrieve  data from a data base and do some verification. 
 # It is often useful if one of these is set:
 #    data.names$pad =T
 #    data.names$pad.end = T
 data <- freeze(data.names)
 #   check that data has not been rebased or otherwise messed up.
 if (0 != (nseriesInput(data)))
   {s <-startInput(verification.data)
    e <-endInput(verification.data)
    error <- inputData(verification.data) -
               tfwindow(inputData(data),start=s, end=e, warn=FALSE)
    if (fuzz < max(abs(error)) )
      {warning(paste("Retrieved input variables do not compare with the verification data.",
       "  Differences occur at ", sum(abs(error)>fuzz), " data points. ",
       " The maximum error is ", max(abs(error))))
       key<-as.character(parse(prompt="plot discrepancy?  y/n: "))
       if (key=="y" | key=="Y")
         {z <- TSdata(input=error)
          tfplot(z, 
              select.inputs=(1:nseriesInput(z))[apply(error,2,any)],
              select.outputs= 0)
         }
      key<-as.character(parse(prompt="plot data and verification data?  y/n: "))
       if (key=="y" | key=="Y")
         {graph.data <- data
          outputData(graph.data) <-tfwindow(outputData(data),start=s,end=e)
          inputData(graph.data)  <-tfwindow(inputData(data), start=s,end=e)
          tfplot(verification.data, graph.data,
            select.inputs=(1:nseriesInput(data))[apply(error,2,any)], select.outputs= 0)
         }
      }
   }
 s <-startOutput(verification.data)
 e <-endOutput(verification.data)
 error <-  outputData(verification.data) -
              tfwindow(outputData(data),start=s,end=e, warn=FALSE)
 if (fuzz < max(abs(error))  )
   {warning(paste("Retrieved output variables do not compare with the verification data.",
    "  Differences occur at ", sum(abs(error)>fuzz), " data points. ",
    " The maximum error is ", max(abs(error))))
    key<-as.character(parse(prompt="plot discrepancy?  y/n: "))
    if (key=="y" | key=="Y")
      {z <- TSdata(output=error)
       tfplot(z, select.inputs=0,
            select.outputs= (1:nseriesOutput(z))[apply(error,2,any)])
      }
    key<-as.character(parse(prompt="plot data and verification data?  y/n: "))
    if (key=="y" | key=="Y")
      {graph.data <- data
       graph.outputData(data) <-tfwindow(outputData(data),start=s,end=e)
       graph.inputData(data)  <-tfwindow(inputData(data), start=s,end=e)
       tfplot(verification.data, graph.data, select.inputs=0,
            select.outputs= (1:nseriesOutput(data))[apply(error,2,any)])
      }
   }
data
}


#######################################################################

#    TS PADI interface tests (from Brief User's Guide)   <<<<<<<<<<
#        moved to tests subdirectory
#######################################################################

