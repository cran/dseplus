
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
                         stop.on.error=T, warn=T)
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
    server.process = NULL, cleanup.script = NULL, stop.on.error = T, 
    warn = T)
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
                        output=NA,             input=NA,
                        output.server=NA,      input.server=NA,
                        output.db=NA,          input.db=NA,
                        output.transforms=NA,  input.transforms=NA, 
                        output.names=NA,       input.names=NA,
                         start=NA, end=NA, frequency=NA, 
                         pad.start=NA, pad.end=NA,
                         server=NA, db=NA, start.server=NA, 
                         server.process=NA, cleanup.script=NA,
                         stop.on.error=NA, warn=NA)
  {
   if( (!all(is.na(c(input, input.server, input.db, input.transforms))))  |
        !all(is.na(c(start, end, frequency, pad.start, pad.end, server, db, 
           start.server, server.process, cleanup.script, stop.on.error, warn))))
    input.data(obj) <-  modify.tfPADIdata(input.data(obj),
      series=input, transforms=input.transforms, names=input.names, 
      start=start, end=end, frequency=frequency,
      pad.start=pad.start, pad.end=pad.end, 
      server=input.server, db=input.db, start.server=start.server, 
      server.process=server.process,  cleanup.script=cleanup.script,
      stop.on.error=stop.on.error, warn=warn)

   if( (!all(is.na(c(output, output.server, output.db, output.transforms))))  |
        !all(is.na(c(start, end, frequency, pad.start, pad.end, server, db, 
           start.server, server.process, cleanup.script, stop.on.error, warn))))
    output.data(obj) <-  modify.tfPADIdata(output.data(obj),
      series=output, transforms=output.transforms, names=output.names, 
      start=start, end=end, frequency=frequency,
      pad.start=pad.start, pad.end=pad.end, 
      server=input.server, db=input.db, start.server=start.server, 
      server.process=server.process,  cleanup.script=cleanup.script,
      stop.on.error=stop.on.error, warn=warn)
    obj
   }



############################################################

#     methods for TSPADIdata class objects <<<<<<<<<<

############################################################

print.TSPADIdata <- function(x, ...) {print.default(x) }

is.TSPADIdata <- function(obj) {inherits(obj, "TSPADIdata") }

# TSdata methods should work for start, end, frequency


tsp.TSPADIdata <- function(x)
  {i <- tsp( input.data(x))
   o <- tsp(output.data(x))
   if (is.null(o)) return(i)
   if (is.null(i)) return(o)
   if (!all(i == o)) 
      warning("tsp results differ for input and output data. Using output")
   o
}

input.periods.TSPADIdata <- function(x) periods( input.data(x))  
output.periods.TSPADIdata <- function(x) periods(output.data(x))  
periods.TSPADIdata <- function(x) periods(output.data(x))

 
input.data.TSPADIdata <- function(x, series=seq(length=input.dimension(x)))
{if(is.null(x$input))  NULL else  x$input[ , series, drop=FALSE]}

output.data.TSPADIdata <- function(x,series=seq(length=output.dimension(x)))
{if(is.null(x$output)) NULL else  x$output[ , series, drop=FALSE]}


#  default should work
# input.dimension.TSPADIdata <- function(x) {nseries( input.data(x))}
#output.dimension.TSPADIdata <- function(x) {nseries(output.data(x))}

# input.series.names, output.series.names default should work

identifiers.TSPADIdata <- function(obj) 
	{list(input=identifiers(obj$input), output=identifiers(obj$output))}
sourcedb.TSPADIdata <- function(obj) 
	{list(input=sourcedb(obj$input), output=sourcedb(obj$output))}
sourceserver.TSPADIdata <- function(obj) 
	{list(input=sourceserver(obj$input), output=sourceserver(obj$output))}
source.info.TSPADIdata <- function(obj) 
	{list(input=source.info(obj$input), output=source.info(obj$output))}


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

   if (! test.equal(attr(x$input,  "start") , attr(x$output, "start")))
         warning("input and output start values do no match. Using outputs.")
   if (! test.equal(attr(x$input,  "end") , attr(x$output, "end")))
         warning("input and output end values do no match. Using outputs.")
   if (! test.equal(attr(x$input,  "frequency"), attr(x$output, "frequency")))
        warning("input and output frequency values do no match. Using outputs.")

   if(attr(x$input, "pad.start") != attr(x$output, "pad.start") |
      attr(x$input, "pad.end")   != attr(x$output, "pad.end")   )
      warning ("input and output padding attibutes do not match. Using outputs")
      
   if(!(test.equal(attr(x$input, "use.tframe"),    attr(x$output, "use.tframe"))  &
        test.equal(attr(x$input, "start.server"),  attr(x$output, "start.server"))&
        test.equal(attr(x$input, "server.process"),attr(x$output, "server.process"))&
        test.equal(attr(x$input, "cleanup.script"),attr(x$output,"cleanup.script"))&
        test.equal(attr(x$input, "stop.on.error"), attr(x$output,"stop.on.error"))&
        test.equal(attr(x$input, "warn"),          attr(x$output, "warn") )))
      warning ("input and output server attibutes do not match. Using outputs")

  r <- freeze(modify(x$output, # output first so attributes are used
         append=list(series=x$input[1,],server=x$input[2,],db=x$input[3,],
	             transforms=x$input[4,],names=series.names(x$input))))
  r <- TSdata(input=r, output=r)
  r$source <- x
  input.data(r)  <-  input.data(r, series=ncol(x$output)+seq(length=ncol(x$input)))
  output.data(r) <- output.data(r, series=seq(length=ncol(x$output)))
  r
}


availability.TSPADIdata <- function(obj, verbose=T, timeout=60)  
{# Indicate  dates for which data is available. 

 i <- if (0 ==  input.dimension(obj)) NULL
      else availability( input.data(obj), verbose=verbose)
 o <- if (0 == output.dimension(obj)) NULL
      else availability(output.data(obj), verbose=verbose)
 if (is.null(i) & is.null(o)) stop("No data.")
 invisible(list(start = rbind(i$start, o$start),
                end   = rbind(i$end, o$end),
                frequency=c(i$frequency, o$frequency),
                series=c(i$series, o$series)))
}


putpadi.TSdata <- function(data, dbname, server=Sys.info()[["nodename"]], 
                   start.server=T, server.process=padi.server.process(), 
                   cleanup.script=padi.cleanup.script(),
                   series=series.names(data),
                   user=Sys.info()[["user"]], passwd="",
                   stop.on.error=T, warn=T){   
   #dbname and server can be a single string in which case it is applied to
   # all series. Otherwise it should be a structure like series: a list with
   # elements input and output, each vectors with a string for each series.

   # This function uses tfputpadi and returns an TSPADIdata object which can 
   #  be used to fetch the data from the database. tfputpadi in turn 
   #  uses putpadi.default.

   m <-input.dimension(data)
   p <-output.dimension(data)

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
     {if(all (1 == start(input.data(data))))
         warning("Fame may choke on a start date of 1,1")
      mat <- tframed(input.data(data), list(start = start(input.data(data)), 
                frequency=frequency(input.data(data))))

      i <- tfputpadi(mat, server=server$input, dbname=dbname$input, 
         series=series$input,
         start.server = start.server, server.process = server.process, 
         cleanup.script = cleanup.script,
         user=user, passwd=passwd, stop.on.error=stop.on.error, warn=warn)   
     }
   if (p == 0) o <- NULL else
     {if(all (1 == start(output.data(data))))
         warning("Fame may choke on a start date of 1,1")
      mat <- tframed(output.data(data), list(start = start(output.data(data)),
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



set.TSPADIdata <- function()
 {# prompt for input and output series identifiers, sets class, etc.
  cat("This function prompts for the names and database locations for\n")
  cat("input and output series, until an empty line is entered.\n")
  cat("If your model has no input or no output then return an empty line.\n\n")
  cat("Input (exogenous) variables...\n")
  i <- set.tfPADIdata(preamble=F)
  cat("Output (endogenous) variables...\n")
  o <- set.tfPADIdata(preamble=F)
  data <- classed(list(input=i, output=o),   # bypass constructor (set.TSPADIdata)
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
 if (0 != (input.dimension(data)))
   {s <-input.start(verification.data)
    e <-input.end(verification.data)
    error <- input.data(verification.data) -
               tfwindow(input.data(data),start=s, end=e, warn=F)
    if (fuzz < max(abs(error)) )
      {warning(paste("Retrieved input variables do not compare with the verification data.",
       "  Differences occur at ", sum(abs(error)>fuzz), " data points. ",
       " The maximum error is ", max(abs(error))))
       key<-as.character(parse(prompt="plot discrepancy?  y/n: "))
       if (key=="y" | key=="Y")
         {z <- TSdata(input=error)
          tfplot(z, 
              select.inputs=(1:input.dimension(z))[apply(error,2,any)],
              select.outputs= 0)
         }
      key<-as.character(parse(prompt="plot data and verification data?  y/n: "))
       if (key=="y" | key=="Y")
         {graph.data <- data
          output.data(graph.data) <-tfwindow(output.data(data),start=s,end=e)
          input.data(graph.data)  <-tfwindow(input.data(data), start=s,end=e)
          tfplot(verification.data, graph.data,
            select.inputs=(1:input.dimension(data))[apply(error,2,any)], select.outputs= 0)
         }
      }
   }
 s <-output.start(verification.data)
 e <-output.end(verification.data)
 error <-  output.data(verification.data) -
              tfwindow(output.data(data),start=s,end=e, warn=F)
 if (fuzz < max(abs(error))  )
   {warning(paste("Retrieved output variables do not compare with the verification data.",
    "  Differences occur at ", sum(abs(error)>fuzz), " data points. ",
    " The maximum error is ", max(abs(error))))
    key<-as.character(parse(prompt="plot discrepancy?  y/n: "))
    if (key=="y" | key=="Y")
      {z <- TSdata(output=error)
       tfplot(z, select.inputs=0,
            select.outputs= (1:output.dimension(z))[apply(error,2,any)])
      }
    key<-as.character(parse(prompt="plot data and verification data?  y/n: "))
    if (key=="y" | key=="Y")
      {graph.data <- data
       graph.output.data(data) <-tfwindow(output.data(data),start=s,end=e)
       graph.input.data(data)  <-tfwindow(input.data(data), start=s,end=e)
       tfplot(verification.data, graph.data, select.inputs=0,
            select.outputs= (1:output.dimension(data))[apply(error,2,any)])
      }
   }
data
}


#######################################################################

#    TS PADI interface tests (from Brief User's Guide)   <<<<<<<<<<
#        moved to tests subdirectory
#######################################################################

