###########################################################################

# Combination forecasting  functions.                       <<<<<<<<<<<<

###########################################################################


combineAndForecast <- function(model, new.data,  
                      overlapping.period.forecast.tag="g", forecast.tag="f") 

{# model should be a TSmodel.
 # new data should be a list with $data and $overriding.data.
 # It can also contain elements data.tag and overriding.data.tag, character string
 #   tags which are passed along to construct.data.to.override.horizon.
 # $overriding.data is used in place of data and model forecasts to the horizon
 # for which it is available. $overriding.data should also include any input (policy)
 # variables to the forecast horzon.
 # best.guess in the result is a combination of available data, overriding.data,
 # and predictions. 
 # first splice and fill with model predictions.
 con.data <- construct.data.to.override.horizon(new.data, model, plot=FALSE, 
                      forecast.tag=overlapping.period.forecast.tag) 
		      
 pred <-l(model, con.data ,predictT=dim(con.data$input)[1])$estimates$pred 
   # do residual analysis ?
# forecast<-forecast(l(model, con.data), percent=c(80,100,120), horizon=6, plot=F)
#   pchange<-percentChange(forecast[[1]],forecast[[2]],forecast[[3]], base=base,lag=12,cumulate=T,e=T)
 best.guess <-splice.tagged(con.data$output, pred, 
                  tag1=con.data$output.tags,tag2=forecast.tag) 
# the following result could also include con.data and pred, but it should be possible to
#    reconstruct these from the information in the list.

 invisible(list(model=model,
                data=new.data$data,
                overriding.data=new.data$overriding.data, 
                override=con.data$override,
                best.guess=best.guess))
}

reconstruct.combined.forecast <- function(combined.forecast) 
{# use the result of combineAndForecast to re-do and verify results
 con.data <- construct.data.to.override.horizon(combined.forecast, 
                   combined.forecast$model, plot=FALSE)
 pred <-l(combined.forecast$model, con.data ,predictT=dim(con.data$input)[1])$estimates$pred 
 best.guess <-splice.tagged(con.data$output, pred) 
 all(combined.forecast$best.guess==best.guess)
}

tfplot.combined.forecast <- function(x, 
       start=tfstart(x$data$output), end=tfend(x$data$output),
       select.inputs=NULL, select.outputs=NULL,
       Title="Projection", xlab=NULL, ylab=NULL, 
       graphs.per.page=5, mar=par()$mar, verbose=FALSE )
{ 
   if (verbose)
     {tfplot(x$data, xlab=xlab, ylab=ylab, graphs.per.page=graphs.per.page,
            start=start, end=end, 
	    mar=mar, Title="Data and combined forecast")
      tfplot(x$pred, xlab=xlab, ylab=ylab, graphs.per.page=graphs.per.page,
            start=start, end=end, 
	    mar=mar, Title="Model predictions (one step ahead for history)")
     }
   graph.data <- x$data
   graph.data$output <- x$best.guess
   if (is.null(select.inputs))  select.inputs  <- seq(dim(graph.data$input)[2])
   if (is.null(select.outputs)) select.outputs <- seq(dim(graph.data$output)[2])
   tfplot(graph.data, xlab=xlab, ylab=ylab, graphs.per.page=graphs.per.page,
           start=start, end=end, mar=mar, Title="Projection", 
           select.inputs=select.inputs, select.outputs=select.outputs)
#   tfplot(x$forecast[[2]],x$forecast[[1]],
#         x$forecast[[3]], start=start,
#         Title="Projection using future policy=most recent value and 20% higher and lower")
#   tfplot(x$pchange[[2]],x$pchange[[1]],
#         x$pchange[[3]],start=start, Title=
#    "Year over year percent change using future policy=most recent value and 20% higher and lower")
   invisible()
}

###########################################################################

# functions for misc. data retrieval, checking, and transformation <<<<<<<<<<<<

###########################################################################


construct.data.to.override.horizon <- function(new.data, model, 
         plot=TRUE, forecast.tag="f")
{# model should be a TSmodel.
 # new.data should be a list with $data and $overriding.data.
 # $overriding.data is used in place of $data and model forecasts to 
 # the horizon for which it is available. 
 #  Splice together $data and $overriding.data and if necessary
 #  calculate predictions for $overriding.data period and use them where $overriding.data
 #  or $data are not available, then return complete data set 
 #  to the end of the $overriding.data horizon, along with input data.
 #    (Typically the end of $overriding.data$output determines the periods
 #     for which forecast are combined and the end of $overriding.data$input
 #     determines how far into the future the model is used to extend the
 #     combined forecast. )
 #  Note that the $overriding.data is used in place of data in the 
 #  returned data set to allow for over-riding with anticipated data revisions.
 #  However, for any predictions during the combined.forecast period (ie. to augment
 #  $data and $overriding.data as returned by this function),  
 #  only $data is used and only to the last period for which observations
 #  for all variables are available.

 # if $overriding.data and $data overlap indicate override locations in 
 #     logical matrix dup:

 # tbind aligns the matrices
 dup <- tbind(outputData(new.data$data), outputData(new.data$overriding.data))
 if (!is.null(dup))
  {p <- nseriesOutput(new.data$data)
   dup <- (!is.na(dup[,1:p,drop=FALSE])) & (!is.na(dup[,(p+1):(2*p),drop=FALSE]))
  }

    # This can be used to provide a warning. eg
    #if (any(dup))
    #  {cat("WARNING:overriding is being used where data is available as follows:\n")
    #   print(dup)
    #  }

 z <- trimNA(outputData(new.data$data), endNAs=FALSE)
 zz <- new.data$overriding.data$output
 z <- splice(zz,z)
 start <- tfstart(z)
 if (is.null(new.data$data$input)) z.in <-NULL
 else
   {# note that $overriding.data does not override actual data for $input, but 
    #  that could be changed by reversing the order in the next line. (extra  
    #  code is necessary to give a warning.)
    z.in <-trimNA(splice(inputData(new.data$data),
                          inputData(new.data$overriding.data)))
    start <- latestStart(z, z.in)
    z.in <- tfwindow(z.in, start=start, warn=FALSE)
    if (any(is.na(z.in)))
       stop(paste("Input (exogenous) series data cannot be specified as NA. (note ",
                  "differenced data requires an overlap of one period at the end of ",
                  "historical data and the beginning of monitoring overriding data.)"))
   }
 z <- tfwindow(z, start=start, warn=FALSE)
 con.data <- TSdata(output=z,  input=z.in)

 # now augment $data and $overriding.data with model predictions for 
 #  the combined forecast period if necessary.
 if (any(is.na(outputData(con.data))))    
   {z <- TSdata(input = inputData(con.data),
                output= trimNA(outputData(new.data$data)))
    pred <- l(model,z, predictT= periodsOutput(con.data))$estimates$pred
    z <-splice.tagged(outputData(con.data),pred, 
                    tag1=con.data$output.tags, tag2=forecast.tag)
    outputData(con.data) <- z
   }

 con.data<- freeze(con.data)
 con.data$override <- dup
 if (plot &&  dev.cur() != 1 ) 
    {tfplot(con.data,start=(tfend(outputData(data))-c(1,0)), 
           Title="Historical and overriding data data")
    }
  invisible(con.data)
}

get.overriding.data <- function(file="overriding.data", 
 first.input="",first.output="", second.output="", m=1, p=10)
{#Get other data eg(monitoring or other forecast data) 
  #   N.B. This cues on series names in the file
  # m is the number of input series
  # p is the number of output series
  z  <- dsescan(file=file,what=character())
  first.in   <- (1:length(z))[z == first.input & !is.na(z)] 
  if (0== length(first.in))
     stop(paste("Cannot find keying string:", first.input," in file", file))
  first.out  <- (1:length(z))[z == first.output & !is.na(z)] 
  if (0== length(first.out))
     stop(paste("Cannot find keying string:", first.output," in file", file))
  second.out <- (1:length(z))[z == second.output & !is.na(z)] 
  if (0== length(second.out))
     stop(paste("Cannot find keying string:", second.output," in file", file))
  periodsInput <- (first.out-(first.in+m))/m     
  zz <- matrix(z[first.in:(first.out-1)],(periodsInput+1),m)
  input.names <- zz[1,]
  input <-  matrix( as.numeric(zz[2:(1+periodsInput),]), periodsInput,m)
  dimnames(input) <- list(NULL,input.names)
  input <- tframed(input, list(start=as.integer(z[1:2]),frequency=12))
  periodsOutput<- second.out-(first.out+1)
  zz <- matrix(z[first.out:length(z)],(periodsOutput+1),p)
  output.names <- zz[1,]
  output <-  matrix( as.numeric(zz[2:(1+periodsOutput),]), periodsOutput,p)
  dimnames(output) <- list(NULL,output.names)
  output <- tframed(output, list(start=as.integer(z[1:2]),frequency=12))
  TSdata(input=input , output=output)
}


restrict.overriding.data <- function(data, overriding.horizon=0)  
{#This function truncates overriding.data$output if it extends 
 #  overriding.horizon periods beyond the present. 
 year.mo <- c(dateParsed()$y,dateParsed()$m) - c(0,1) + c(0,overriding.horizon)
#check this - also note NAs should not be nec in overriding fame data
 data$output <-tfwindow(data$output, end=year.mo, warn=FALSE )
 invisible(data)
}

###########################################################################

# functions for e-mail of results of combination forecasting <<<<<<<<<<<<

###########################################################################

combinationMonitoring <- function(model, data.names,
   previous.data=NULL,
   overriding.data.names=NULL, 
   restrict.overriding.data=TRUE, overriding.horizon=0,
   mail.list=NULL,
   error.mail.list=NULL,
   message.title="Combination Monitoring",
   message.subject="Combination Monitoring",
   message.footnote=NULL,
   show.start= c(0,-3),
   show.end  = c(0,12),    
   report.variables=seriesNames(data.names),
   data.sub.heading=NULL,
   data.tag=" ",
   future.inputData.tag="p",
   overriding.data.tag="m",
   overlapping.period.forecast.tag="g",
   forecast.tag="f",
   run.again=FALSE,
   save.as=NULL)

{ # Step 0 - prepare message files and error checking
    error.message <- c(message.title, paste(dateParsed(), collapse="."),
              "An error condition occurred running combinationMonitoring.",
              "The message.file at the time of the error follows:") 
    message <- ""     
    on.exit(Sys.mail(error.mail.list, subject=paste("error ", message.subject),
                 body= c(error.message, message)))
    if ( dseclass(model)[1] == "TSestModel" ) model <- model$model
    if (!is.null(data.names$pad.end))
       {if(!data.names$pad.end)
          warning("pad.end in data definition may disable retrieving all data.")
       } 
    else if (!is.null(data.names$pad))
       {if(!data.names$pad)
          warning("pad in data definition may disable retrieving all data.")
       } 

# The following line can be removed if the code works reliably
   Sys.mail(error.mail.list,subject=paste("checking ",message.subject),
                           body=paste(dateParsed(), collapse="."))

 # Step 1 - retrieve & check for updated data  or
 #            initialize system and if previous.data is NULL
    if (is.null(previous.data))
      {data <- freeze(data.names)
       message <- "Initializing combination monitoring:"   
       status <- "Combination monitoring initialized."   
      }
    else if (run.again)
      {data <-previous.data$data  
       overriding.update <- previous.data$overriding.data
       status <- "Combination monitoring re-run."   
      }
    else
      {updated.data<-checkForValueChanges(data.names,
                           verification.data=previous.data$data,
                           discard.current=TRUE)
       if (is.null(overriding.data.names)) overriding.update<-list(FALSE)
       else overriding.update<-checkForValueChanges(overriding.data.names,
                           verification.data=previous.data$overriding.data)
       if(updated.data[[1]] | overriding.update[[1]])
         {status <- "Combination monitoring updated."     
          if(updated.data[[1]])
            {data <-updated.data$data
             message <- c("data updates: ", 
                 seriesNames(data)$input[updated.data$input],
                 seriesNames(data)$output[updated.data$output])
            }
          if(overriding.update[[1]])
            {overriding.data <- overriding.update$data
             if(restrict.overriding.data & (!is.null(overriding.data$output))) 
                overriding.data <- restrict.overriding.data(overriding.data, 
                                 overriding.horizon=overriding.horizon)
             message <- c(message,"monitoring data updates: ",
             seriesNames(overriding.data)$input[ overriding.update$input],
             seriesNames(overriding.data)$output[overriding.update$output])
            }
         }
       else
         {on.exit()
          return(invisible(list(data=previous.data, 
                status="Combination monitoring updates not necessary.")))
         }
      }

 # Step 2 - check data and overriding data
   # Sometimes data is available as soon as there are any days in a month (with
   #   ignore on in Fame). The following 4 lines trim these, but that may not be
   #   the best way to handle them.
   year.mo <- c(dateParsed()$y,dateParsed()$m) - c(0,1)
   data  <- tfwindow(data,  end=year.mo, warn=FALSE )
   fr <- c(tffrequency(data), 1)
      
   # check matching of starting date with end of available data.
   #   period for which all data is available in data
   end.ets <- tfend(trimNA(outputData(data))) 
   if (!is.null(overriding.data))
    {if (is.null(overriding.data$output))
     {overriding.data$output <- ts(matrix(NA, 1, nseriesOutput(data)),
                           end=tfend(data$output), 
                           frequency=tffrequency(data$output), 
                           names=dimnames(data$output)[[2]])
      if (!is.null(data$output.names))
         overriding.data$output.names <- data$output.names
     }
   else
     {if (!( (1+fr %*% end.ets) >= (fr %*%tfstart(overriding.data$output))))
        stop(paste("Monitoring data (or NAs) must be indicated after ", end.ets))
      if (1== latestEndIndex(outputData(data), outputData(overriding.data)))
         warning(paste("Overriding data file does not appear to be updated.",
         "True data is available past the end of the overriding data."))
    }}   

    if (is.null(overriding.data.names)) overriding.data <- NULL
    else
       overriding.data <- tagged(overriding.data,
          input.tags=future.inputData.tag, output.tags=overriding.data.tag)
    data <- tagged(data, tags=list(input=data.tag, output=data.tag))

 # Step 3 - run forecast
   # warnings from this should be mailed!!!!
    combined.forecast<-combineAndForecast(model, list(data, overriding.data),
           overlapping.period.forecast.tag=overlapping.period.forecast.tag, 
           forecast.tag=forecast.tag) 

 # Step 4 - write and mail files
    message <- c(message, "Projections are conditioned on forecast of ",
                            seriesNames(updated.data$data)$input, 
                          "                        with tranformation ",
                           data.names$input.transformations,
                          "The forecasts are now:")
    #starting and end period for plots & printing:
    start<-(tfend(combined.forecast$data$output)+show.start) 
    end  <-(tfend(combined.forecast$data$output)+show.end)
    # this can be too long if sufficient input data is not provided, so: 
    if ((fr %*% tfend(combined.forecast$best.guess)) < ((end-c(0,1)) %*% fr))
       end  <-tfend(combined.forecast$best.guess)

    report.variables$input<- 
            (report.variables$input == seriesNames(data.names)$input)
    report.variables$output<- 
            (report.variables$output == seriesNames(data.names)$output)


    rv <- tagged(
              combined.forecast$best.guess[,report.variables$output, drop=FALSE],
              tags= (attr(combined.forecast$best.guess,"tags")
                             ) [,report.variables$output, drop=FALSE])
    tframe(rv) <- tframe(combined.forecast$best.guess)
    inp <- splice(combined.forecast$data$input, 
                  combined.forecast$overriding.data$input,
                  tag1=data.tag, tag2=future.inputData.tag)
    rv <-tfwindow(cbind(inp,rv), start=start, end=end, warn=FALSE) 
    message <- c(message,fprint(rv, digits=5, sub.title=data.sub.heading)) 

    if (any(combined.forecast$override))
       {message <- c(message, "WARNING: overriding data is being used where historical data is available as follows:",
              combined.forecast$override)
       }

#    print(tfwindow(tbind(combined.forecast$data$input, combined.forecast$best.guess), 
#      start=start), digits=print.digits)

# The following needs a postscipt viewer like gv or pageview
#    postscript(file=graphics.file, width=7, height=8, pointsize=14,
#        horizontal=F, onefile=F, print.it=F, append=FALSE)
#    graph.combined.forecast(combined.forecast, start=start)
#    dev.off()
#    message <- c(message,"For graphic (in OpenWindows) type:\n    pageview ")
#    if("/" == substring(graphics.file,1,1) )
#             message <- c(message,graphics.file)
#    else
#      {pwd <- getwd()
#       if("/tmp_mnt" == substring(pwd,1,8)) pwd <-substring(pwd,9)
#       message <- c(message,paste(pwd,"/",graphics.file, sep=""))
#      }
#    message <- c(message," in a command tool window. (Be patient. It takes a few seconds.)")

    if (!is.null(message.footnote)) message <-c(message, message.footnote)
    Sys.mail(mail.list, subject=message.subject, body= message)


 # Step 4 - clean-up
    if (!is.null(save.as)) 
      {assign(save.as, combined.forecast, where=1)
#       file.copy( graphics.file, save.as)   # save graph
      } 
    if (updated.data[[1]] ) previous.data$data   <- updated.data$data
    if ( overriding.update[[1]])
       previous.data$overriding.data<- overriding.update$data
    on.exit()
    #return latest data for comparison next time
    invisible(list(data=previous.data, status=status, message=message)) 
}

