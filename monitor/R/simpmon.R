###########################################################################

# Simple monitoring functions and data checking        <<<<<<<<<<<<

###########################################################################


checkForValueChanges <- function(data.names, verification.data,
     discard.current=FALSE, ignore.before= NULL, fuzz=1e-10)
  { # Check if data is modified or if more data is available.
    # data.names is an object of class c("TSPADIdata","TSdata").
    # verification.data is an object of class TSdata.
    # T is returned for any series which has a modified value.
    #   NA in one series and not in the other is counted as modified.
    # If data are not the same length then series are padded with NA
    #  so new NAs on the end will not count as a change.
    # It is assumed that changes happen at the end (not the beginning) of
    #   the data. The data is trimmed at the beginning to the point where
    #   all series are available. (this simplifies padding the end with NA)
    # If ignore.before is not NULL it should indicate a year and period
    #   before which data is trimmed, no comparison is performed and the
    #   data before is not returned. If there are NAs at the beginning then
    #   trimming as described above may make the data even shorter than
    #   indicated by ignore.before.
    # discard.current controls whether current period data is considered.
    #  (some series are available for a period from the first day of the
    #   period, and are updated daily. Usually these should be discarded
    #   by setting discard.current=T)

   data <- freeze(data.names) 
   if (discard.current)
     {year.mo <- c(dateParsed()$y,dateParsed()$m) - c(0,1)
      data  <- tfwindow( data,  end=year.mo, warn=FALSE )
     }
   if (!is.null(ignore.before)) 
     {data <- tfwindow(data, start= ignore.before)
      verification.data <-tfwindow(verification.data, start= ignore.before)
     }
   data <-trimNA(data, start.=TRUE, end.=FALSE)
   verification.data <-trimNA(verification.data, start.=TRUE, end.=FALSE)
   # which series are changed:
   if (is.null(seriesNamesInput(data.names))) in.up <- NULL
   else
     {ld <-periodsInput(data)
      lv <-periodsInput(verification.data)
      l <- max(ld, lv)
      if (ld < l)
        inputData(data) <- ts(rbind(inputData(data),  
                                     matrix(NA,l-ld,nseriesInput(data))),
                     start=start(inputData(data)),  frequency=frequency(data))
      if (lv < l)
        inputData(verification.data) <- ts(rbind(inputData(verification.data),
                                        matrix(NA,l-lv, nseriesInput(data))),
                     start=start(inputData(verification.data)),
                     frequency=frequency(verification.data))
      z <- (is.na(inputData(data)) & is.na(inputData(verification.data)))   # both NA
    # next fixes an Splus bug (IMHO) that the dim is dropped for col matrix
      if (!is.matrix(z)) z <- array(z, dim(inputData(data)))
      z <- (abs(inputData(data) - inputData(verification.data)) <= fuzz) | z
      z <- z & !is.na(z)
      in.up <- !apply(z,2, all)
     }
   if (is.null(seriesNamesOutput(data.names))) out.up <- NULL
   else
     {ld <-periodsOutput(data)
      lv <-periodsOutput(verification.data)
      l <- max(ld, lv)
      if (ld < l)
        outputData(data) <- ts(rbind(outputData(data), 
                                      matrix(NA,l-ld, nseriesOutput(data))),
                         start=start(data), frequency=frequency(data))
      if (lv < l)
        outputData(verification.data) <- ts(
                                rbind(outputData(verification.data), 
                                      matrix(NA,l-lv, nseriesOutput(data))),
                     start=start(outputData(verification.data)),
                     frequency=frequency(verification.data))
      z <- ( is.na(outputData(data)) & is.na(outputData(verification.data)))    # both NA
    # next fixes an Splus bug (IMHO) that the dim is dropped for col matrix
      if (!is.matrix(z)) z <- array(z, dim(outputData(data)))
      z <- (abs(outputData(data) - outputData(verification.data)) <= fuzz) | z
      z <- z & !is.na(z)
      out.up <- !apply(z,2, all)
     }
   list(any(c(in.up,out.up)), input=in.up, output=out.up, data=data)   
  }

checkForFileDateChanges <- function(data.names, verification.dates)
  {# check file dates against previous dates
   # It is preferable to do file date checking with a Unix shell script rather 
   #   than in S, and then start S for further checks only when the time stamp
   #   on the database files has changed.
   up.in <-NULL
   if (!is.null(seriesNamesInput(data.names)))
    {for (f in data.names$input$db) up.in <- c(up.in, file.date.info(f))
     inT <-any(verification.dates$input != up.in)
    }
   up.out <-NULL
   for (f in data.names$output$db) up.out <- c(up.out,file.date.info(f))
   outT <-any(verification.dates$output != up.out)
   list( any(c(inT,outT)), input=inT, output=outT, 
         dates=list(input=up.in, output=up.out))
  }



simpleMonitoring <- function(model, data.names, 
   previous.data=NULL,
   mail.list=NULL,
   error.mail.list=Sys.info()[["user"]],
   message.title="Simple Monitoring",
   message.subject="Simple Monitoring",
   message.footnote=NULL,
   show.start= c(0,-3),
   show.end  = c(0,12),    
   report.variables= seriesNames(data.names),
   data.sub.heading=NULL,
   data.tag=" ",
   forecast.tag="f",
   run.again=FALSE,
   save.as=NULL)

{# Step 0 -  prepare message files and error checking
    error.message <- c(message.title, paste(dateParsed(), collapse="."),
              "An error condition occurred running simpleMonitoring.",
              "The message.file at the time of the error follows:") 
    message <- ""     
    on.exit(Sys.mail(error.mail.list,
                 subject=paste("error ",message.subject),
                 body= c(error.message, message)))
    if ( dseclass(model)[1] == "TSestModel" ) model <- TSmodel(model)
    if (!is.null(data.names$pad.end))
       {if(!data.names$pad.end)
          warning("pad.end in data definition may disable retrieving all data.")
       } 
    else if (!is.null(data.names$pad))
       {if(!data.names$pad)
          warning("pad in data definition may disable retrieving all data.")
       } 

# The following line is useful for debugging
#Sys.mail(error.mail.list, subject=paste("checking ",message.subject), 
#                         body=paste(dateParsed(), collapse="."))

 # Step 1 - retrieve & check for updated data  or
 #            initialize system and if previous.data is NULL
    if (is.null(previous.data))
      {data <- freeze(data.names)
       message <- "Initializing simple monitoring:"   
       status <- "Simple monitoring initialized."   
      }
    else if (run.again)
      {data <-previous.data  
       status <- "Simple monitoring re-run."   
      }
    else
      {updated.data<-checkForValueChanges(data.names,
                           verification.data=previous.data,
                           discard.current=TRUE)
       if(updated.data[[1]])
         {data <-updated.data$data
	  # updated.data$input & $output are logical vector so don't use
	  #    inputData() and outputData() in next
          message <- c("data updates: ", 
               seriesNamesInput(data)[updated.data$input],
              seriesNamesOutput(data)[updated.data$output])
          status <- "Simple monitoring updated."   
         }
       else
         {on.exit()
          return(invisible(list(data=previous.data, 
                status="Simple monitoring updates not necessary.")))
         }
      }

 # Step 2 - check data
   # Sometimes data is available as soon as there are any days in a month (with
   #   ignore on in Fame). The following 4 lines trim these, but that may not be
   #   the best way to handle them.
   year.mo <- c(dateParsed()$y,dateParsed()$m) - c(0,1)
   data  <- tfwindow(data,  end=year.mo, warn=FALSE )

 # Step 3 - run forecast
   pred<-forecast(model, data)$forecast[[1]]
   pred <-splice.tagged(outputData(data), pred, tag1=data.tag,tag2=forecast.tag) 
 
 # Step 4 - generate report and mail
    message <-c(message,"The forecasts are now:")
    #starting and end period for plots & printing:
    start.<-(endOutput(data)+show.start) 
    end.  <-(endOutput(data)+show.end)

    report.variables$input<- 
            (report.variables$input == seriesNamesInput(data.names))
    report.variables$output<- 
            (report.variables$output == seriesNamesOutput(data.names))
#    rv <- tagged(selectSeries(pred, series=report.variables$output),
#                 tags= (attr(pred,"tags")) [,report.variables$output, drop=FALSE])
    rv <- selectSeries(pred, series=report.variables$output)
#    tframe(rv) <- tframe(pred)
#    inp <-tagged(inputData(data)[,report.variables$input, drop=FALSE],tags= data.tag)
#   should probably have data <- tagged(data, tags=data.tag) at the beginning
#    but 
    inp <-tagged(inputData(data),tags= data.tag)
    inp <-selectSeries(inp, report.variables$input)
#    tframe(inp) <-  tframe(inputData(data))
    rv <- tfwindow(tbind( inp, rv), start=start., end=end., warn=FALSE)   
    message <- c(message,fprint(rv, digits=5, sub.title=data.sub.heading)) 

    if (!is.null(message.footnote)) message <-c(message, message.footnote)
    if (!is.null(mail.list)) 
        Sys.mail(mail.list, subject=message.subject, body= message)

 # Step 4 - clean-up
    if (!is.null(save.as)) 
       assign(save.as,list(model=model, data=data, pred=pred), where=1)
    on.exit()
    #return latest data for comparison next time. Note - the forecast is NOT
    # returned (but may be saved above).
    invisible(list(data=data, status=status, message=message)) 
}

watch.data <- function(data.names, 
   previous.data=NULL,
   mail.list="gilp",
   error.mail.list=NULL,
   message.title="Data Monitor\n",
   message.subject="Data Monitor",
   message.footnote=NULL)

{# monitors data bases and check series for changes with e-mail of results.
 # this should be used with a script which watches for file date changes.
 #  ( see example in the file watch.data.readme)
 # data.names is a TSdata (names) object.
 # mail.list and error.mail.list should be single strings (not vectors)
 # If mail.list is null then mail is not sent (useful for testing).
 #  but the string can contain multiple user ids for mail
 # previous.data must normally be supplied. If it is not (ie. if it is NULL)
 # then the system will be initiallized and the returned result will be
 # the previous.data for the next time the function is called.

 # Step 0 - prepare message files 
    error.message <- c(message.title, paste(dateParsed(), collapse="."),
              "An error condition occurred running watch.data.",
              "The message.file at the time of the error follows:") 
    message <- ""     
    on.exit(Sys.mail(error.mail.list, subject=paste("error ",message.subject),
                 body= c(error.message, message)))

 # Step 1 - retrieve & check for updated data 

    data.names <- modify.TSPADIdata(data.names, pad.end=TRUE)
    #  Initialize system and exit if previous.data is NULL
    if (is.null(previous.data))
      {current.data <- freeze(data.names)
       on.exit()
       #return latest data for comparison next time
       return(invisible(list(data=current.data,
           status="System watch.data initialized."))) 
      }
    update<-checkForValueChanges(data.names,
                           verification.data=previous.data$data,
                           discard.current=FALSE)
    if (!update[[1]] )
        {on.exit()
         return(invisible(list(data=previous.data$data, 
             status="No data updates.")))
        }
    else
       message <- c(message, "data updates: ", 
              seriesNamesOutput(update$data)[update$output],)

 # Step 2 - mail 
    if(!is.null(message.footnote)) message <- c(message,message.footnote)
    Sys.mail(mail.list, subject=message.subject, body= message)

 # Step 3 - clean-up
    on.exit()
    #return latest data for comparison next time
    invisible(list(data=update$data, status="Data has been updated.")) 
}

