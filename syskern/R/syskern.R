##########

# Functions for S should be loaded from the S/ directory before 
#    functions in this file.

###########

#  Functions are for identifying S or R and flavours. 
#  Functions depending on differences between S and R
    is.S <- is.Svanilla <- is.Splus <- is.Splus.pre3.3 <- function(){FALSE}

    file.date.info <- function(file)
     	  {# format of the result could be better. 
	   x <- as.POSIXlt(file.info(file)$mtime)
	   c(1+x$mon,x$mday,x$hour,x$sec) # as previously used. should be improved
     	  }
   
    date.parsed <- function() 
          {d <- as.POSIXlt(Sys.time())
           list(y = 1900 + d$year,   # as previously used. should be improved
              m = 1+ d$mon,
              d = d$mday,
              H = d$hour,
              M = d$min,
              S = d$sec,
              tz = attr(d, "tzone"))
          }
        #  syskern.rm() was previously unlink() which is now supported in R.
        #  Unfortunately the argument recursive is not supported in S and the
        #    R 1.2 default value of FALSE is a change from previous Unix versions of R
        #    and from S and causes problems the way it is used in DSE.
#    syskern.rm <- function(file) unlink(file, recursive = TRUE)
	
    .SPAWN <- FALSE

Sys.mail <- function(address = Sys.info()$user,
                     ccaddress  = NULL,
                     bccaddress = NULL,
		     subject    = NULL,
		     body= "no body",
                     method = getOption("mailer")) {
   if(!(is.character(address) && nchar(address)>0))
      stop("A character string address must be specified.")

   # The arguments to mail, mailx, and Mail are all the same, but a different
   # mailer will require that this first part be re-organized under the
   # specific method.
   file <- tempfile()
   on.exit(unlink(file))
   cat(body, file=file, append=FALSE, sep="\n")
   cmdargs <- paste(address, "<", file, "2>/dev/null")
	
   if(is.character(ccaddress) && nchar(ccaddress)>0) 
            cmdargs <- paste(" -c '", ccaddress, "' ",  cmdargs)

   if(is.character(bccaddress) && nchar(bccaddress)>0) 
            cmdargs <- paste(" -b '", bccaddress, "' ",  cmdargs)

   if(is.character(subject) && nchar(subject)>0) 
            cmdargs <- paste(" -s '", subject, "' ",  cmdargs)

   status <- 1
   if(method == "mailx") status <- system(paste("mailx", cmdargs)) else
   if(method == "mail") status <- system(paste("mail", cmdargs))   else 
   if(method == "Mail") status <- system(paste("Mail", cmdargs))   else {
	warning(paste("mail method ", method, " not supported.\n"))
	return(FALSE)
	}
   if(status > 0) {
	     warning(paste("sending email with ", method, " failed.\n"))
	     return(FALSE)
	     }
   TRUE
   }

 
