
#  Functions in this file are defined for S (and differently in R).

is.S <- function(){is.Svanilla() | is.Splus() }
is.Svanilla <- function(){!exists("version")}
is.Splus <- function(){exists("version") && is.null(version$language)}
is.Splus.pre3.3 <- function()
   { ## <= 3.2
    if( !exists("version")) 
         { #-- `Vanilla' S (i.e. here "S version 4")
           #- this now works for  S version 4  (this is not S-plus 4.0, maybe 
           #             part of S-plus 5.0 !):
           lv <- nchar(Sv <- Sversion())
           r <- list(
	      major = substring(Sv, 1,1),
	      minor = substring(Sv, lv,lv))
         }
   else 
     {r <- version
      r$minor <- as.numeric(r$minor)
      r$major <- as.numeric(r$major)
     }
    is.Splus() &&  ((r$major+.1*r$minor) < 3.3)
   }


   is.R <- function() {F}

    # faking R compatibility:
    
    require <- function(package, ...)
      {if(package == "padi")
           {if("" == Sys.getenv("PADI_LDLIB")) return(FALSE)
	    else if (! exists("load.padi"))    return(FALSE)
	    load.padi()  
	    return(TRUE)}
       else TRUE 
      }
        
    Sys.getenv <- getenv

    Platform <- function()
      {list( OS.type = if(charmatch("MS Windows", version$os, nomatch=0)) "windows" else
       if(exists("unix")) "unix"  else "unkown")
      }
   
    Sys.info <- function()
     {if( !exists("version")) 
         { #-- `Vanilla' S (i.e. here "S version 4")
           #- this now works for  S version 4  (this is not S-plus 4.0, maybe 
           #             part of S-plus 5.0 !):
           lv <- nchar(Sv <- Sversion())
           r <- list(
	      major = substring(Sv, 1,1),
	      minor = substring(Sv, lv,lv))
         }
       else 
         {r <- version
          r$minor <- as.numeric(r$minor)
          r$major <- as.numeric(r$major)
         }
       if      (is.Splus())    r$language <- "default"
       else if (is.Svanilla()) r$language <- "S"
       r$sysname <- system("uname -s", intern=T) 
       r$release <- system("uname -r ", intern=T)
       r$version <- system("uname -v", intern=T)
       r$nodename  <- system("uname -n", intern=T)
       r$machine  <- system("uname -m", intern=T)
       r$login     <- system("whoami", intern=T) # not always correct
       if(0 == length(r$login)) r$login <- "unknown"
       r$user      <- Sys.getenv("USER")
       if(0 == length(r$user)) r$user <- "unknown"
       class(r) <- "Sys.info"
       r
      }
      

    .Syskern.options <- list(mailer = if ( "unix" == Platform()$OS.type )
                               "mailx" else "none" )
    getOption <- function(x) {.Syskern.options[[x]]}

   
    .SPAWN <- TRUE

   if (!exists("unlink.Splus")) unlink.Splus <- unlink
   unlink <- function(x, recursive = FALSE) unlink.Splus(x)

  # in S these only work in Unix so far.
    if("unix" == Platform()$OS.type) {
      system <- function(cmd, intern=F) {unix(cmd, output.to.S = intern)}
                 
      file.exists <- function(file) {1 == unix(paste("if [ -f ", file, 
           " ] ; then (echo 1) ; else (echo 0); fi"))}
      dateParsed <- function() 
        {d <- parse(text=system("date '+%Y %m %d %H %M %S'", intern=T),white=T)
         list(y=  eval(d[1]),
              m=eval(d[2]),
              d= eval(d[3]),
              H= eval(d[4]),
              M= eval(d[5]),
              S= eval(d[6]),
              tz=system("date '+%Z'", intern=T))
        }

      getwd <- function(){system("pwd", intern=T)} 
      Sys.sleep <- function(n) {system(paste("sleep ", n))} 
      file.copy <- function(from, to)system(paste("cp ", from, to)) 
      file.date.info <- function(file.name)
     	{# This could be a lot better. It will fail for files older than a year.
      	# Also, a returned format like date() below would be better.
      	mo <- (1:12)[c("Jan","Feb","Mar","Apr","May", "Jun","Jul","Aug", "Sep",
         "Oct","Nov","Dec") ==
          substring(system(paste("ls -l ",file), intern=T),33,35)]
      	day <- as.integer(substring(system(paste("ls -l ",file.name), intern=T),37,38))
      	hr  <- as.integer(substring(system(paste("ls -l ",file.name), intern=T),40,41))
      	sec <- as.integer(substring(system(paste("ls -l ",file.name), intern=T),43,44))
      	c(mo,day,hr,sec)
     	}
      }

    if("windows" == Platform()$OS.type)
      {file.date.info <- Sys.sleep <- system <-
         getwd <- file.copy <- function(cmd) 
         {stop("syskern system calls not implemented for this operating system.")}
      }
    
data <- function(x,  package = NULL, lib.loc = .path.package(package)) 
  {# possibly getenv("SHOME")
   #this is a bit of a kludge
   source(paste(lib.loc, "/data/", x, ".R", sep=""))
   if (exists(x))
     {# this is a kludge which works only if x is also the name of the object
      # assigned in the sourced file.
      r <- get(x)
      if (!is.null(class(r))) if(any(class(r) == "ts")) class(r) <- NULL
      if (!is.null(class(r$output)))
        if(any(class(r$output) == "ts")) class(r$output) <- NULL
      if (!is.null(class(r$input)))
        if(any(class(r$input)  == "ts")) class(r$input)  <- NULL
      assign(x,r, where=1)
     }
   invisible()
  }

.path.package <- function(pac) 
  {#this is a bit of a kludge
   paste(getenv("S_LIBS"),"/", pac, sep="")
  }

La.eigen <- eigen # fake La.eigen
La.svd <- function(x, nu = min(n, p), nv = min(n, p)) {
       # fake La.svd using svd. This will not give exactly the same result as
       # R's La.svd. It is probably closer to R's svd. 
       #  See R documention for details.
       x <- as.matrix(x)
       n <- nrow(x)
       p <- ncol(x)
       r <- svd(x, nu=nu, nv=nv)
       r$vt <- if (is.complex(r$v)) Conj(t(r$v)) else t(r$v)
       r$v <- NULL
       r
       }

if (!exists("par.Splus")) 
   {par.Splus <- par
    par <- function(no.readonly=FALSE, ...) par.Splus(...)
   }

help.start.DSE <- function(browser = "netscape") 
   {# this is really defunct, and not used in R, but I don't yet have
    #  an alternative for html help in S (and this may be broken too)
    #browser = "netscape -mono"
    #browser = "mosaic"
    #redirecting standard input and output is necessary in order to return to S 
    if(is.R()) help.start() else if ("unix" == Platform()$OS.type())
    system(paste(browser, " ", .path.package("dse1"),"/help/dsehome.htm </dev/null >/dev/null &", sep=""))
    invisible()
   }
