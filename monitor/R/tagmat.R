###########################################################################

# tagged data class  (matrix with a "tags" attribute)       <<<<<<<<<<<<

###########################################################################

tags <- function(x) {attr(x, "tags")}
 
"tags<-" <- function(x, value)   
  {if (is.null(value))
       {attr(x, "tags") <- NULL
        if (is.object(x)) class(x) <- class(x)[ class(x) != "tagged"]
        return(x)
       }
   if (length(value) == 1) value <- array(value, dim(x))
    # drop any extra attributes
   attributes(value) <- list(dim=attributes(value)$dim)
   attr(x, "tags") <- value 
   classed(x, c("tagged", dseclass(x))) # constructor ("tags<-")
  }

tagged <- function(x, tags)UseMethod("tagged")

tagged.default <- function(x, tags){
  tags(x) <- tags; x }


tagged.TSdata <- function(x, tags)
  {if(0 != nseriesInput(x))   inputData(tags(x)) <- tags$input
   if(0 != nseriesOutput(x)) outputData(tags(x)) <- tags$output
   x
  }

"tframe<-.tagged" <- function(x, value){
  cls <- dseclass(x)
  x <- classed(x, cls[-1])
  tframe(x) <- value 
  # may not have class back the way it should be ???  
  classed(x, cls)
}

selectSeries.tagged <- function(x, series=seq(ncol(x)))
     {names <- seriesNames(x)
      if (is.character(series)) series <- match(names,series, nomatch=0) > 0
      tagged(selectSeries.default(x, series=series),
             selectSeries.default(tags(x), series=series))
     }

tbind.tagged <- function(mat1, mat2)
{# aline and bind ts matrices and tags
 if (is.tagged(mat1)) tag1 <- tags(mat1)
 else                 tag1 <- array("mat1", dim(mat1))
 if (is.tagged(mat2)) tag2 <- tags(mat2)
 else                 tag2 <- array("mat2", dim(mat2))
 tframe(tag1) <- tframe(mat1)
 tframe(tag2) <- tframe(mat2)
 cls <- dseclass(mat1)
 # this should use NextMethod
 dseclass(mat1) <- dseclass(mat1)[-1]  # otherwise tbind calls this tbind
 if (0 == length(dseclass(mat1))) dseclass(mat1) <- NULL
 dseclass(mat2) <- dseclass(mat2)[-1]  # otherwise tbind calls this tbind
 if (0 == length(class(mat2))) dseclass(mat2) <- NULL
 tagged(classed(tbind(mat1, mat2), cls),  tbind(tag1,tag2)) 
}

is.tagged <- function(obj)  {inherits(obj,"tagged")}

testEqual.tagged <- function(obj1, obj2)
{ testEqual.matrix(obj1, obj2) & 
  testEqual.matrix(tags(obj1), tags(obj2))
}



fprint <- function(x, super.title=NULL, sub.title=NULL, 
        digits=options()$digits, space=" ", file=NULL, append=FALSE)
   {UseMethod("fprint")}

fprint.tagged <- function(x, super.title=NULL, sub.title=NULL, 
        digits=options()$digits, space=" ", file=NULL, append=FALSE) 
 {# Formattted print of a matrix of class tagged.
  # Corresponding characters are printed after matrix numbers.
  # A character matrix (out) is returned invisibly.
  # If file is not NULL then elements of out are printed to lines of the file.
  tags <- tags(x)
  out <- NULL
  f <- frequency(x)
  s <- start(x)
  s <- s[1] + (s[2]-1)/f
  if (12 ==f) p <- c("Jan","Feb","Mar","Apr","May", "Jun","Jul","Aug", "Sep",
         "Oct","Nov","Dec")
  else if (4 == f) p <- c("Q1","Q2","Q3","Q4")
  else if (52 == f) p <- format(1:52)
  else p <-NULL
  pre.space <- paste(rep(" ",nchar(format(s))+nchar(p[1])),collapse="")
  if (!is.null(super.title))  out <- paste(pre.space, super.title, sep="")
  names <- format(seriesNames(x), digits=digits)
  if (!is.null(names))
    {ot <- pre.space
     for (i in seq(length(names)))
        ot <- paste(ot, space,names[i],sep="")
     out <- c(out, ot)
    }
  if (!is.null(sub.title)) out <- c(out,paste(pre.space, sub.title,sep=""))
  m <- format(signif(x[,], digits=digits))
  for (i in seq(nrow(m))) 
    {d <- (s+(i-1)/f) +.Options$ts.eps # +eps or trunc sometimes gets wrong year
     ot <- paste(trunc(d)," ", p[round(1+f*(d%%1))]," ", sep ="")
     for (j in seq(ncol(m))) 
       {ot <-paste(ot,space, m[i,j], sep="")
        if (!is.null(tags)) ot <- paste(ot,tags[i,j], sep="")
       }
      out <- c(out, ot)
    }
  if (!is.null(file)) write(out, file=file, append=append)
  invisible(out)
 }


splice.tagged <- function(mat1, mat2, tag1=tags(mat1), tag2=tags(mat2), ...)
{#  (... further arguments, currently disregarded)
 # splice together 2 time series matrices as with splice.ts.
 # If data  is provided in both for a given period then mat1 takes priority.
 # The frequencies should be the same.
 # tag1 and tag2 are taken from mat1 and mat2 unless
 #   they are specified in the argument. If specified they
 #   should be single character strings or matrices of character 
 #   strings of same dimension as mat1 and mat2. This second is useful for multiple
 # applications of the function. The result is the
 # resulting spliced matrix of class "tagged"  
 # (suitable for use with fprint).
 # In the case tags are not available and are not specified 
 #   in the argument then they are set to "mat1" and "mat2".
 cls <- dseclass(mat1)
 if (is.null(tag1)) tag1 <- "mat1"
 if (is.null(tag2)) tag2 <- "mat2"
 if (length(tag1) == 1) tag1 <- array(tag1, dim(mat1))
 if (length(tag2) == 1) tag2 <- array(tag2, dim(mat2))
 if (is.null(mat1) & is.null(mat2)) return(NULL)
 if (is.null(mat2)) return(tagged(mat1, tag1))
 if (is.null(mat1)) return(tagged(mat2, tag2))
 freq <- frequency(mat1)
 if (freq != frequency(mat2)) stop("frequencies must be the same.\n")
 p <- dim(mat1)[2]
 if (p != dim(mat2)[2]) stop("number of series must be the same.\n")
 tframe(tag1) <- tframe(mat1)
 tframe(tag2) <- tframe(mat2)

 fr <- c(freq,1)
 st <- min(fr %*% start(mat1), fr %*% start(mat2))
 strt <- c(st %/% freq, st %% freq)
 en <- max(fr %*% end(mat1), fr%*% end(mat2))
 tf <- list(start=strt, frequency=freq)
 if (fr %*% start(mat1) > st) 
    {tag1 <-tframed(rbind(matrix("?", fr %*% start(mat1) -st, p), tag1),tf)
     mat1 <-tframed(rbind(matrix(NA,  fr %*% start(mat1) -st, p), mat1), tf)
    }
 if (fr %*%   end(mat1) < en) 
    {tag1 <-tframed(rbind(tag1, matrix("?", en - fr %*% end(mat1), p)), tf)
     mat1 <-tframed(rbind(mat1, matrix(NA,  en - fr %*% end(mat1), p)), tf)
    }
 if (fr %*% start(mat2) > st) 
    {tag2 <-tframed(rbind(matrix("?", fr %*% start(mat2) -st, p), tag2), tf)
     mat2 <-tframed(rbind(matrix(NA,  fr %*% start(mat2) -st, p), mat2), tf)
    }
 if (fr %*%   end(mat2) < en) 
    {tag2 <-tframed(rbind(tag2,matrix("?", en - fr %*% end(mat2), p)), tf)
     mat2 <-tframed(rbind(mat2, matrix(NA, en - fr %*% end(mat2), p)), tf)
    }
 na <- is.na(mat1)
#browser()
 mat1[na]  <- mat2[na]
 tag1[na] <- tag2[na]
 dimnames(mat1) <-list(round(time(mat1),digits=3),dimnames(mat1)[[2]])
 tags(mat1) <- tag1
 classed(mat1, cls )
}

tfwindow.tagged <- function(x, start.=NULL, end.=NULL, tf=NULL, warn=TRUE)
{# window a ts matrix of class "tagged".
 # With the default warn=T warnings will be issued if no truncation takes
 #  place because start or end is outside the range of data.
 tags <- tags(x)
 dseclass(x) <- dseclass(x)[-1]
 if (0 == length(dseclass(x))) dseclass(x) <- NULL
 # The next line converts scalars tags to a matrix.
 if (length(tags) == 1) tags <- array(tags, dim(x))
 # The next lines converts missing tags to a matrix.
 if (length(tags) == 0)
   {tags <- array("", dim(x))
    if (warn) warning("missing tags converted to empty string.")
   }
 tframe(tags) <- tframe(x)
 # The following is complicated by the fact that some versions of window
 #    look for missing arguments.
 if (is.null(start.))
   {x   <- tfwindow(x  , end.=end., tf=tf, warn=warn)
    tags<- tfwindow(tags,end.=end., tf=tf, warn=warn)
   }
 else if (is.null(end.))
   {x   <- tfwindow(x  , start.=start., tf=tf, warn=warn)
    tags<- tfwindow(tags,start.=start., tf=tf, warn=warn)
   }
 else
   {x   <- tfwindow(x,   start.=start., end.=end., tf=tf, warn=warn)
    tags<- tfwindow(tags,start.=start., end.=end., tf=tf, warn=warn)
   }
 tagged(x, tags)
}


