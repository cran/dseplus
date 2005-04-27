
# This files contains routines for concentrating (prcomp/cancor) and 
#  reconstituting series which should contain similar information. 

# The object "TSdataconcentrate" is a TSdata object plus a reduction 
#  (proj) estimated with estProjection.


############################################################################

#  functions for concentrating and reconstituting data 
# (using principal components canon. corr) 

############################################################################


# prcomp does orthogonal proj while cancor gives orthonormal result
#   (so proj changes scaling and 
#reverse has to be inverse not just transpose !!!!, thus square


######################################################################

#    concentrate data

######################################################################


# Note that concentrated data is really the original data with information for
#  concentrating it. That way the original data is not lost and is available
#  for comparison (as in tfplot). The truely concentrated data is produced on
#  demand by concentrateOnly.

estProjection <- function(data, center=TRUE, scale=TRUE, ...) UseMethod("estProjection")

estProjection.default <- function(data, center=TRUE, scale=TRUE, n=1, ...)
 {#  (... further arguments, currently disregarded)
  # data should be a matrix. n is number of components to keep.
  data <- freeze(data)
  if (ncol(data) < n) stop("n cannot exceed columns of data.")
  center <- if (center) colMeans(data) else rep(0, ncol(data))
  scale  <- if (scale) sqrt(apply(data,2,var)) else  rep(1, ncol(data))
  pr<- prcomp(data, retx = FALSE, center=center, scale=scale,
              tol=sqrt(.Machine$double.eps))
  classed(list(sdev=pr$sdev,  proj=pr$rotation[,1:n, drop=FALSE],
     center=center, scale=scale), "concentrator") # constructor
 }

estProjection.TSdata <- function(data, center=TRUE, scale=TRUE, m=1,p=1, ...)
 {#  (... further arguments, currently disregarded)
  # use principle components if there is just input or just output, otherwise
  # otherwise use canonical correlation.
  cancorrelation <- function(x, y, xcenter=TRUE, ycenter=TRUE, xscale=TRUE, yscale=TRUE)
   {# scaling does not affect the result but is useful for backward calculation?
    if (! require("stats")) stop("cancorrelation requires library stats.")
    # above was mva rather than stats prior to R 1.9.1
    if (xcenter) xcenter <- colMeans(x)
    else         xcenter <- rep(0, dim(x)[2])
    if (xscale)  xscale  <- sqrt(apply(x,2,var))
    else         xscale  <- rep(1, dim(x)[2])
    if (ycenter) ycenter <- colMeans(y)
    else         ycenter <- rep(0, dim(y)[2])
    if (yscale)  yscale  <- sqrt(apply(y,2,var))
    else         yscale  <- rep(1, dim(y)[2])
    s <- cancor(sweep(sweep(x,2, xcenter),2, xscale, FUN="/"),
                sweep(sweep(y,2, ycenter),2, yscale, FUN="/"),
                xcenter=FALSE, ycenter=FALSE)
    r <- list(xcoef=s$xcoef, ycoef=s$ycoef, cor=s$cor,
              xcenter=xcenter, ycenter=ycenter, xscale=xscale, yscale=yscale)
    r
   }
  data <- freeze(data)
  conc <- list()
  if ((0!= nseriesInput(data)) & (0!=nseriesOutput(data)))
    {if (nseriesInput(data) < m) stop("m cannot exceed input data dimension.")
     if (nseriesOutput(data)< p) stop("p cannot exceed output data dimension.")
     cc <- cancorrelation(inputData(data), outputData(data),
                  xcenter=center, ycenter=center, xscale=scale, yscale=scale)
     conc$input$proj    <- cc$xcoef[ , 1:m, drop=FALSE] 
     conc$input$center  <- cc$xcenter
     conc$input$scale   <- cc$xscale
     dseclass(conc$input) <- "concentrator" # constructor
     conc$output$proj   <- cc$ycoef[ , 1:p, drop=FALSE]
     conc$output$center <- cc$ycenter
     conc$output$scale  <- cc$yscale
     dseclass(conc$output) <- "concentrator" # constructor
    }
  else if (0!= nseriesInput(data))
    {if (nseriesInput(data) < m) stop("m cannot exceed input data dimension.")
     conc$input <- estProjection( inputData(data),center=center, scale=scale, n=m)
    }
  else if (0!=nseriesOutput(data))
    {if (nseriesOutput(data) < p) stop("p cannot exceed output data dimension.")
     conc$output<- estProjection(outputData(data),center=center, scale=scale, n=p)
    }
  classed(conc, "TSdataconcentrator") # constructor
 }


concentrate <- function(d, conc=NULL, center=TRUE, scale=TRUE, ...)
   UseMethod("concentrate")

concentrate.default <- function(d, conc=NULL, center=TRUE, scale=TRUE, n=1, ...)
 {#  (... further arguments, currently disregarded)
  #conc=estProjection(d, center=center, scale=scale, n=n)) misses freeze
  # conc (projection) as returned by prcomp (see estProjection),
  d <- freeze(d)
  if (is.null(conc)) conc <- estProjection(d, center=center, scale=scale, n=n)
  else conc <- concentrator(conc)
  attr(d, "concentrator") <- conc
#  attr(newd, "original") <- d
  classed(d, c("concentrate", dseclass(d))) #constructor concentrate.default
 }

concentrate.TSdata <- function(d, conc=NULL, center=TRUE, scale=TRUE, m=1, p=1, ...)
 {#  (... further arguments, currently disregarded)
  #conc=estProjection(d, center=center, scale=scale, m=m, p=p)) misses freeze
  d <- freeze(d)
  if (is.null(conc)) conc <- estProjection(d, center=center, scale=scale, m=m, p=p)
  else conc <- concentrator(conc)
  if (0!= nseriesInput(d))
      inputData(d) <- concentrate(inputData(d), conc=conc$input)
  if (0!=nseriesOutput(d))
     outputData(d) <- concentrate(outputData(d), conc=conc$output)
  classed(d, c("TSdataconcentrate", "TSdata")) # constructor (concentrate.TSdata)
 }

is.TSdataconcentrate <- function(x) {inherits(x, "TSdataconcentrate")}
is.TSmodelconcentrate <- function(x) {inherits(x, "TSmodelconcentrate")}
is.concentrate <- function(x) {inherits(x, "concentrate")}
is.TSdataconcentrator <- function(x) {inherits(x, "TSdataconcentrator")}
is.concentrator <- function(x) {inherits(x, "concentrator")}

print.concentrate <- function(x, ...)
  {cat("Original data:")     ; print(concentrateOriginal(x), ...)
   cat("Concentrated data:") ; print(concentrateOnly(x, ...))
   invisible(x)
  }



concentrateOnly <- function(d) UseMethod("concentrateOnly")

concentrateOnly.concentrate <- function(d)
 {#return concentrate (as simple data) with original and concentrator removed
  conc <- concentrator(d)
  newd <-  sweep(sweep(d,2, conc$center), 2, conc$scale, FUN="/") %*% conc$proj
  tframe(newd) <- tframe(d) # newd will not (necessarily) have the class of d
#  seriesNames(newd) <- paste("concentrate", seq(ncol(d)))
  seriesNames(newd) <- concentratedSeriesNames(d)
  attr(newd, "concentrator") <- conc
  attr(newd, "original") <- NULL
#  classed(newd, dseclass(d)[-1])  # newd may not have correct structure for
#                                    this (e.g. tsp attr for a ts)
  newd
 } 
 
concentrateOnly.TSdataconcentrate <- function(d) 
 {# beware infinite recursion. Do not call input.* or output.*
  TSdata( input=if (is.null(d$input))  NULL else concentrateOnly(d$input),
         output=if (is.null(d$output)) NULL else concentrateOnly(d$output))
 }

concentrateOnly.TSdatareconstitute <- function(d) # deconstructor
 {# beware infinite recursion. Do not call input.* or output.*
  TSdata( input=if (is.null(d$input))  NULL else concentrateOnly(d$input),
         output=if (is.null(d$output)) NULL else concentrateOnly(d$output))
 }

concentrateOnly.TSmodelconcentrate <- function(d) 
 {d$conc <- NULL
  classed(d, dseclass(d)[-1]) # deconstructor
 }

concentrateOnly.TSestModel <- function(d) 
 {if (!is.TSmodelconcentrate(TSmodel(d)) | !is.TSdataconcentrate(TSdata(d)))
     stop("model and data must be concentrates.")
  TSmodel(d) <- concentrateOnly(TSmodel(d))
  TSdata(d)  <- concentrateOnly(TSdata(d))
  d
 }


concentrateOriginal <- function(d) {UseMethod("concentrateOriginal") }
concentrateOriginal.concentrate <- function(d)
 {attr(d, "concentrator") <- NULL
  classed(d, dseclass(d)[-1]) #deconstructor
 }   

concentrateOriginal.TSdataconcentrate <- function(d) # deconstructor
 {# beware infinite recursion. Do not call input.* or output.*
  TSdata( input=if (is.null(d$input))  NULL else concentrateOriginal(d$input),
         output=if (is.null(d$output)) NULL else concentrateOriginal(d$output))
 }

concentrateOriginal.TSdatareconstitute <- function(d) # deconstructor
 {# beware infinite recursion. Do not call input.* or output.*
  TSdata( input=if (is.null(d$input))  NULL else concentrateOriginal(d$input),
         output=if (is.null(d$output)) NULL else concentrateOriginal(d$output))
 }




concentrator <- function(d) {UseMethod("concentrator") }#extract conc 
concentrator.concentrate <- function(d) {attr(d, "concentrator")} 
concentrator.concentrator <- function(d) {d} 
concentrator.TSdataconcentrator <- function(d) {d} 

concentrator.TSdata <- function(d)
 {classed(list( input=if(is.null(d$input))  NULL else concentrator(d$input),
               output=if(is.null(d$output)) NULL else concentrator(d$output)),
	"TSdataconcentrator")} 

concentrator.TSmodelconcentrate <- function(d) {d$conc} 


######################################################################

#    reconstitute data

######################################################################



TSdata.TSdataconcentrate <- function(data, names=NULL, ...)
    reconstitute(data, conc=concentrator(data), names=names) 


reconstitute <- function(d, conc=NULL, names=NULL) {UseMethod("reconstitute")}

reconstitute.default <- function(d, conc=NULL, names=seriesNames(d))
 {# this actually concentrates the data and then reconstitutes it.
  if(is.null(conc)) stop("conc argument must be supplied to reconstitute.default.")
  conc <- concentrator(conc)
  inv <- La.svd(conc$proj) 
  newd <- freeze(d)
  newd <- newd %*% Conj(t(inv$v)) %*% sweep(t(inv$u), 1, 1/inv$d, "*")
  newd <-  sweep(sweep(newd,2, conc$scale, FUN="*"),2, -conc$center)
  tframe(newd) <- tframe(d)
  if (!is.null(names)) seriesNames(newd) <- paste("recon.", names)
  attr(newd, "concentrator") <- conc
  classed(newd, c("reconstitute", dseclass(d))) #constructor reconstitute.default
 }
  

reconstitute.concentrate <- function(d, conc=concentrator(d),
                                    names=seriesNames(d))
 {# d as returned by concentrate.
  newd <- reconstitute.default(concentrateOnly(d), conc=conc, names=names)
  attr(newd, "original") <- concentrateOriginal(d)
  newd #constructor reconstitute.concentrate
 }

reconstitute.TSdataconcentrate <- function(d, conc=concentrator(d),
                                    names=seriesNames(d))
 {# d as returned by concentrate. A different conc (proj) can be used.
  # don't use inputData(d) or inputData(concentrateOnly(d)) here as they both
  # obscure the fact that d is concentrate. ??but so what ? this comment may be old
  if (0!= nseriesInput(d))
    inputData(d) <- reconstitute(d$input,  conc=conc$input,  names=names$input) 
  if (0!=nseriesOutput(d)) 
   outputData(d) <- reconstitute(d$output, conc=conc$output, names=names$output)
  classed(d, c("TSdatareconstitute", "TSdata"))# constructor reconstitute.TSdataconcentrate)
 }


is.TSdatareconstitute <- function(x) {inherits(x, "TSdatareconstitute")}


canonical.prediction <- function(d, conc=concentrator(d),
      q=min(concentrated.nseriesInput(d), concentrated.nseriesOutput(d)))
 {if (all(q == 0))
    stop("Both input and output data must exist for canonical.prediction.")
  if (max(q) >
      min(concentrated.nseriesInput(d), concentrated.nseriesOutput(d)))
     stop(paste("number of canonical coordinates cannot exceed min. of input and output canonical components (",
    min(concentrated.nseriesInput(d), concentrated.nseriesOutput(d)),")"))
  prj <- conc$output
  if (1 == length(q)) q <- seq(q)
  if (all(q != seq(concentrated.nseriesOutput(d))))
    {prj$u <- prj$u[ , q, drop=FALSE]
     prj$v <- prj$v[q,  , drop=FALSE] 
     prj$d <- prj$d[q]
    }
  pred <- reconstitute(inputData(concentrateOnly(d), series=q), conc=prj,
                           names=seriesNamesOutput(concentrateOriginal(d)))
  attr(pred, "original") <- outputData(concentrateOriginal(d))
  classed(pred, "TScanonicalPrediction") #constructor (canonical.prediction)
 }

is.TScanonicalPrediction <- function(x) {inherits(x, "TScanonicalPrediction")}

concentrateOriginal.TScanonicalPrediction <- function(d)
    {attr(d, "original")}   

tfplot.TScanonicalPrediction <- function(x,
         tf=NULL, start=tfstart(tf), end=tfend(tf),
	 series=seq(nseries(x)),
	 Title=NULL, xlab=NULL, ylab=NULL,
         graphs.per.page=5, mar=par()$mar, reset.screen=TRUE, ...)
{#  (... further arguments, currently disregarded)
 # plot actual data and data reconstituted from canonical.prediction.
 z <- concentrateOriginal(x)
 seriesNames(z) <- seriesNames(x) # used on plot
 tfplot(z, x, start=start, end=end, series=series,
	 Title=Title, xlab=xlab, ylab=ylab,
	 graphs.per.page=graphs.per.page, mar=mar, reset.screen=reset.screen)
 invisible()
}

	 	 
start.TScanonicalPrediction <- function(x, ...){start(concentrateOriginal(x), ...)}
#  (... further arguments, currently disregarded)
end.TScanonicalPrediction <- function(x, ...){end(concentrateOriginal(x), ...)}
#  (... further arguments, currently disregarded)
periods.TScanonicalPrediction <- function(x){periods(concentrateOriginal(x))}
frequency.TScanonicalPrediction <- function(x, ...){frequency(concentrateOriginal(x), ...)}
#  (... further arguments, currently disregarded)


percentChange.TScanonicalPrediction <-function (obj, base=NULL, lag=1,
    cumulate=FALSE, e=FALSE, ...) {
	pchange <- percentChange(classed(obj, dseclass(obj)[-1]),
	    base=base, lag=lag, cumulate=cumulate, e=e)
	attr(pchange, "original") <- percentChange(attr(obj, "original"),
	    base=base, lag=lag, cumulate=cumulate, e=e)
	classed(pchange, "TScanonicalPrediction") #re constructor (percentChange)
}
  

######################################################################

#    model estimation

######################################################################


estConcentratedModel <- function(data, estimation="estVARXls",
                                         estimation.args=NULL, ...)
     UseMethod("estConcentratedModel")

estConcentratedModel.TSdata <- function(data, estimation="estVARXls",  
              estimation.args=NULL, m=1, p=1, center=TRUE, scale=TRUE, ...) 
  {#  (... further arguments, currently disregarded)
   d <- concentrate(data, conc=estProjection(data, m=m, p=p,
                                        center=center, scale=scale))
   estConcentratedModel(d, estimation=estimation, 
                             estimation.args=estimation.args)
  }

estConcentratedModel.TSdataconcentrate <- function(data, 
     estimation="estVARXls", 
     estimation.args=NULL, warn=TRUE, ...) 
  {#  (... further arguments, currently disregarded)
   d <- concentrateOnly(data)
   m <-TSmodel(do.call(estimation, append(list(d), estimation.args)))
# next should be   concentrator(m) <- concentrator(data)
   m$conc <- concentrator(data)
   dseclass(m) <- c("TSmodelconcentrate", dseclass(m))
   seriesNames(m) <- seriesNames(data)
   l(m, data, warn=warn)
  }


is.TSmodelconcentrate <- function(x) {inherits(x, "TSmodelconcentrate")}

l.TSmodelconcentrate <- function(obj1, obj2, sampleT=nrow(outputData(obj2)), 
                                  predictT=sampleT, result=NULL, warn=TRUE, ...)
  {#obj1 should be TSmodelconcentrate and obj2 should be TSdataconcentrate
   if (!is.TSdataconcentrate(obj2)) stop("obj2 should be TSodataconcentrate.")
   if (!is.TSmodelconcentrate(obj1)) stop("obj1 should be TSmodelconcentrate.")
   pred  <- l(concentrateOnly(obj1), concentrateOnly(obj2),
              sampleT=sampleT, predictT=predictT)$estimates$pred
   pred  <- reconstitute(pred, conc=concentrator(obj2)$output,
        names=seriesNamesOutput(obj2))

   if((!is.null(result)) && (result == "pred")) return(pred)
   r <- residualStats(pred, outputData(concentrateOriginal(obj2)),
                       sampleT=sampleT, warn=warn)
   if(is.null(result)) return(classed(list(estimates = r, data = obj2,
                                  model = obj1), "TSestModel"))
   else if(result == "like") return(r$like[1])# neg.log.like.
   else return(r[[result]]) 
   stop("should never get to here.")
  }
  

checkConsistentDimensions.TSmodelconcentrate <- function(obj1, obj2=NULL)
  {m <- classed(obj1, dseclass(obj1)[-1]) # deconstructor 
   checkConsistentDimensions(m)
   if(!is.null(obj2))
       {if(nseriesOutput(obj1) != nseriesOutput(obj2))
	    stop("Model and data output dimensions do not correspond.")
        if( nseriesInput(obj1) != nseriesInput(obj2))
	    stop("Model and data input dimensions do not correspond.")
       }
   return(TRUE)
  }

     
   
##########################################################

#    methods for generics

##########################################################




nseriesInput.TSmodelconcentrate <- function(x)
     {d <-nrow(concentrator(x)$input$proj);  if(is.null(d)) 0 else d}
nseriesOutput.TSmodelconcentrate <- function(x)
     {d <-nrow(concentrator(x)$output$proj); if(is.null(d)) 0 else d}



concentratedDimension <- function(x){UseMethod("concentratedDimension")}
   
concentratedDimension.concentrate <- function(x) {ncol(concentrator(x)$proj)}



concentrated.nseriesInput <- function(x) {concentratedDimension(x$input)}
concentrated.nseriesOutput <- function(x) {concentratedDimension(x$output)}



concentratedSeriesNames <- function(x){UseMethod("concentratedSeriesNames")}
   
concentratedSeriesNames.concentrate <- function(x)
  { paste("concentrate", seq(concentratedDimension(x))) }

concentratedSeriesNames.TSdata <- function(x)
   {list(input = concentratedSeriesNamesInput(x),
        output = concentratedSeriesNamesOutput(x)) }

concentratedSeriesNamesInput <- function(x)
 {if(is.null(x$input)) NULL else concentratedSeriesNames(x$input)}

concentratedSeriesNamesOutput <- function(x)
 {if(is.null(x$output)) NULL else concentratedSeriesNames(x$output)}


"tframe<-.concentrate" <- function(x, value){
  cls <- dseclass(x)
  x <- classed(x, cls[-1])
  tframe(x) <- value 
  # may not have class back the way it should be ???  
  classed(x, cls)
}

tframed.concentrate <- function(x, tf=NULL, names = NULL) 
{cls <- dseclass(x)
 classed(tframed(classed(x, cls[-1]), tf=tf, names=names), cls)
}

#settf.concentrate <- function(value, x)
# {# this would not be necessary if tframe inheritence was more cleanly 
#  # separated from classes of x. See "pure" comments in settf.default
#  cls <- dseclass(x)
#  classed(x, settf(value, classed(x, cls[-1])), cls)
# }


tfprint.concentrate <- function(x, ...)
  {cat("Original data:")                  ; tfprint(concentrateOriginal(x), ...)
   cat("Data reconstituted from concentrate:") ; tfprint(reconstitute(x), ...)
   invisible(x)
  }



tfwindow.concentrate <- function(x, tf=NULL, start=tfstart(tf), end=tfend(tf), warn=TRUE)
 {conc <- concentrator(x)
  cls <- dseclass(x)
  x <- tfwindow(classed(x, cls[-1]), tf=tf, start=start, end=end, warn=warn)  # NextMethod might? work
  attr(x, "concentrator") <- conc  # kludge Rbug attr gets lost ??
  classed(x, cls)
 }


selectSeries.concentrate <-function (x, series = seq(nrow(concentrator(x)$proj))) 
 {conc <- concentrator(x)
  cls <- dseclass(x)
  orig <- selectSeries(concentrateOriginal(x), series=series)
  # look at selectSeries.default to get series right
  conc$sdev <- conc$sdev[series]
  conc$proj <- conc$proj[series,,drop=FALSE]
  conc$center <- conc$center[series]
  conc$scale <- conc$scale[series]
  attr(x, "concentrator") <- conc  
  classed(x, cls)
 }


concentrated.tfplot <- function(x, ...) {tfplot(concentrateOnly(x), ...)}

   
tfplot.concentrate <- function(x, tf=NULL, start=tfstart(tf), end=tfend(tf),
         series=seq(nseries(x)), 
	 Title=NULL, xlab=NULL, ylab=NULL,
	 graphs.per.page=5, mar=par()$mar, reset.screen=TRUE, ...)
{#  (... further arguments, currently disregarded)
 # plot actual data and data reconstituted from concentrate.
 tfplot(concentrateOriginal(x), reconstitute(x),  
        start=start, end=end, series=series, 
        Title=Title, xlab=xlab, ylab=ylab, 
	graphs.per.page=graphs.per.page, mar=mar,reset.screen=reset.screen)
}

tfplot.TSdataconcentrate <- function(x, 
    tf=NULL, start=tfstart(tf), end=tfend(tf),
    select.inputs  = seq(length = nseriesInput(x)),
    select.outputs = seq(length = nseriesOutput(x)), 
    Title = NULL, xlab = NULL, ylab = NULL, 
    graphs.per.page = 5, mar=par()$mar, reset.screen = TRUE, ...)
{#  (... further arguments, currently disregarded)
 # plot actual data and data reconstituted from concentrate.
 #was tfplot.TSdata( 
 tfplot(concentrateOriginal(x), reconstitute(x), start=start,end=end, 
    select.inputs  = select.inputs, select.outputs = select.outputs, 
    Title = Title, xlab = xlab, ylab = ylab, 
    graphs.per.page = graphs.per.page, mar=mar, reset.screen = reset.screen)
 invisible()
}




tfplot.TSdatareconstitute <- function(x, ...)
{# plot data reconstituted from concentrate.
 tfplot( as.TSdata(x),  ...)
 invisible()
}


# this doesn't need to be with juice?

plot2by2 <- function(data, ...) {UseMethod("plot2by2")}

plot2by2.TSdata <- function(data, ...)
  {if (0==nseriesInput(data)) plot2by2(outputData(data))
   else  plot2by2(tbind(inputData(data), outputData(data)))
  }

plot2by2.default <- function(data, pch=".", ...)
{  p <- ncol(data)
   old.par <- par(mfrow = c(p, p), mar = c(2.1, 4.1, 3.1, 0.1), no.readonly=TRUE)
   on.exit(par(old.par))
   for (i in 1:p) 
     for (j in 1:p)  tfplot(data[,i], data[,j], pch=pch, ...)
   invisible()
}

checkResiduals.TSdataconcentrate <- function (obj, ...) 
   {warning("This residual is the difference between original and reconstituted data")
    invisible(checkResiduals.TSdata(
      TSdata(output=outputData(reconstitute(obj)) - 
                    outputData(concentrateOriginal(obj))), ...))
   }
    
checkResiduals.TSdatareconstitute <- function (obj, ...) 
   {invisible(checkResiduals.TSdata(as.TSdata(obj), ...)) }


checkResiduals.concentrated <- function (obj, ...) 
        {stop("defunct. Use concentrated.checkResiduals")}

concentrated.checkResiduals <- function (data, ...) 
        {checkResiduals(concentrateOnly(data), ...)}


