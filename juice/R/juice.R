
# This files contains routines for concentrating (prcomp/cancor) and 
#  reconstituting series which should contain similar information. 

# The object "TSdataconcentrate" is a TSdata object plus a reduction 
#  (proj) estimated with est.projection.


############################################################################

#  functions for concentrating and reconstituting data 
# (using principal components canon. corr) 

############################################################################

# A copy of prcomponents may be added to R's mva??
prcomponents <- function(x, center=TRUE, scale=TRUE, N=nrow(x)-1)
   {if (center) center <- apply(x,2,mean)
    else        center <- rep(0, ncol(x))
    if (scale)  scale  <- sqrt(apply(x,2,var)) 
    else        scale  <- rep(1, ncol(x))
    s <- svd(sweep(sweep(as.matrix(x),2, center),2, scale, FUN="/"))
    # remove anything corresponding to effectively zero singular values.
    rank <- sum(s$d > (s$d[1]*sqrt(.Machine$double.eps)))
    if (rank < ncol(x)) s$v <- s$v[,1:rank, drop=FALSE]
    s$d <- s$d/sqrt(N)
    
#   r <- list(sdev=s$d, proj=s$v,x=x %*% s$v, center=center, scale=scale)
    r <- list(sdev=s$d, proj=s$v, center=center, scale=scale)
    r
}

cancorrelation <- function(x, y, xcenter=TRUE, ycenter=TRUE, xscale=TRUE, yscale=TRUE)
   {# scaling does not affect the result but is useful for backward calculation?
    if (! require("mva")) stop("cancorrelation requires library mva.")
    if (xcenter) xcenter <- apply(x,2,mean)
    else         xcenter <- rep(0, dim(x)[2])
    if (xscale)  xscale  <- sqrt(apply(x,2,var))
    else         xscale  <- rep(1, dim(x)[2])
    if (ycenter) ycenter <- apply(y,2,mean)
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

est.projection <- function(data, ...) {UseMethod("est.projection")}

est.projection.default <- function(data, center=TRUE, scale=TRUE, n=1)
 {# data should be a matrix. n is number of components to keep.
  data <- freeze(data)
  if (ncol(data) < n) stop("n cannot exceed columns of data.")
  conc<- prcomponents(data,center=center, scale=scale)
  conc$proj <- conc$proj[,1:n, drop=FALSE]
  classed(conc, "concentrator") # constructor
 }

est.projection.TSdata <- function(data, center=TRUE, scale=TRUE, m=1,p=1)
 {# use principle components if there is just input or just output, otherwise
  # otherwise use canonical correlation.
  data <- freeze(data)
  conc <- list()
  if ((0!= input.dimension(data)) & (0!=output.dimension(data)))
    {if (input.dimension(data) < m) stop("m cannot exceed input data dimension.")
     if (output.dimension(data)< p) stop("p cannot exceed output data dimension.")
     cc <- cancorrelation(input.data(data), output.data(data),
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
  else if (0!= input.dimension(data))
    {if (input.dimension(data) < m) stop("m cannot exceed input data dimension.")
     conc$input <- est.projection( input.data(data),center=center, scale=scale, n=m)
    }
  else if (0!=output.dimension(data))
    {if (output.dimension(data) < p) stop("p cannot exceed output data dimension.")
     conc$output<- est.projection(output.data(data),center=center, scale=scale, n=p)
    }
  classed(conc, "TSdataconcentrator") # constructor
 }



est.concentrate <- function(data, m=1,p=1, center=TRUE, scale=TRUE)
 {warning("est.concentrate is depreciated. Use concentrate(data, conc=est.projection(...))")
  concentrate(data, conc=est.projection(data, m=m, p=p, center=center, scale=scale))
  }
 
concentrate <- function(d, ...) {UseMethod("concentrate")}

concentrate.default <- function(d, center=TRUE,   scale=TRUE,  n=1, conc=NULL)
 {#conc=est.projection(d, center=center, scale=scale, n=n)) misses freeze
  # conc (projection) as returned by prcomp (see est.projection),
  d <- freeze(d)
  if (is.null(conc)) conc <- est.projection(d, center=center, scale=scale, n=n)
  else conc <- concentrator(conc)
  attr(d, "concentrator") <- conc
#  attr(newd, "original") <- d
  classed(d, c("concentrate", dseclass(d))) #constructor concentrate.default
 }

concentrate.TSdata <- function(d, center=TRUE, scale=TRUE,  m=1, p=1, conc=NULL)
 {#conc=est.projection(d, center=center, scale=scale, m=m, p=p)) misses freeze
  d <- freeze(d)
  if (is.null(conc)) conc <- est.projection(d, center=center, scale=scale, m=m, p=p)
  else conc <- concentrator(conc)
  if (0!= input.dimension(d))
      input.data(d) <- concentrate(input.data(d), conc=conc$input)
  if (0!=output.dimension(d))
     output.data(d) <- concentrate(output.data(d), conc=conc$output)
  classed(d, c("TSdataconcentrate", "TSdata")) # constructor (concentrate.TSdata)
 }

is.TSdataconcentrate <- function(x) {inherits(x, "TSdataconcentrate")}
is.TSmodelconcentrate <- function(x) {inherits(x, "TSmodelconcentrate")}
is.concentrate <- function(x) {inherits(x, "concentrate")}
is.TSdataconcentrator <- function(x) {inherits(x, "TSdataconcentrator")}
is.concentrator <- function(x) {inherits(x, "concentrator")}

print.concentrate <- function(x, ...)
  {cat("Original data:")     ; print(concentrateOriginal(x, ...))
   cat("Concentrated data:") ; print(concentrateOnly(x, ...))
   invisible(x)
  }



concentrateOnly <- function(d, ...) {UseMethod("concentrateOnly") }

concentrateOnly.concentrate <- function(d)
 {#return concentrate (as simple data) with original and concentrator removed
  conc <- concentrator(d)
  newd <-  sweep(sweep(d,2, conc$center), 2, conc$scale, FUN="/") %*% conc$proj
  tframe(newd) <- tframe(d)
#  series.names(newd) <- paste("concentrate", seq(ncol(d)))
  series.names(newd) <- concentrated.series.names(d)
  attr(newd, "concentrator") <- conc
  attr(newd, "original") <- NULL
  classed(newd, dseclass(d)[-1])  # deconstructor
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


concentrateOriginal <- function(d, ...) {UseMethod("concentrateOriginal") }
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




concentrator <- function(d, ...) {UseMethod("concentrator") }#extract conc 
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




reconstitute <- function(d, conc=NULL, names=NULL) {UseMethod("reconstitute")}

TSdata.TSdataconcentrate <- reconstitute 

reconstitute.default <- function(d, conc, names=series.names(d))
 {# this actually concentrates the data and then reconstitutes it.
  conc <- concentrator(conc)
  inv <- svd(conc$proj) 
  newd <- freeze(d)
  newd <- newd %*% inv$v %*% sweep(t(inv$u), 1, 1/inv$d, "*")
  newd <-  sweep(sweep(newd,2, conc$scale, FUN="*"),2, -conc$center)
  tframe(newd) <- tframe(d)
  if (!is.null(names)) series.names(newd) <- paste("recon.", names)
  attr(newd, "concentrator") <- conc
  classed(newd, c("reconstitute", dseclass(d))) #constructor reconstitute.default
 }
  

reconstitute.concentrate <- function(d, conc=concentrator(d),
                                    names=series.names(d))
 {# d as returned by concentrate.
  newd <- reconstitute.default(concentrateOnly(d), conc=conc, names=names)
  attr(newd, "original") <- concentrateOriginal(d)
  newd #constructor reconstitute.concentrate
 }

reconstitute.TSdataconcentrate <- function(d, conc=concentrator(d),
                                    names=series.names(d))
 {# d as returned by concentrate. A different conc (proj) can be used.
  # don't use input.data(d) or input.data(concentrateOnly(d)) here as they both
  # obscure the fact that d is concentrate. ??but so what ? this comment may be old
  if (0!= input.dimension(d))
    input.data(d) <- reconstitute(d$input, conc$input, names=names$input) 
  if (0!=output.dimension(d)) 
   output.data(d) <- reconstitute(d$output, conc$output,names=names$output)
  classed(d, c("TSdatareconstitute", "TSdata"))# constructor reconstitute.TSdataconcentrate)
 }


is.TSdatareconstitute <- function(x) {inherits(x, "TSdatareconstitute")}


canonical.prediction <- function(d, conc=concentrator(d),
      q=min(concentrated.input.dimension(d), concentrated.output.dimension(d)))
 {if (all(q == 0))
    stop("Both input and output data must exist for canonical.prediction.")
  if (max(q) >
      min(concentrated.input.dimension(d), concentrated.output.dimension(d)))
     stop(paste("number of canonical coordinates cannot exceed min. of input and output canonical components (",
    min(concentrated.input.dimension(d), concentrated.output.dimension(d)),")"))
  prj <- conc$output
  if (1 == length(q)) q <- seq(q)
  if (all(q != seq(concentrated.output.dimension(d))))
    {prj$u <- prj$u[ , q, drop=FALSE]
     prj$v <- prj$v[q,  , drop=FALSE] 
     prj$d <- prj$d[q]
    }
  pred <- reconstitute(input.data(concentrateOnly(d))[,q], prj,
                           names=output.series.names(concentrateOriginal(d)))
  attr(pred, "original") <- output.data(concentrateOriginal(d))
  classed(pred, "TScanonical.prediction") #constructor (canonical.prediction)
 }

is.TScanonical.prediction <- function(x) {inherits(x, "TScanonical.prediction")}

concentrateOriginal.TScanonical.prediction <- function(d)
    {attr(d, "original")}   

tfplot.TScanonical.prediction <- function(x, xlab=NULL, ylab=NULL,
         graphs.per.page=5, start.=NULL, end.=NULL, series=seq(nseries(x)))
{# plot actual data and data reconstituted from canonical.prediction.
 z <- concentrateOriginal(x)
 series.names(z) <- series.names(x) # used on plot
 tfplot(z, x, xlab=xlab, ylab=ylab, graphs.per.page=graphs.per.page, 
         start.=start., end.=end., series=series)
 invisible()
}

	 	 
end.TScanonical.prediction <- function(obj){end(concentrateOriginal(obj))}


percent.change.TScanonical.prediction <-function (mat,...) {
	pchange <- percent.change(classed(mat, dseclass(mat)[-1]), ...)
	attr(pchange, "original") <- percent.change(attr(mat, "original"), ...)
	classed(pchange, "TScanonical.prediction") #re constructor (percent.change)
}
  

######################################################################

#    model estimation

######################################################################


est.concentrated.model <- function(data, estimation="est.VARX.ls",
                                         estimation.args=NULL, ...)
     {UseMethod("est.concentrated.model")}

est.concentrated.model.TSdata <- function(data, estimation="est.VARX.ls",  
                                                estimation.args=NULL,
					m=1,p=1, center=TRUE, scale=TRUE) 
  {d <- concentrate(data, conc=est.projection(data, m=m, p=p,
                                        center=center, scale=scale))
   est.concentrated.model(d, estimation=estimation, 
                             estimation.args=estimation.args)
  }

est.concentrated.model.TSdataconcentrate <- function(data, 
     estimation="est.VARX.ls", 
     estimation.args=NULL, warn=TRUE) 
  {d <- concentrateOnly(data)
   m <-TSmodel(do.call(estimation, append(list(d), estimation.args)))
# next should be   concentrator(m) <- concentrator(data)
   m$conc <- concentrator(data)
   dseclass(m) <- c("TSmodelconcentrate", dseclass(m))
   series.names(m) <- series.names(data)
   l(m, data, warn=warn)
  }


is.TSmodelconcentrate <- function(x) {inherits(x, "TSmodelconcentrate")}

l.TSmodelconcentrate <- function(model,data, sampleT=nrow(output.data(data)), 
                                  predictT=sampleT, result=NULL, warn=TRUE)
  {#model should be TSmodelconcentrate and data should be TSdataconcentrate
   if (!is.TSdataconcentrate(data)) stop("data should be TSdataconcentrate.")
   if (!is.TSmodelconcentrate(model)) stop("model should be TSmodelconcentrate.")
   pred  <- l(concentrateOnly(model), concentrateOnly(data),
              sampleT=sampleT, predictT=predictT)$estimates$pred
   pred  <- reconstitute(pred, concentrator(data)$output,
        names=output.series.names(data))

   if((!is.null(result)) && (result == "pred")) return(pred)
   r <- residual.stats(pred, output.data(concentrateOriginal(data)),
                       sampleT=sampleT, warn=warn)
   if(is.null(result)) 
       return(classed(list(estimates = r, data = data, 
           model = model), "TSestModel"))
   else if(result == "like") return(r$like[1])# neg.log.like.
   else return(r[[result]]) 
   stop("should never get to here.")
  }
  

check.consistent.dimensions.TSmodelconcentrate <- function(model, data=NULL)
  {m <- classed(model, dseclass(model)[-1]) # deconstructor 
   check.consistent.dimensions(m)
   if(!is.null(data))
       {if(output.dimension(model) != output.dimension(data))
	    stop("Model and data output dimensions do not correspond.")
        if( input.dimension(model) != input.dimension(data))
	    stop("Model and data input dimensions do not correspond.")
       }
   return(T)
  }

     
   
##########################################################

#    methods for generics

##########################################################




input.dimension.TSmodelconcentrate <- function(m)
     {d <-nrow(concentrator(m)$input$proj);  if(is.null(d)) 0 else d}
output.dimension.TSmodelconcentrate <- function(m)
     {d <-nrow(concentrator(m)$output$proj); if(is.null(d)) 0 else d}



concentrated.dimension <- function(x){UseMethod("concentrated.dimension")}
   
concentrated.dimension.concentrate <- function(x) {ncol(concentrator(x)$proj)}



concentrated.input.dimension <- function(x) {concentrated.dimension(x$input)}
concentrated.output.dimension <- function(x) {concentrated.dimension(x$output)}



concentrated.series.names <- function(x){UseMethod("concentrated.series.names")}
   
concentrated.series.names.concentrate <- function(x)
  { paste("concentrate", seq(concentrated.dimension(x))) }

concentrated.series.names.TSdata <- function(x)
   {list(input = concentrated.input.series.names(x),
        output = concentrated.output.series.names(x)) }

concentrated.input.series.names <- function(x)
 {if(is.null(x$input)) NULL else concentrated.series.names(x$input)}

concentrated.output.series.names <- function(x)
 {if(is.null(x$output)) NULL else concentrated.series.names(x$output)}




settf.concentrate <- function(value, x)
 {# this would not be necessary if tframe inheritence was more cleanly 
  # separated from classes of x. See "pure" comments in settf.default
  cls <- dseclass(x)
  classed(x, settf(value, classed(x, cls[-1])), cls)
 }


tfprint.concentrate <- function(x, ...)
  {cat("Original data:")                  ; tfprint(concentrateOriginal(x, ...))
   cat("Data reconstituted from concentrate:") ; tfprint(reconstitute(x), ...)
   invisible(x)
  }



tfwindow.concentrate <- function(x, ...)
 {conc <- concentrator(x)
  cls <- dseclass(x)
  x <- tfwindow(classed(x, cls[-1]), ...)  # NextMethod might? work
  attr(x, "concentrator") <- conc  # kludge Rbug attr gets lost ??
  classed(x, cls)
 }


select.series.concentrate <-function (x, series = seq(nrow(concentrator(x)$proj))) 
 {conc <- concentrator(x)
  cls <- dseclass(x)
  orig <- select.series(concentrateOriginal(x), series=series)
  # look at select.series.default to get series right
  conc$sdev <- conc$sdev[series]
  conc$proj <- conc$proj[series,,drop=F]
  conc$center <- conc$center[series]
  conc$scale <- conc$scale[series]
  attr(x, "concentrator") <- conc  
  classed(x, cls)
 }





tfplot.concentrated <- function(x,...){stop("defunct. Use concentrated.tfplot")}

concentrated.tfplot <- function(x, ...) {tfplot(concentrateOnly(x), ...)}

   
tfplot.concentrate <- function(x, xlab=NULL, ylab=NULL,
         graphs.per.page=5, start.=NULL, end.=NULL, series=seq(nseries(x)))
{# plot actual data and data reconstituted from concentrate.
 tfplot(concentrateOriginal(x), reconstitute(x),  
         xlab=xlab, ylab=ylab, graphs.per.page=graphs.per.page, 
         start.=start., end.=end., series=series)
}

tfplot.TSdataconcentrate <- function(x, start. = NULL, end. = NULL,
    Title = "", reset.screen = T, 
    select.inputs  = seq(length = input.dimension(x)),
    select.outputs = seq(length = output.dimension(x)), 
    graphs.per.page = 5, ylab = NULL)
{# plot actual data and data reconstituted from concentrate.
 tfplot.TSdata(concentrateOriginal(x), reconstitute(x), start.=start.,end.=end., 
    Title = Title, reset.screen = reset.screen, 
    select.inputs  = select.inputs, select.outputs = select.outputs, 
    graphs.per.page = graphs.per.page, ylab = ylab)
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
  {if (0==input.dimension(data)) plot2by2(output.data(data))
   else  plot2by2(tbind(input.data(data), output.data(data)))
  }

plot2by2.default <- function(data, pch=".", ...)
{  p <- ncol(data)
   old.par <- par(mfrow = c(p, p), mar = c(2.1, 4.1, 3.1, 0.1))
   on.exit(par(old.par))
   for (i in 1:p) 
     for (j in 1:p)  tfplot(data[,i], data[,j], pch=pch, ...)
   invisible()
}

check.residuals.TSdataconcentrate <- function (data, ...) 
   {warning("This residual is the difference between original and reconstituted data")
    invisible(check.residuals.TSdata(
      TSdata(output=output.data(reconstitute(data)) - 
                    output.data(concentrateOriginal(data))), ...))
   }
    
check.residuals.TSdatareconstitute <- function (data, ...) 
   {invisible(check.residuals.TSdata(as.TSdata(data), ...)) }


check.residuals.concentrated <- function (data, ...) 
        {stop("defunct. Use concentrated.check.residuals")}

concentrated.check.residuals <- function (data, ...) 
        {check.residuals(concentrateOnly(data), ...)}


