\name{tfplot.concentrate}
\alias{tfplot.concentrate}
\alias{tfplot.concentrated}
\alias{tfplot.TScanonical.prediction}
\alias{tfplot.TSdataconcentrate}
\alias{tfplot.TSdatareconstitute}
\title{tfplot Specific Methods}
\description{See the generic function description.}
\usage{
    \method{tfplot}{concentrate}(x, xlab=NULL, ylab=NULL,
         graphs.per.page=5, start.=NULL, end.=NULL, 
	 series=seq(nseries(x)), mar=par()$mar)
    \method{tfplot}{TScanonical.prediction}(x, xlab=NULL, ylab=NULL,
         graphs.per.page=5, start.=NULL, end.=NULL, 
	 series=seq(nseries(x)), mar=par()$mar)
    \method{tfplot}{TSdataconcentrate}(x, start. = NULL, end. = NULL,
       Title = "", reset.screen = T, 
       select.inputs  = seq(length = input.dimension(x)),
       select.outputs = seq(length = output.dimension(x)), 
       graphs.per.page = 5, ylab = NULL, mar=par()$mar)
    \method{tfplot}{TSdatareconstitute}(x, ...)
}

\seealso{ \code{\link{tfplot}} }

%\keyword{DSE}
\keyword{ts}

