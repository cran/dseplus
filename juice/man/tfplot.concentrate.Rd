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
       select.inputs  = seq(length = nseriesInput(x)),
       select.outputs = seq(length = nseriesOutput(x)), 
       graphs.per.page = 5, ylab = NULL, mar=par()$mar)
    \method{tfplot}{TSdatareconstitute}(x, ...)
}

\arguments{
    \item{x}{an object to plot.}
    \item{start.}{see the generic \code{tfplot}.}
    \item{end.}{see the generic \code{tfplot}.}
    \item{series}{see the generic \code{tfplot}.}
    \item{select.inputs}{see the generic \code{tfplot}.}
    \item{select.outputs}{see the generic \code{tfplot}.}
    \item{xlab}{see the generic \code{tfplot}.}
    \item{ylab}{see the generic \code{tfplot}.}
    \item{mar}{see the generic \code{tfplot}.}
    \item{Title}{see the generic \code{tfplot}.}
    \item{graphs.per.page}{see the generic \code{tfplot}.}
    \item{reset.screen}{see the generic \code{tfplot}.}
    \item{...}{arguments to be passed to other methods.}
}

\seealso{ \code{\link[tframe]{tfplot}} }

%\keyword{DSE}
\keyword{ts}

