\name{tfplot.concentrate}
\alias{tfplot.concentrate}
\alias{tfplot.TScanonicalPrediction}
\alias{tfplot.TSdataconcentrate}
\alias{tfplot.TSdatareconstitute}

\title{tfplot Specific Methods}
\description{See the generic function description.}
\usage{
    \method{tfplot}{concentrate}(x, 
         tf=NULL, start=tfstart(tf), end=tfend(tf),
         series=seq(nseries(x)), 
	 Title=NULL, xlab=NULL, ylab=NULL,
	 graphs.per.page=5, mar=par()$mar, reset.screen=TRUE, ...)
    \method{tfplot}{TScanonicalPrediction}(x,
         tf=NULL, start=tfstart(tf), end=tfend(tf),
	 series=seq(nseries(x)),
	 Title=NULL, xlab=NULL, ylab=NULL,
         graphs.per.page=5, mar=par()$mar, reset.screen=TRUE, ...)
    \method{tfplot}{TSdataconcentrate}(x,
         tf=NULL, start=tfstart(tf), end=tfend(tf),
         select.inputs  = seq(length = nseriesInput(x)),
         select.outputs = seq(length = nseriesOutput(x)), 
         Title = NULL, xlab = NULL, ylab = NULL, 
         graphs.per.page = 5, mar=par()$mar, reset.screen =TRUE, ...)
    \method{tfplot}{TSdatareconstitute}(x, ...)
}

\arguments{
    \item{x}{an object to plot.}
    \item{start}{see the generic \code{tfplot}.}
    \item{end}{see the generic \code{tfplot}.}
    \item{tf}{see the generic \code{tfplot}.}
    \item{series}{see the generic \code{tfplot}.}
    \item{select.inputs}{see the generic \code{tfplot}.}
    \item{select.outputs}{see the generic \code{tfplot}.}
    \item{Title}{see the generic \code{tfplot}.}
    \item{xlab}{see the generic \code{tfplot}.}
    \item{ylab}{see the generic \code{tfplot}.}
    \item{graphs.per.page}{see the generic \code{tfplot}.}
    \item{mar}{see the generic \code{tfplot}.}
    \item{reset.screen}{see the generic \code{tfplot}.}
    \item{...}{arguments to be passed to other methods.}
}

\seealso{ \code{\link[tframe]{tfplot}} }

\concept{DSE}
\keyword{ts}

