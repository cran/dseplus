\name{tfplot.combined.forecast}
\alias{tfplot.combined.forecast}
\title{Specific Methods for tfplot}
\description{See the generic function description.}
\usage{
    \method{tfplot}{combined.forecast}(x, 
       start.=start(x$data$output), end.=end(x$data$output),
       select.inputs=NULL, select.outputs=NULL,
       Title="Projection", xlab=NULL, ylab=NULL, 
       graphs.per.page=5, mar=par()$mar, verbose=FALSE )
}
\arguments{
    \item{x}{object to be plotted.}
    \item{start.}{start of plot. (passed to tfwindow)}
    \item{end.}{end of plot. (passed to tfwindow)}
    \item{select.inputs}{vector of integers or strings indicating inputs 
       to be plotted.}
    \item{select.outputs}{vector of integers or strings indicating outputs 
       to be plotted.}
    \item{Title}{title for plot.)}
    \item{xlab}{xlab for plot.)}
    \item{ylab}{ylab for plot.)}
    \item{mar}{margins passed to plot. See par.}
    \item{graphs.per.page}{integer indicating number of graphs to place on a
       page.}
    \item{verbose}{logical indicating if additional information is provided.}
}

\seealso{ \code{\link[tframe]{tfplot}} }
%\keyword{DSE}
\keyword{ts}

