\name{tfplot.combined.forecast}
\alias{tfplot.combined.forecast}
\title{Specific Methods for tfplot}
\description{See the generic function description.}
\usage{
    \method{tfplot}{combined.forecast}(obj,verbose=F, 
       start.=start(obj$data$output),
       Title="Projection", select.inputs  = NULL, select.outputs = NULL,
       mar=par()$mar)
}
\arguments{
    \item{obj}{object to be plotted.}
    \item{verbose}{logical indicating if additional information is provided.}
    \item{start.}{start of plot. (passed to tfwindow)}
    \item{Title}{title for plot.)}
    \item{select.inputs}{vector of integers or strings indicating inputs 
       to be plotted.}
    \item{select.outputs}{vector of integers or strings indicating outputs 
       to be plotted.}
    \item{mar}{margins passed to plot. See par.}
}

\seealso{ \code{\link[tframe]{tfplot}} }
%\keyword{DSE}
\keyword{ts}

