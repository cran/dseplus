\name{curvature.TSestModel}
\alias{curvature.TSestModel}
\title{Curvature for a TSestModel}
\description{Calculate curvature for a TSestModel.}
\usage{
    \method{curvature}{TSestModel}(obj, warn=T, ...)
}

\arguments {
    \item{emodel}{A TSestModel.}
    \item{warn}{Print additional warning messages.}
}

\details{
    None yet.
}
\seealso{
    \code{\link{curvature}}
}
\examples{
   if(is.R()) data("eg1.DSE.data.diff", package="dse1")
   model <- est.VARX.ls(TSdata( output=output.data(eg1.DSE.data.diff, series=1:2)),
                     max.lag=2)
   curvature(genD(model)) 
}
%\keyword{DSE}
\keyword{ts}

