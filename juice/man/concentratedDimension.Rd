\name{concentratedDimension}
\alias{concentratedDimension}
\alias{concentratedDimension.concentrate}

\title{Concentrated Dimension}
\description{
    The dimension (number of series) in concentrate data. This is the
    dimension onto which the original series has been projected.
}
\usage{
    concentratedDimension(x)
    \method{concentratedDimension}{concentrate}(x)
}
\arguments{
    \item{x}{a concentrated data object.}
}
\value{
    Depends on the argument. For a simple concentrated data object 
    the result is an integer.
}

\seealso{
    \code{\link{concentrated.nseriesInput}}
    \code{\link{concentrated.nseriesOutput}}
    \code{\link{concentrate}}
}
\examples{
    data("eg1.DSE.data", package="dse1")
    require("stats")
    z <- concentrate(eg1.DSE.data)
    concentratedDimension(outputData(z))
    concentrated.nseriesOutput(z)
}
\concept{DSE}
\keyword{ts}

