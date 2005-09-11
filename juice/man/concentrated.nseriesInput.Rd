\name{concentrated.nseriesInput}
\alias{concentrated.nseriesInput}
\alias{concentrated.nseriesOutput}
\title{Concentrated Dimension of TSdata}
\description{
    The dimension (number of series) in concentrate data. This is the
    dimension onto which the original series has been projected.
}
\usage{
    concentrated.nseriesInput(x)
    concentrated.nseriesOutput(x)
}
\arguments{
    \item{x}{A concentrated TSdata object.}
}
\value{An integer.}

\seealso{
    \code{\link{concentratedDimension}}
    \code{\link{concentrate}}
}
\examples{
    data("eg1.DSE.data", package="dse1")
    require("stats")
    z <- concentrate(eg1.DSE.data)
    concentrated.nseriesOutput(z)
}
\concept{DSE}
\keyword{ts}

