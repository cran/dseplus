\name{concentrated.dimension}
\alias{concentrated.dimension}
\alias{concentrated.dimension.concentrate}

\title{Concentrated Dimension}
\description{
    The dimension (number of series) in concentrate data. This is the
    dimension onto which the original series has been projected.
}
\usage{
    concentrated.dimension(x)
    \method{concentrated.dimension}{concentrate}(x)
    concentrated.nseriesInput(x)
    concentrated.nseriesOutput(x)
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
    if(is.R()) {
    data("eg1.DSE.data", package="dse1")
    require("mva")
    }
    z <- concentrate(eg1.DSE.data)
    concentrated.dimension(output.data(z))
    concentrated.nseriesOutput(z)
}
%\keyword{DSE}
\keyword{ts}

