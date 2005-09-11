\name{concentratedSeriesNames}
\alias{concentratedSeriesNames}
\alias{concentratedSeriesNames.concentrate}
\alias{concentratedSeriesNames.TSdata}
\alias{concentratedSeriesNamesInput}
\alias{concentratedSeriesNamesOutput}

\title{Concentrated Series Names}
\description{
    The names of series in concentrate data.
}
\usage{
    concentratedSeriesNames(x)
    \method{concentratedSeriesNames}{concentrate}(x)
    \method{concentratedSeriesNames}{TSdata}(x)
    concentratedSeriesNamesInput(x)
    concentratedSeriesNamesOutput(x)
}
\arguments{
    \item{x}{A concentrated data object.}
}
\value{
    Depends on the argument. For a simple concentrated data object 
    the result is a vector of strings.
}

\seealso{
    \code{\link[tframe]{seriesNames}}
    \code{\link[dse1]{seriesNamesInput}}
    \code{\link[dse1]{seriesNamesOutput}}
    \code{\link{concentratedDimension}}
    \code{\link{concentrate}}
}
\examples{
    data("eg1.DSE.data", package="dse1")
    require("stats")
    z <- concentrate(eg1.DSE.data)
    concentratedSeriesNames(z)
}
\concept{DSE}
\keyword{ts}

