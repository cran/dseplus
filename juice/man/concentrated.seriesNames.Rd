\name{concentrated.seriesNames}
\alias{concentrated.seriesNames}
\alias{concentrated.seriesNames.concentrate}
\alias{concentrated.seriesNames.TSdata}
\alias{concentrated.seriesNamesInput}
\alias{concentrated.seriesNamesOutput}

\title{Concentrated Series Names}
\description{
    The names of series in concentrate data.
}
\usage{
    concentrated.seriesNames(x)
    concentrated.seriesNamesInput(x)
    concentrated.seriesNamesOutput(x)
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
    \code{\link{concentrated.dimension}}
    \code{\link{concentrate}}
}
\examples{
    if(is.R()) {
    data("eg1.DSE.data", package="dse1")
    require("mva")
    }
    z <- concentrate(eg1.DSE.data)
    concentrated.seriesNames(z)
}
%\keyword{DSE}
\keyword{ts}

