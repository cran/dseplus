\name{concentrator}
\alias{concentrator}
\alias{concentrator.concentrate}
\alias{concentrator.concentrator}
\alias{concentrator.TSdata}
\alias{concentrator.TSdataconcentrator}
\alias{concentrator.TSmodelconcentrate}
\alias{is.concentrator}
\alias{is.TSdataconcentrator}

\title{Data Concentrator}
\description{concentrator}
\usage{
    concentrator(d)
    \method{concentrator}{concentrate}(d)
    \method{concentrator}{concentrator}(d)
    \method{concentrator}{TSdata}(d)
    \method{concentrator}{TSdataconcentrator}(d)
    \method{concentrator}{TSmodelconcentrate}(d)
    is.concentrator(x)
    is.TSdataconcentrator(x)
}
\arguments{
    \item{d}{a concentrate or concentrator object.}
    \item{x}{any object.}
}
\value{A concentrator.}
\details{
    The concentrator is extracted from a concentrated data object.
}
\seealso{
    \code{\link{concentrate}}
    \code{\link{concentrateOnly}}
    \code{\link{concentrateOriginal}}
}

\examples{
    data("eg1.DSE.data", package="dse1")
    require("stats")
    z <- concentrate(eg1.DSE.data)
    concentrator(z)
    is.concentrator(concentrator(z))
}
\concept{DSE}
\keyword{ts}

