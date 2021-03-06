\name{concentrateOnly}
\alias{concentrateOnly}
\alias{concentrateOnly.concentrate}
\alias{concentrateOnly.TSdataconcentrate}
\alias{concentrateOnly.TSdatareconstitute}
\alias{concentrateOnly.TSestModel}
\alias{concentrateOnly.TSmodelconcentrate}

\title{Extract Concentrate}
\description{concentrateOnly}
\usage{
    concentrateOnly(d)
    \method{concentrateOnly}{concentrate}(d)
    \method{concentrateOnly}{TSdataconcentrate}(d)
    \method{concentrateOnly}{TSdatareconstitute}(d)
    \method{concentrateOnly}{TSestModel}(d)
    \method{concentrateOnly}{TSmodelconcentrate}(d)
}
\arguments{
    \item{d}{a concentrate object.}
}
\value{A TSdata object.}
\details{
    The concentrated data set is returned as a TSdata object, stripped of
    the fact that it is a concentrate.
}
\seealso{
    \code{\link{concentrate}}
    \code{\link{concentrator}}
    \code{\link{concentrateOriginal}}
}

\examples{
    data("eg1.DSE.data", package="dse1")
    require("stats")
    z <- concentrate(eg1.DSE.data)
    z <- concentrateOnly(z)
}
\concept{DSE}
\keyword{ts}

