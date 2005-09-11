\name{concentrateOriginal}
\alias{concentrateOriginal}
\alias{concentrateOriginal.TSdataconcentrate}
\alias{concentrateOriginal.TSdatareconstitute}
\alias{concentrateOriginal.concentrate}
\alias{concentrateOriginal.TScanonicalPrediction}

\title{Extract Original Series}
\description{concentrateOriginal}
\usage{
    concentrateOriginal(d)
    \method{concentrateOriginal}{TSdataconcentrate}(d)
    \method{concentrateOriginal}{TSdatareconstitute}(d)
    \method{concentrateOriginal}{concentrate}(d)
    \method{concentrateOriginal}{TScanonicalPrediction}(d)
}
\arguments{
    \item{d}{A concentrate object.}
}
\value{A TSdata object.}
\details{
    The original data set is returned as a TSdata object, stripped of
    the fact that it is a concentrate.
}
\seealso{
    \code{\link{concentrate}}
    \code{\link{concentrator}}
    \code{\link{concentrateOnly}}
}

\examples{
    data("eg1.DSE.data", package="dse1")
    require("stats")
    z <- concentrate(eg1.DSE.data)
    z <- concentrateOriginal(z)
}
\concept{DSE}
\keyword{ts}

