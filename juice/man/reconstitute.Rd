\name{reconstitute}
\alias{reconstitute}
\alias{reconstitute.default}
\alias{reconstitute.concentrate}
\alias{reconstitute.TSdataconcentrate}
\alias{is.TSdatareconstitute}

\title{Reconstitute}
\description{reconstitute}
\usage{
    reconstitute(d, conc=NULL, names=NULL)
    \method{reconstitute}{default}(d, conc=NULL, names=seriesNames(d))
    \method{reconstitute}{concentrate}(d, conc=concentrator(d),
                                    names=seriesNames(d))
    \method{reconstitute}{TSdataconcentrate}(d, conc=concentrator(d),
                                    names=seriesNames(d))
    is.TSdatareconstitute(x)
}
\arguments{
    \item{d}{a concentrated data object.}
    \item{conc}{a concentrator.}
    \item{names}{series names for the result.}
    \item{x}{any object.}
}

\value{Depends on the argument.}
\details{
    A concentrated data object is used to reconstruct the full dimension 
    data. Thus the result has the same dimension as the original data, but
    will not be exactly the same because some information is lost when the
    data is concentrated (unless the concentrate has the full
    dimension of the original data, which would usually be pointless). 
}
\seealso{
    \code{\link{concentrate}}
}
\examples{
    data("eg1.DSE.data", package="dse1")
    require("stats")
    z <- concentrate(eg1.DSE.data)
    z <- reconstitute(z)
    is.TSdatareconstitute(z)
}
\concept{DSE}
\keyword{ts}

