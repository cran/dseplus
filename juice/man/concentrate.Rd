\name{concentrate}
\alias{concentrate}
\alias{concentrate.default}
\alias{concentrate.TSdata}
\alias{is.concentrate}
\alias{is.TSdataconcentrate}
\alias{is.TSmodelconcentrate}

\title{Concentrate Series in a TSdata Object}
\description{
    Calculate a reduced dimesion version of the data using principal 
    components (or cannonical correlation for TSdata with input and output).
}
\usage{
   concentrate(d, conc=NULL, center=TRUE, scale=TRUE, ...)
   \method{concentrate}{default}(d, conc=NULL, center=TRUE, scale=TRUE, n=1, ...)
   \method{concentrate}{TSdata}(d,conc=NULL, center=TRUE, scale=TRUE, m=1, p=1, ...)
   is.concentrate(x)
   is.TSdataconcentrate(x)
   is.TSmodelconcentrate(x)
   }
\arguments{
    \item{d}{a matrix or TSdata object.}
    \item{...}{arguments to be passed to other methods.}
    \item{conc}{object containing the concentrator (projection) matrix 
       used for the reduction}
    \item{center}{center the observations to mean zero first (passed to estProjection).}
    \item{scale}{scale the observations to SD one first (passed to estProjection).}
    \item{n}{dimension of the concentrated series (passed to estProjection).}
    \item{m}{dimension of the concentrated input series (passed to estProjection).}
    \item{p}{dimension of the concentrated output series (passed to estProjection).}
    \item{x}{any object.}
}
\value{A matrix or TSdata object.}
\seealso{
    \code{\link{estProjection}}
    \code{\link{reconstitute}}
    \code{\link{prcomp}}
}
\examples{
    if(is.R()) {
       data("eg1.DSE.data", package="dse1")
       require("mva")
    }
    z <- concentrate(eg1.DSE.data)
    is.concentrate(z)
}
%\keyword{DSE}
\keyword{ts}

