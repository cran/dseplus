\name{est.projection}
\alias{est.projection}
\alias{est.projection.default}
\alias{est.projection.TSdata}
\title{Calculate Projection from Concentrating Series.}
\description{
Calculate the projection to use for a reduced dimesion version of the data 
using principal components (or
cannonical correlation for TSdata with input and output).
}
\usage{
   est.projection(data, center=TRUE, scale=TRUE, ...)     
   \method{est.projection}{default}(data, center=TRUE, scale=TRUE, n=1, ...) 
   \method{est.projection}{TSdata}(data, center=TRUE, scale=TRUE, m=1,p=1, ...) 
   }
\arguments{
    \item{data}{a matrix for the default method or TSdata object.}
    \item{n}{dimension of the concentrated series.}
    \item{m}{dimension of the concentrated input series.}
    \item{p}{dimension of the concentrated output series.}
    \item{center}{logical indicating center the observations to mean zero first.}
    \item{scale}{logical indicating scale the observations to SD one first.}
    \item{...}{arguments passed to other methods.}
}
\value{An object containing matrix (conc) to use to 
 concentrate the data.}
\seealso{
\code{\link{est.concentratedModel}}
\code{\link{concentrate}}
\code{\link{reconstitute}}
\code{\link{prcomp}}
}
\examples{
if(is.R()) {
   data("eg1.DSE.data", package="dse1")
   require("mva")
  }
z <- est.projection(eg1.DSE.data)
}
%\keyword{DSE}
\keyword{ts}

