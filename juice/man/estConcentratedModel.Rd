\name{estConcentratedModel}
\alias{estConcentratedModel}
\alias{estConcentratedModel.TSdata}
\alias{estConcentratedModel.TSdataconcentrate}

\title{Estimate a Concentrated Model}
\description{estConcentratedModel}
\usage{
    estConcentratedModel(data, estimation="estVARXls",
                           estimation.args=NULL, ...)
    \method{estConcentratedModel}{TSdata}(data, estimation="estVARXls",  
               estimation.args=NULL, m=1, p=1, center=TRUE, scale=TRUE, ...)
    \method{estConcentratedModel}{TSdataconcentrate}(data, 
                estimation="estVARXls", estimation.args=NULL, warn=TRUE, ...)
}
\arguments{
    \item{data}{A TSdata or TSdataconcentrate object.}
    \item{estimation}{Estimation method.}
    \item{estimation.args}{Estimation method arguments.}
    \item{m, p}{dimension of the concentrated series.}
    \item{center}{center the observations to mean zero first.}
    \item{scale}{scale the observations to SD one first.}
    \item{warn}{logical indicating if certain warning messsages should be
      printed.}
    \item{...}{arguments to be passed to other methods.}
}
\value{A TSmodelconcentrate.}
\details{
    A concentrated version of the data (reduced dimension) is used to 
    estimate a reduced dimesion model. The projections for concentrating the
    data are retained so that model predictions can be expanded to the full
    dimension data space.
    
    If data is TSdataconcentrate then the concentrator with that data is 
    used and m, p, center and scale are not used. For TSdata these arguments
    are used to first estimate a concentrated version of the data.
}
\seealso{
    \code{\link{concentrate}}
    \code{\link{concentrator}}
    \code{\link{estProjection}}
}

\examples{
    data("eg1.DSE.data.diff", package="dse1")
    model <- estConcentratedModel(eg1.DSE.data.diff)
}
\concept{DSE}
\keyword{ts}

