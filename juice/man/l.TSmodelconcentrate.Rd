\name{l.TSmodelconcentrate}
\alias{l.TSmodelconcentrate}
\title{Specific Methods for l}
\description{See the generic function description.}
\usage{
    \method{l}{TSmodelconcentrate}(model, data, sampleT=nrow(output.data(data)), 
                                  predictT=sampleT, result=NULL, warn=TRUE)
}
\arguments{
    \item{model}{a TSmodelconcentrate object.}
    \item{data}{a TSdataconcentrate object.}
    \item{sampleT}{an integer indicating the number of periods of data to use.}
    \item{predictT}{an integer to what period forecasts should be extrapolated.}
    \item{result}{
      if non-NULL then the returned value is only the sub-element indicated by 
      result. result can be a character string or integer.}
    \item{warn}{if FALSE then certain warning messages are turned off.}
}
\seealso{
   \code{\link[dse1]{l}} 
   \code{\link[dse1]{l.ARMA}} 
   \code{\link[dse1]{l.SS}} 
}

%\keyword{DSE}
\keyword{ts}

