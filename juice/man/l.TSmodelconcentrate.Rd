\name{l.TSmodelconcentrate}
\alias{l.TSmodelconcentrate}
\title{Specific Methods for l}
\description{See the generic function description.}
\usage{
    \method{l}{TSmodelconcentrate}(obj1, obj2, sampleT=nrow(outputData(obj2)), 
                                  predictT=sampleT, result=NULL, warn=TRUE, ...)
}
\arguments{
    \item{obj1}{a TSmodelconcentrate model object.}
    \item{obj2}{a TSdataconcentrate data object.}
    \item{sampleT}{an integer indicating the number of periods of data to use.}
    \item{predictT}{an integer to what period forecasts should be extrapolated.}
    \item{result}{
      if non-NULL then the returned value is only the sub-element indicated by 
      result. result can be a character string or integer.}
    \item{warn}{if FALSE then certain warning messages are turned off.}
    \item{...}{arguments passed to other methods.}
}
\seealso{
   \code{\link[dse1]{l}} 
   \code{\link[dse1]{l.ARMA}} 
   \code{\link[dse1]{l.SS}} 
}

%\keyword{DSE}
\keyword{ts}

