\name{print.tfPADIdata}
\alias{print.tfPADIdata}

\title{Print a tfPADIdata Object}
\description{See the generic function description.}

\usage{
    \method{print}{tfPADIdata}(x, ...)
}
\arguments{
    \item{x}{a tfPADIdata object.}
    \item{...}{arguments to be passed to other methods.}
}
\section{Side Effects}{The object is printed.}
\seealso{
    \code{\link[base]{print}}
    \code{\link{tfPADIdata}}
}
\examples{
if(require("padi") && checkPADIserver("ets")) 
    print(tfPADIdata(c("series1", "series2"), server="ets"))
}
%\keyword{tframe}
\keyword{ts}

