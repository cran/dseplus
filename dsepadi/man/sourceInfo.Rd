\name{sourceInfo}
\alias{sourceInfo}
\alias{sourceInfo.default}
\alias{sourceInfo.tfPADIdata}
\title{Get source information from a data object}
\description{Get source information from an object}
\usage{
    sourceInfo(obj)
    \method{sourceInfo}{default}(obj)
    \method{sourceInfo}{tfPADIdata}(obj)
}
\arguments{
    \item{obj}{An object which contains source series information.}
}
\value{
A list with elements input and output which are strings indicating the 
input and output source.
}
\seealso{
\code{\link[dse2]{TSPADIdata}} 
\code{\link[dse2]{identifiers}}
}
\examples{
if(require("padi") && checkPADIserver("ets")) {
   d <- tfPADIdata("D1", server="ets")
   d <- freeze(d)
   sourceInfo(d)}
}
%\keyword{tframe}
\keyword{ts}

