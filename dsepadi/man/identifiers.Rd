\name{identifiers}
\alias{identifiers}
\alias{identifiers.default}
\alias{identifiers.tfPADIdata}
\title{Get Identifiers from a data object}
\description{Get Identifiers from a data object}
\usage{
    identifiers(obj)
    \method{identifiers}{default}(obj)
    \method{identifiers}{tfPADIdata}(obj)
}
\arguments{
    \item{obj}{An object which contains source series identifier information.}
}
\value{
A list with elements input and output which are strings indicating the 
input and output series identifiers.
}
\seealso{ 
\code{\link{tfPADIdata}}
\code{\link[dse2]{TSPADIdata}}
\code{\link{sourceInfo}}
\code{\link{sourcedb}}
}
\examples{
 if(require("padi") && checkPADIserver("ets")) {
  d <- tfPADIdata("D1", server="ets")
  identifiers(d)
  }
 }
%\keyword{tframe}
\keyword{ts}

