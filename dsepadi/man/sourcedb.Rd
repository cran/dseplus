\name{sourcedb}
\alias{sourcedb}
\alias{sourcedb.default}
\alias{sourcedb.tfPADIdata}
\title{Get sourcedb from a data object}
\description{Get sourcedb from a data object}
\usage{
    sourcedb(obj)
    \method{sourcedb}{default}(obj)
    \method{sourcedb}{tfPADIdata}(obj)
}
\arguments{
    \item{obj}{An object which contains source data base information.}
}
\value{
A list with elements input and output which are strings indicating the 
input and output sourcedb.
}
\seealso{
\code{\link[dse2]{TSPADIdata}} 
\code{\link{identifiers}}
\code{\link{source.info}}
}
\examples{
if(require("padi")) {
   d <- tfPADIdata("D1", server="ets")
   d <- freeze(d)
   sourcedb(d) }
}
%\keyword{tframe}
\keyword{ts}

