\name{availability}
\alias{availability}
\alias{availability.default}
\alias{availability.tfPADIdata}
\title{Check Data Availability}
\description{Check the dates for which date is available.}
\usage{
    availability(obj, ...)
    \method{availability}{default}(obj, names=NULL, server="ets", dbname="",
                       verbose=TRUE, timeout=60, stop.on.error=TRUE, warn=TRUE)
    \method{availability}{tfPADIdata}(obj, verbose=TRUE, timeout=60)
}
\arguments{
    \item{obj}{a vector of series identifiers or an object of 
       class \code{tfPADIdata}. See TSdata.TSPADIdata.}
    \item{...}{arguments to be passed to other methods.}
    \item{verbose}{TRUE or FALSE indicating if the results should be printed.}
    \item{timeout}{an integer indicating the number of seconds to wait before
       concluding that the server is not available.}
    \item{names}{
    A character vector of names to be associated with the ids in printed results.}
    \item{server}{
    A character vector of servers to be associated with the ids.}
    \item{dbname}{
    A character vector of databases to be associated with the ids.}
    \item{stop.on.error}{logical indicating if the function should stop if any
       series produces an error, or continue with other series.}
    \item{warn}{logical indicating if warning messages should be supressed.}
}
\value{A list with elements start, end, and frequency.}
\details{
If verbose is TRUE then the start, end, and frequency are printed for each 
series in data.id. A list is return (invisibly) with the same information. 
The default method works for a character vector argument. An argument of 
class TSPADIdata may supply some additional information, and several of 
the indicated optional arguments do not apply to that method.
}
\seealso{
\code{\link{tfPADIdata}},
\code{\link{TSPADIdata}},
\code{\link{set.TSPADIdata}},
\code{\link{retrieve.and.verify.data}}
}
\examples{
if(require("padi")) {
   library()  # debugging to see where padi is found if this fails
   d <- tfPADIdata("D1", server="ets")
   availability(d)}
}
%\keyword{tframe}
\keyword{ts}

