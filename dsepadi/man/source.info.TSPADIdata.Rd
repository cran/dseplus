\name{source.info.TSPADIdata}
\alias{source.info.TSPADIdata}
\alias{sourcedb.TSPADIdata}
\alias{sourceserver.TSPADIdata}
\alias{availability.TSPADIdata}
\alias{identifiers.TSPADIdata}
\alias{periods.TSPADIdata}

\title{TSPADIdata Specific Methods}
\description{See the generic function description.}
\usage{
    \method{source.info}{TSPADIdata}(obj)
    \method{sourcedb}{TSPADIdata}(obj)
    \method{sourceserver}{TSPADIdata}(obj)
    \method{availability}{TSPADIdata}(obj, verbose=T, timeout=60)
    \method{identifiers}{TSPADIdata}(obj)
    \method{periods}{TSPADIdata}(x)
}
\arguments{
    \item{obj}{a TSPADIdata object.}
    \item{x}{a TSPADIdata object.}
    \item{verbose}{a logical indicating if extra information should be printed.}
    \item{timeout}{an integer indicating the number of seconds to wait before
       concluding that the server is not available.}
}\seealso{
    \code{\link{source.info}} 
}
%\keyword{DSE}
\keyword{ts}

