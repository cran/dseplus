\name{putpadi}
\alias{putpadi}
\alias{putpadi.default}

\title{PADI Functions}
\description{see details}
\usage{
  putpadi(data,  server=Sys.info()[["nodename"]], dbname="",
                   series=seriesNames(data),
                   start.server=TRUE, server.process=PADIserverProcess(),
                   cleanup.script=PADIcleanupScript(),
                   user=Sys.info()[["user"]], passwd="",
                   stop.on.error=TRUE, warn=TRUE, timeout=60)
  \method{putpadi}{default}(data,  server=Sys.info()[["nodename"]], dbname="",
                   series=seriesNames(data),
                   start.server=TRUE, server.process=PADIserverProcess(),
                   cleanup.script=PADIcleanupScript(),
                   user=Sys.info()[["user"]], passwd="",
                   stop.on.error=TRUE, warn=TRUE, timeout=60)
}
\arguments{
    \item{data}{a tfPADIdata object with data to put on the server.}
    \item{server}{string or vector of strings indicating the server where
       the series will be found. Scalar values are expanded to a 
       vector of appropriate length.}
    \item{dbname}{string or vector of strings indicating additional information for 
       the server. Scalar values are expanded to a vector of appropriate length.}
    \item{series}{vector of strings to use for series identifiers.}
    \item{start.server}{logical indicating if a (local) server should be
       started.}
    \item{server.process}{string indicating the name to be used to start a
       server process.}
    \item{cleanup.script}{string indicating the name to be used to shut down a
       server process.}
    \item{user}{an optional string used by the server to check permission.}
    \item{passwd}{an optional string used by the server to check permission.}
    \item{stop.on.error}{logical indicating if the function should stop if any
       series produces an error, or continue with other series.}
    \item{warn}{logical indicating if warning messages should be supressed.}
    \item{timeout}{an integer indicating the number of seconds to wait before
       concluding that the server is not available.}
}
\value{depends}
\details{
 The function putpadi writes data to a specified databases on a specified 
      server (default local).

 If start.server is FALSE then there will be no attempt to start a server and the 
      function will stop if a server is not running.
      If start.server is TRUE (the default) then if server==Sys.info()[["nodename"]] (the
      default) and a server is not already running, there will be an attempt
      to start a server using the argument server.process. The default is
      a  "local mode" Fame server.

 The argument server.process is only used if it is necessary to start a server.

 If specified, server would typically be a single string, though there is
    some attempt to handle vectors of strings (indicating different servers
    for each series).

 dbname should typically be supplied. (but some servers might accept an empty
    string ("") or ignore any string supplied.
    If the server is being started then it will run in the Unix pwd and local
    path names should work, BUT in general there is no guarantee that the 
    server is running in the pwd and complete path names may be required.
 data can be a vector, matrix, time series or time series matrix.
    Dates are used when available. If data is a vector or single time 
 series, then server, dbname,  and series should have length 1. 
 If data is a matrix or time series matrix with more than one series
 then server, dbname,  and series should be character strings
 of the same length or, server and dbname can be length 1, in
 which case the string will be repeated for all series.
 If they do not exist then an error is indicated.
 The strings user and passwd are passed to the server, but may or may not be
   used, depending on the server implementation.
 If stop.on.error is TRUE (default) the function stops if there is an error
   writing any series. Otherwise, the result is TRUE or FALSE for each series, 
   indicating success or failure.
}
\seealso{
   \code{\link{checkPADIserver}}
   \code{\link{getpadi}}
}
%\examples{
%}
\concept{DSE}
\keyword{ts}


