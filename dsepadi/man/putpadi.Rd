\name{putpadi}
\alias{putpadi}
\alias{putpadi.TSdata}

\title{PADI Functions}
\description{see details}
\usage{
  putpadi(data, server=Sys.info()[["nodename"]],  dbname="",
                   series=seriesNames(data),
                   start.server=TRUE, server.process=PADIserverProcess(),
                   cleanup.script=PADIcleanupScript(),
                   user=Sys.info()[["user"]], passwd="",
                   stop.on.error=TRUE, warn=TRUE, timeout=60)
    \method{putpadi}{TSdata}(data, server=Sys.info()[["nodename"]], dbname, 
                   series=seriesNames(data),
                   start.server=TRUE, server.process=padi.server.process(), 
                   cleanup.script=padi.cleanup.script(),
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
 See putpadi.default.
}
\seealso{
%%   \code{\link[padi]{checkPADIserver}}
   \code{\link{getpadi}}
%%   \code{\link[padi]{putpadi.default}}
   \code{\link{putpadi.TSdata}}
}

\concept{DSE}
\keyword{ts}

