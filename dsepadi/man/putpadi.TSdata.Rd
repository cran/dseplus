\name{putpadi.TSdata}
\alias{putpadi.TSdata}
\title{Specific Methods for putpadi}
\description{See the generic function description.}
\usage{
    \method{putpadi}{TSdata}(data, dbname, server=Sys.info()[["nodename"]], 
                   start.server=T, server.process=padi.server.process(), 
                   cleanup.script=padi.cleanup.script(),
                   series=seriesNames(data),
                   user=Sys.info()[["user"]], passwd="",
                   stop.on.error=T, warn=T)
}
\arguments{
    \item{data}{a TSdata object.}
    \item{dbname}{a string used by the server to indicate where the 
       data should go.}
    \item{server}{see tfPADIdata descriptions.}
    \item{start.server}{see tfPADIdata descriptions.}
    \item{server.process}{see tfPADIdata descriptions.}
    \item{cleanup.script}{see tfPADIdata descriptions.}
    \item{series}{see tfPADIdata descriptions.}
    \item{user}{an optional string used by the server to check permission.}
    \item{passwd}{an optional string used by the server to check permission.}
    \item{stop.on.error}{see tfPADIdata descriptions.}
    \item{warn}{see tfPADIdata descriptions.}
}
\seealso{
    \code{\link[tframe]{tfputpadi}}
    \code{\link[padi]{putpadi}}
}

%\keyword{DSE}
\keyword{ts}


