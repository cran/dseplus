\name{tfputpadi}
\alias{tfputpadi}

\title{Write Data to a Data Server}
\description{Write data to a server.}
\usage{
    tfputpadi(data,  
         server = Sys.info()[["nodename"]],
         dbname = "", 
         series = series.names(data),
         start.server = T,
         server.process = padi.server.process(), 
         cleanup.script = padi.cleanup.script(),
         user = Sys.info()[["user"]], passwd= "",
         stop.on.error = T, warn = T)
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
}

\value{A tfPADIdata object (as would be suitable for retrieving the data).}

\details{
    This function writes data to 
    a database using the TS PADI programs available 
    at \url{www.bank-banque-canada.ca/pgilbert} . 
}

\seealso{ 
    \code{\link{tfPADIdata}}
    \code{\link[padi]{putpadi}}
}
%\keyword{tframe}
\keyword{ts}


