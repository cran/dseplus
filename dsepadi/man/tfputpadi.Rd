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
\value{A tfPADIdata object.}

\details{
    This function writes data to 
    a database using the TS PADI programs available 
    at \url{www.bank-banque-canada.ca/pgilbert} . 
}

\seealso{ 
    \code{\link{tfPADIdata}}
    \code{\link[padi]{putpadi}}
}
\keyword{tframe}
\keyword{time series data}

