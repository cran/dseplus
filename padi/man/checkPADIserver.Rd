\name{checkPADIserver}
\alias{checkPADIserver}
\alias{load.padi}

\title{PADI Functions}
\description{see details}
\usage{
  load.padi(from = Sys.getenv("PADI_LDLIB"))
  checkPADIserver(server=Sys.info()[["nodename"]],
        user=Sys.info()[["user"]], timeout=50)
}
\arguments{
    \item{from}{string indicating the location of compiled PADI code to be
         loaded in S. This argument is ignored in R."}
    \item{server}{string indicating the network name of the server."}
    \item{user}{optional string indicating user ID used by the server to 
       start the process."}
    \item{timeout}{an integer indicating the number of seconds to wait before
       concluding that the server is not available.}
}
\value{depends}
\details{
 The function load.padi loads the associated compiled C object code.
 The function checkPADIserver checks if a server is running.
}
\seealso{
   \code{\link{getpadi}}
   \code{\link{putpadi}}
}
%\examples{
%}
\concept{DSE}
\keyword{ts}

