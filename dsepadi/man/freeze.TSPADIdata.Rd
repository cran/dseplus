\name{freeze.TSPADIdata}
\alias{freeze.TSPADIdata}
\title{Get fixed data snapshot using TSPADI}
\description{
Get a TSdata structure from TSPADI Database Interface}
\usage{
    \method{freeze}{TSPADIdata}(data, timeout=60, ...)
}
\arguments{
    \item{data}{An object of class TSPADIdata.}
    \item{timeout}{an integer indicating the number of seconds to wait before
       concluding that the server is not available.}
    \item{...}{Arguments passed to other methods.}
}
\value{A TSdata object taken from the database.
}
\details{
This function extracts data using the TS PADI interface
to a database. See TSPADIdata, tfPADIdata, and freeze for more details. 
}
\seealso{
  
\code{\link[tframe]{freeze}}
\code{\link[dse1]{TSdata}}
\code{\link{TSPADIdata}}
}
\examples{
   inflation.sa.names <- TSPADIdata2(
    output = list(
         c("ets","", "b820678",  "ytoypc", "price services"),
         c("ets","", "b800938",  "ytoypc", "total unit labour costs")
           ), stop.on.error=TRUE, warn=TRUE )
   if(require("padi") & require("dse2") && checkPADIserver("ets"))
      z <- freeze(inflation.sa.names)
}
%\keyword{DSE}
\keyword{ts}

