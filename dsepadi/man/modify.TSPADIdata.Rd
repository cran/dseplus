\name{modify.TSPADIdata}
\alias{modify.TSPADIdata}
\title{Modify a TSPADIdata Object}
\description{Modify a TSPADIdata object with new information}
\usage{
    \method{modify}{TSPADIdata}(obj,
                        output=NA,             input=NA,
                        output.server=NA,      input.server=NA,
                        output.db=NA,          input.db=NA,
                        output.transforms=NA,  input.transforms=NA, 
                        output.names=NA,       input.names=NA,
                         start=NA, end=NA, frequency=NA, 
                         pad.start=NA, pad.end=NA,
                         server=NA, db=NA, start.server=NA, 
                         server.process=NA, cleanup.script=NA,
                         stop.on.error=NA, warn=NA)
  
}

\arguments{
    \item{see TSPADIdata}{As for TSPADIdata.}
}
\value{An object of class "TSPADIdata" "TSdata".}
\details{This is the way to modify a TSPADIdata object.}
\seealso{
\code{\link{TSPADIdata}}
\code{\link{modify.tfPADIdata}}
\code{\link{freeze.TSPADIdata}}
\code{\link{freeze.TSPADIdata}}
\code{\link{TSPADIdata}}
}
\examples{
inflation.sa.names <- TSPADIdata2(
    output = list(
         c("ets","", "b820678",  "ytoypc", "price services"),
         c("ets","", "b800938",  "ytoypc", "total unit labour costs")
           ), stop.on.error=T, warn=T )
z <- modify.TSPADIdata(inflation.sa.names, server="newserver")
}
%\keyword{DSE}
\keyword{ts}

