\name{modify.TSPADIdata}
\alias{modify.TSPADIdata}
\title{Modify a TSPADIdata Object}
\description{Modify a TSPADIdata object with new information}
\usage{
    \method{modify}{TSPADIdata}(obj,
                         start=NA, end=NA, frequency=NA, 
                         pad=NA, pad.start=pad, pad.end=pad,
                         server=NA, db=NA, start.server=NA, 
                         server.process=NA, cleanup.script=NA,
                         stop.on.error=NA, warn=NA,
			 append=NA, use.tframe=NA,
                        output=NA,             input=NA,
                        output.server=NA,      input.server=NA,
                        output.db=NA,          input.db=NA,
                        output.transforms=NA,  input.transforms=NA, 
                        output.names=NA,       input.names=NA,
			...)  
}

\arguments{
    \item{obj}{a tfPADIdata object.}
    \item{output}{see tfPADIdata description of \code{series}.}
    \item{output.server}{see tfPADIdata descriptions.}
    \item{output.db}{see tfPADIdata descriptions.}
    \item{output.transforms}{see tfPADIdata descriptions.}
    \item{output.names}{see tfPADIdata descriptions.}
    \item{input}{see tfPADIdata descriptions of \code{series}}
    \item{input.server}{see tfPADIdata descriptions.}
    \item{input.db}{see tfPADIdata descriptions.}
    \item{input.transforms}{see tfPADIdata descriptions.}
    \item{input.names}{see tfPADIdata descriptions.}

    \item{start}{see tfPADIdata descriptions.}
    \item{end}{see tfPADIdata descriptions.}
    \item{frequency}{see tfPADIdata descriptions.}
    \item{pad}{see tfPADIdata descriptions.}
    \item{pad.start}{see tfPADIdata descriptions.}
    \item{pad.end}{see tfPADIdata descriptions.}
    \item{start.server}{see tfPADIdata descriptions.}
    \item{server.process}{see tfPADIdata descriptions.}
    \item{server}{see tfPADIdata descriptions.}
    \item{db}{see tfPADIdata descriptions.}
    \item{cleanup.script}{see tfPADIdata descriptions.}
    \item{stop.on.error}{see tfPADIdata descriptions.}
    \item{warn}{see tfPADIdata descriptions.}
    \item{append}{see tfPADIdata descriptions.}
    \item{use.tframe}{see tfPADIdata descriptions.}
    \item{...}{arguments to be passed to other methods.}
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
           ), stop.on.error=TRUE, warn=TRUE )
z <- modify(inflation.sa.names, server="newserver")
}
%\keyword{DSE}
\keyword{ts}

