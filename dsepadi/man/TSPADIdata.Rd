\name{TSPADIdata}
\alias{TSPADIdata}
\alias{is.TSPADIdata}

\title{Constructor for TSPADIdata object}
\description{Make a TSPADIdata object from identifier information}
\usage{
    TSPADIdata(output=NULL,           input=NULL,
                        output.server=server,  input.server=server,
                        output.db=db,          input.db=db,
                        output.transforms="",  input.transforms="", 
                        output.names=NULL,     input.names=NULL,
                         start=NA, end=NA, frequency=NA, 
                         pad=FALSE, pad.start=pad, pad.end=pad,
                         server="", db="", start.server=NULL, 
                         server.process=NULL, cleanup.script=NULL,
                         stop.on.error=TRUE, warn=TRUE)
  
    is.TSPADIdata(obj)
}
\arguments{
    \item{output}{see tfPADIdata description of \code{series}.}
    \item{output.server}{see tfPADIdata descriptions.}
    \item{output.db}{see tfPADIdata descriptions.}
    \item{output.transforms}{see tfPADIdata descriptions.}
    \item{output.names}{see tfPADIdata descriptions.}
    \item{input}{see tfPADIdata description of \code{series}.}
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
    \item{server}{see tfPADIdata descriptions.}
    \item{db}{see tfPADIdata descriptions.}
    \item{start.server}{see tfPADIdata descriptions.}
    \item{server.process}{see tfPADIdata descriptions.}
    \item{cleanup.script}{see tfPADIdata descriptions.}
    \item{stop.on.error}{see tfPADIdata descriptions.}
    \item{warn}{see tfPADIdata descriptions.}
    \item{obj}{any object.}
}

\value{An object of class "TSPADIdata" "TSdata".}
\details{
This is the way to produce an object which can be used for PADI database access.
Input and output values are passed to tfPADIdata so the TSdata object return is 
a list of tfPADIdata objects called input and output (but either of these can 
be NULL).
}
\seealso{
\code{\link{freeze.TSPADIdata}}
\code{\link{freeze.TSPADIdata}}
\code{\link{tfPADIdata}}
\code{\link{TSPADIdata2}}
\code{\link{setTSPADIdata}}
}
\examples{
    z <- TSPADIdata(output=c("seriesA", "seriesB"), server="myserver")
    is.TSPADIdata(z)
}
%\keyword{DSE}
\keyword{ts}

