\name{modify}
\alias{modify}
\alias{modify.tfPADIdata}
\title{modify a database descriptor object}
\description{Modify a database descriptor object with new information}
\usage{
    modify(obj, start=NA, end=NA, frequency=NA, 
                         pad=NA, pad.start=pad, pad.end=pad,
                         server=NA, db=NA, start.server=NA, 
                         server.process=NA, cleanup.script=NA,
                         stop.on.error=NA, warn=NA,
			 append=NA, use.tframe=NA, ...)
    \method{modify}{tfPADIdata}(obj,
                         start=NA, end=NA, frequency=NA, 
                         pad=NA, pad.start=pad, pad.end=pad,
                         server=NA, db=NA, start.server=NA, 
                         server.process=NA, cleanup.script=NA,
                         stop.on.error=NA, warn=NA,
			 append=NA, use.tframe=NA,
                      series=NA, transforms=NA, names=NA, ...)
}
\arguments{
    \item{obj}{database descriptor object to be modified.}
    \item{append}{if specified, append should be a list of (series, 
       server, db, transforms, names) as specified for tfPADIdata. These 
       are appended as additional series to the object.}
    \item{series}{see tfPADIdata descriptions.}
    \item{server}{see tfPADIdata descriptions.}
    \item{db}{see tfPADIdata descriptions.}
    \item{transforms}{see tfPADIdata descriptions.}
    \item{start}{see tfPADIdata descriptions.}
    \item{end}{see tfPADIdata descriptions.}
    \item{frequency}{see tfPADIdata descriptions.}
    \item{names}{see tfPADIdata descriptions.}
    \item{pad}{see tfPADIdata descriptions.}
    \item{pad.start}{see tfPADIdata descriptions.}
    \item{pad.end}{see tfPADIdata descriptions.}
    \item{use.tframe}{see tfPADIdata descriptions.}
    \item{start.server}{see tfPADIdata descriptions.}
    \item{server.process}{see tfPADIdata descriptions.}
    \item{cleanup.script}{see tfPADIdata descriptions.}
    \item{stop.on.error}{see tfPADIdata descriptions.}
    \item{warn}{see tfPADIdata descriptions.}
    \item{...}{arguments to be passed to other methods.}
}
\value{
   A object of class "tfPADIdata" which can be used to retrieve a 
   matrix time series object.
}
\details{Any specified (non NA) optional arguments are used to modify 
the object. Values replace existing values (except in the case of append).
}
\seealso{
\code{\link{tfPADIdata}} 
\code{\link[tframe]{freeze}}
\code{\link[dse2]{modify.TSPADIdata}} 
}
\examples{
if(require("padi")) {
  z <- tfPADIdata( c("seriesA", "seriesB"), server="myserver")
  z <- modify(z, pad.end=TRUE )
  }
 }
%\keyword{tframe}
\keyword{ts}

