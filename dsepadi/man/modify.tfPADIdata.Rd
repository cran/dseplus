\name{modify.tfPADIdata}
\alias{modify.tfPADIdata}
\title{modify a database descriptor object}
\description{Modify an object with new information}
\usage{
    \method{modify}{tfPADIdata}(obj, append=NA, 
           series=NA, server=NA, db=NA, transforms=NA, 
           start=NA, end=NA, frequency=NA, names=NA, 
           pad=NA, pad.start=NA, pad.end=NA,
           use.tframe=NA,
           start.server=NA, server.process=NA, cleanup.script=NA,
           stop.on.error=NA, warn=NA)
}
\arguments{
    \item{obj}{a tfPADIdata object.}
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
\code{\link{freeze}}
}
\examples{
  z <- tfPADIdata( c("seriesA", "seriesB"), server="myserver")}
  z <- modify(z, pad.end=T )
 }
%\keyword{tframe}
\keyword{ts}

