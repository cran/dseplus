\name{checkForValueChanges}
\alias{checkForValueChanges}
\alias{checkForFileDateChanges}
\alias{watch.data}
 
\title{Simple Monitoring Utility Functions}
\description{Utility functions for simple monitoring.}
\usage{
checkForValueChanges(data.names, verification.data,
     discard.current=FALSE, ignore.before= NULL, fuzz=1e-10)
checkForFileDateChanges(data.names, verification.dates)
watch.data(data.names, previous.data=NULL, mail.list="gilp", error.mail.list=NULL,
   message.title="Data Monitor\n", message.subject="Data Monitor", message.footnote=NULL)
}
\value{depends}
\details{Internal functions, not to be called by user.}
\seealso{
    \code{\link{simpleMonitoring}}
}

%\keyword{DSE}
\keyword{internal}
