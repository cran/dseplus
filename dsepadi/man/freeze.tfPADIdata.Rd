\name{freeze.tfPADIdata}
\alias{freeze.tfPADIdata}
\alias{freeze.FAMEdata}
\title{Get fixed data snapshot}
\description{Get a time series matrix structure from a database}
\usage{
    \method{freeze}{tfPADIdata}(data, timeout=60, ...)
    \method{freeze}{FAMEdata}(data, ...)
}
\arguments{
    \item{data}{A description of how to get data.}
    \item{timeout}{an integer indicating the number of seconds to wait before
       concluding that the server is not available.}
    \item{...}{Arguments passed to other methods.}
}
\value{A time series matrix.}
\details{
This function extracts data from 
a database (for example using the TS PADI programs are available 
at \url{www.bank-banque-canada.ca/pgilbert} ). This method is generic.
Typically the argument data is a tfPADIdata or TSPADIdata object identifying 
the source of 
the data. See help for tfPADIdata and TSPADIdata. 
The default method usually just
returns its argument, so freeze has no effect. This way freeze
can be used to write functions which will take a snapshot from the database
when they execute or will work with an already fixed copy of data if that is
what is supplied. The default does allow for a character argument, in which case
it is used to construct a tfPADIdata object using server="ets", then freeze that
object. This allows for a simple syntax to grab a series from the database. The
server="ets" is for convenience at the Bank of Canada and another default server
might be more convenient elsewhere.
}

\seealso{ 
\code{\link{tfPADIdata}}
\code{\link{freeze.TSPADIdata}}
\code{\link{TSPADIdata}}
}
\examples{

if (require("padi") && checkPADIserver("ets")){
 z1 <- freeze("D1")
 z2 <- freeze(tfPADIdata(c("B1630", "B1642"), server="ets"))
  }
}
%\keyword{tframe}
\keyword{ts}

