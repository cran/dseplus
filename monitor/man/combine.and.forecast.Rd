\name{combine.and.forecast}
\alias{combine.and.forecast}
\alias{reconstruct.combined.forecast}
\alias{construct.data.to.override.horizon}
\alias{get.overriding.data}
\alias{restrict.overriding.data}
 
\title{Simple Monitoring Utility Functions}
\description{Utility functions for simple monitoring.}
\usage{
combine.and.forecast(model, new.data,  
      overlapping.period.forecast.tag="g", forecast.tag="f") 
reconstruct.combined.forecast(combined.forecast) 
construct.data.to.override.horizon(new.data, model, plot=TRUE, forecast.tag="f")
get.overriding.data(file="overriding.data", 
    first.input="",first.output="", second.output="", m=1, p=10)
restrict.overriding.data(data, overriding.horizon=0)  
}
\value{depends}
\details{Internal functions, not to be called by user.}
\seealso{
    \code{\link{combination.monitoring}}
}

%\keyword{DSE}
\keyword{internal}
