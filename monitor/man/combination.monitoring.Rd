\name{combination.monitoring}
\alias{combination.monitoring}
\alias{combine.and.forecast}
\alias{reconstruct.combined.forecast}
\alias{construct.data.to.override.horizon}
\alias{get.overriding.data}
\alias{restrict.overriding.data}

\title{Combination Monitoring}
\description{Automatic monitoring with e-mail of results}
\usage{
    combination.monitoring(model, data.names,
       previous.data=NULL,
       overriding.data.names=NULL, 
       restrict.overriding.data=T, overriding.horizon=0,
       mail.list=NULL,
       error.mail.list=NULL,
       message.title="Combination Monitoring",
       message.subject="Combination Monitoring",
       message.footnote=NULL,
       show.start= c(0,-3),
       show.end  = c(0,12),    
       report.variables=series.names(data.names),
       data.sub.heading=NULL,
       data.tag=" ",
       future.input.data.tag="p",
       overriding.data.tag="m",
       overlapping.period.forecast.tag="g",
       forecast.tag="f",
       run.again=F,
       save.as=NULL)
}
\arguments{
    \item{model}{see \code{simple.monitoring}.}
    \item{data.names}{see \code{simple.monitoring}.}
    \item{previous.data}{see \code{simple.monitoring}.}
    
    \item{mail.list}{see \code{simple.monitoring}.}
    \item{error.mail.list}{see \code{simple.monitoring}.}
    \item{message.title}{see \code{simple.monitoring}.}
    \item{message.subject}{see \code{simple.monitoring}.}
    \item{message.footnote}{see \code{simple.monitoring}.}
    \item{show.end}{see \code{simple.monitoring}.}
    \item{show.start}{see \code{simple.monitoring}.}
    \item{report.variables}{see \code{simple.monitoring}.}
    \item{data.sub.heading}{see \code{simple.monitoring}.}
    \item{data.tag}{see \code{simple.monitoring}.}

    \item{future.input.data.tag}{.}
    \item{overriding.data.tag}{.}
    \item{overlapping.period.forecast.tag}{.}

    \item{forecast.tag}{see \code{simple.monitoring}.}
    \item{run.again}{see \code{simple.monitoring}.}
    \item{save.as}{.}
}

\value{
    Invisibly return latest data for next comparison. This function is
    run mainly for its side effects.}
\details{
    This function allows for 
       -combining forecasts (ie. monitoring or other forecast data)
       -input (policy) projections
    If these feature are not need see simple.monitoring.
    \code{mail.list} and \code{error.mail.list} should be single strings (not 
    vectors) but the string can contain multiple user ids for mail.
    If \code{overriding.data.names=NULL} then no overriding data is used.
    report.variables indicates output variables which are reported. If NULL,
    then all outputs are reported.
    show.end is min of this and \code{overriding.data$input} if needed.

    This functions allow for the use of over-riding data
    which may come from other forecasts or monitoring and 
    can be used to augment (and replace) actual data.
    Overriding data is used in place of data and model forecasts to 
    the horizon for which it is available. 
    Also, input (policy) variable forecasts can be used.
    NB. The combination is not in the sense of averaging
    together forecasts.    

}
\seealso{
    \code{\link{simple.monitoring}}
}

%\keyword{DSE}
\keyword{ts}


