\name{combinationMonitoring}
\alias{combinationMonitoring}

\title{Combination Monitoring}
\description{Automatic monitoring with e-mail of results}
\usage{
    combinationMonitoring(model, data.names,
       previous.data=NULL,
       overriding.data.names=NULL, 
       restrict.overriding.data=TRUE, overriding.horizon=0,
       mail.list=NULL,
       error.mail.list=NULL,
       message.title="Combination Monitoring",
       message.subject="Combination Monitoring",
       message.footnote=NULL,
       show.start= c(0,-3),
       show.end  = c(0,12),    
       report.variables=seriesNames(data.names),
       data.sub.heading=NULL,
       data.tag=" ",
       future.inputData.tag="p",
       overriding.data.tag="m",
       overlapping.period.forecast.tag="g",
       forecast.tag="f",
       run.again=FALSE,
       save.as=NULL)
}
\arguments{
    \item{model}{see \code{simpleMonitoring}.}
    \item{data.names}{see \code{simpleMonitoring}.}
    \item{previous.data}{see \code{simpleMonitoring}.}
    \item{overriding.data.names}{a TSdata (names) object.}
    \item{restrict.overriding.data}{a logical indicating 
        if \code{restrict.overriding.data} should be used to truncate the
	resticion.}
    \item{overriding.horizon}{an integer indicating the horizon of the
        restriction.}

    \item{mail.list}{see \code{simpleMonitoring}.}
    \item{error.mail.list}{see \code{simpleMonitoring}.}
    \item{message.title}{see \code{simpleMonitoring}.}
    \item{message.subject}{see \code{simpleMonitoring}.}
    \item{message.footnote}{see \code{simpleMonitoring}.}
    \item{show.end}{see \code{simpleMonitoring}.}
    \item{show.start}{see \code{simpleMonitoring}.}
    \item{report.variables}{see \code{simpleMonitoring}.}
    \item{data.sub.heading}{see \code{simpleMonitoring}.}
    \item{data.tag}{see \code{simpleMonitoring}.}

    \item{future.inputData.tag}{.}
    \item{overriding.data.tag}{.}
    \item{overlapping.period.forecast.tag}{.}

    \item{forecast.tag}{see \code{simpleMonitoring}.}
    \item{run.again}{see \code{simpleMonitoring}.}
    \item{save.as}{.}
}

\value{
    Invisibly return latest data for next comparison. This function is
    run mainly for its side effects.}
\details{
    This function allows for 
       -combining forecasts (ie. monitoring or other forecast data)
       -input (policy) projections
    If these feature are not need see simpleMonitoring.
    \code{mail.list} and \code{error.mail.list} should be single strings (not 
    vectors) but the string can contain multiple user ids for mail.
    \code{overriding.data.names=NULL} indicates a source for data which should
    be used in place of model forecasts (e.g. preliminary data from a source or
    data from another forecast).
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
    \code{\link{simpleMonitoring}}
}

%\keyword{DSE}
\keyword{ts}


