\name{curvature.stats}
\alias{curvature.stats}

\title{Curvature Statistics Utility Used by Curvature and Project}
\description{xxx}
\usage{
    curvature.stats(cur, n, signif=0.05)
}
\arguments{
    \item{cur}{the relative curvature array. See \code{rel.curvature}.}
    \item{n}{\code{n-p} is denominator degrees of freedom for F statistic, 
       where \code{p} is \code{ncol(cur)}.}
    \item{signif}{The significance level for F test (passed as
       \code{1 - signif} to \code{qf}).}
}
\value{x}
\details{
 ...
}
\references{ 
Bates and Watts(1983), 'Nonlinear Regression Analysis and Its Applications.'
}

\seealso{
    \code{\link{curvature}}
    \code{\link{project}}
    \code{\link{effective.curvature}}
    \code{\link{curvature.Darray}}
    \code{\link{rel.curvature}}
}
%\keyword{DSE}
\keyword{internal}


