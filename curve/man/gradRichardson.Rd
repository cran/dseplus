\name{gradRichardson}
\alias{gradRichardson}
\title{Gradient of a Function}
\description{Calculate the gradient of a function by Richardson improvement.}
\usage{
    gradRichardson(func, x, d=0.01, eps=1e-4, r=6, show.details=FALSE)
}
\arguments{
    \item{func}{A function with a real result.}
    \item{x}{
    A real or real vector argument to func, indicating the point at which the
    gradient is to be calculated.}
    \item{d}{the fraction of x to use for the initial numerical approximation.
        The default is set to 0.01*x or eps if x is 0.0.}
    \item{eps}{Used instead of d for elements of x which are zero.}
    \item{r}{The number of Richardson improvement iterations (repetions with 
        successly smaller d).}
    \item{show.details}{logical indicating if detailed calculations should 
        be shown.}
}
\value{A real or real vector of the calculated gradient.}

\details{
 This function calculates a numerical approximation of the first
   derivative of func at the point x. The calculation
   is done by Richardson's extrapolation (see e.g. Linfield and Penny, 1989)
   he method should be used if accuracy, as opposed to speed, is important.

  GENERAL APPROACH
     --  GIVEN THE FOLLOWING INITIAL VALUES:
	     INTERVAL VALUE D, NUMBER OF ITERATIONS R, AND
	     REDUCED FACTOR V.
      - THE FIRST ORDER aproximation to the DERIVATIVE WITH RESPECT TO Xi IS

	   F'(Xi)={F(X1,...,Xi+D,...,Xn) - F(X1,...,Xi-D,...,Xn)}/(2*D)
       
     --  REPEAT r TIMES  with successively smaller D  and 
	  then apply Richardson extraplolation

}
\references{ 
   Linfield, G. R. and Penny, J. E. T. (1989) \emph{Microcomputers in Numerical 
   Analysis}. New York: Halsted Press.
}
\seealso{\code{\link{gradNumerical}}}
\examples{gradRichardson(sin, pi)}
%\keyword{DSE}
\keyword{ts}

