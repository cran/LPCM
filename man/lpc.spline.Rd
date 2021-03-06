\name{lpc.spline}
\Rdversion{1.1}
\alias{lpc.spline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Representing local principal curves through a  cubic spline.
}
\description{
Fits a natural cubic spline component-wise through the series of local centers of mass. This provides a continuous parametrization in terms of arc length distance, which can be used to compute a projection index for the original or new data points.
}
\usage{
lpc.spline(lpcobject, optimize = TRUE, compute.Rc=FALSE,
     project=FALSE, ...)
}


\arguments{
  \item{lpcobject}{Object of class \code{lpc}.}
  \item{optimize}{Boolean. If TRUE, \code{optimize} is used to find the point on the curve with minimum distance. Otherwise, data points are only projected onto the closest knot.}
  \item{compute.Rc}{Boolean. If TRUE, the goodness-of-fit measure
  \eqn{R_C}{R_c} suggested in [1] is computed  and returned (using the scaled data, if \code{scaled=TRUE} in \code{lpcobject}).  }
  \item{project}{Boolean. If TRUE, projections onto curve are computed.}
  \item{...}{Additional arguments to be passed to \code{lpc.project.spline}}
}


\details{
See  reference [2].
}

\value{
  \item{knots.pi}{LPC parameters (in cubic spline parametrization) at
    position of the knots of the spline function (these are not
    identical to the LPC mass points!)}
\item{knots.coords}{Coordinates of the spline knots. }
\item{closest.pi}{Parameter of the projected data points.}
\item{closest.coords}{Coordinates of projected data points.}
\item{closest.dist}{Euclidean distance between original and projected data point.}
\item{closest.branch}{ID Number of the branch on which the data point
  was 
  projected (the IDs are given in the output of function \code{lpc}).}
\item{Rc}{Value of \eqn{R_C}{R_c}. }
\item{project}{repeats the input value of  \code{project}.}
\item{lpcobject}{returns the provided object \code{lpcobject}.}
\item{splinefun}{returns the cubic spline function (generated by
  \code{lpc.splinefun}).}
}
\references{
[1] Einbeck, J., Tutz, G.,  and Evers, L. (2005). Local principal curves. Statistics and Computing 15, 301-313.

[2] Einbeck, J., Evers, L. & Hinchliff, K. (2010):  Data compression and regression based on local principal curves. In A. Fink, B. Lausen, W. Seidel, and A. Ultsch (Eds),  Advances in Data Analysis, Data Handling, and Business Intelligence, Heidelberg, pp. 701--712, Springer. 
}
\author{
  J. Einbeck and L. Evers
}
\note{
The parametrization of the cubic spline function is not exactly the same as that of the original LPC. The reason is that the latter uses Euclidean distances between centers of masses, while the former uses the arc length along the cubic spline. However, the differences are normally quite small.  
}

\section{Warning}{
Careful with options \code{project} and \code{compute.Rc} - they can take rather long
if the data set is large!
 } 


\seealso{
\code{\link{lpc}}
}

\examples{
data(gvessel)
gvessel.lpc <- lpc(gvessel[,c(2,4,5)],   h=0.11,  x0=c(35, 1870, 6.3))
gvessel.spline  <- lpc.spline(gvessel.lpc)
plot(gvessel.spline, lwd=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smooth}
\keyword{multivariate}% __ONLY ONE__ keyword per line
