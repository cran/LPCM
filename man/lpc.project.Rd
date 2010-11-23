\name{lpc.project}
\Rdversion{1.1}
\alias{lpc.project}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Projection onto LPC
}
\description{
Projects a new observation onto the spline representation of the local principal curve.
}
\usage{
lpc.project(object, newdata, ...)
}

\arguments{
  \item{object}{ Object of class \code{lpc} or \code{lpc.spline}.}
  %\item{optimize}{Boolean. If TRUE, \code{optimize} is used to find the closest point on the (spline respresentation of the) principal curve. Otherwise, the point can only be projected to the nearest knot.
  \item{newdata}{A data frame containing the new data to be projected.}
  \item{...}{Additional arguments to be passed to \code{lpc.project.spline}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{closest.pi}{Projection index of projected point(s) (in cubic spline parametrization).}
\item{closest.or.pi}{Projection index of projected point(s) (in terms of the original LPC parametrization). }
\item{closest.coords}{Coordinates of projected data point(s)}
\item{closest.dist}{Euclidean distance between data point(s) and their projected counterpart(s).} 
\item{closest.branch}{ID of branch onto which the data point was projected (the IDs are allocated in the output of function \code{lpc}).}
}

\references{
Einbeck, J., Evers, L. & Hinchliff, K. (2010):  Data compression and regression based on local principal curves. In A. Fink, B. Lausen, W. Seidel, and A. Ultsch (Eds),  Advances in Data Analysis, Data Handling, and Business Intelligence, Heidelberg, pp. 701--712, Springer. 
}
\author{
  J. Einbeck and L. Evers
}
\note{
The parametrization of the cubic spline function is not exactly the same as that of the original LPC. The reason is that the latter uses Euclidean distances between centers of masses, while the former uses the arc length along the cubic spline. The differences are normally quite small, though.  
}

%\section{Warning}{
%Careful with options project and compute.Rc - they can take rather long!!!
% } 

\seealso{
\code{\link{lpc}}, \code{\link{lpc.spline}}
}

\examples{
data(gvessel)
gvessel.lpc <- lpc(gvessel[,c(2,4,5)], scaled=TRUE,   h=0.11,  x0=c(35, 1870, 6.3))
lpc.project(gvessel.lpc, newdata=data.frame(salg=35,dephtg= 2000,oxyg=6))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}% __ONLY ONE__ keyword per line
