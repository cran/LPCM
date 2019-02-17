\name{print.lpc}
\alias{print.lpc}
\alias{print.lpc.spline}
\alias{print.ms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Printing output for lpc, lpc.spline, and ms objects}
\description{
  Takes an object of class \code{lpc}, \code{lpc.spline}, \code{ms} and displays some standard output.
}

\usage{
\method{print}{lpc}( x, digits = max(3, getOption("digits") - 3), ...) 
\method{print}{lpc.spline}( x, digits = max(3, getOption("digits") - 3), ...) 
\method{print}{ms}( x, digits = max(3, getOption("digits") - 3), ...) 
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{lpc}, \code{lpc.spline}, or \code{ms}.}
  \item{digits}{not yet in use.}
  \item{\dots}{further arguments.}
}
%\details{
% Plots a local principal curve  
%}
\value{
 Some short text.
 }



%\references{
%Einbeck, Tutz, and Evers (2005). Local principal curves. Statistics and Computing 15, 301-313.
%
%Einbeck, J., Evers, L. & Hinchliff, K. (2010):  Data compression and regression based on local principal curves. In A. Fink, B. Lausen, W. Seidel, and A. Ultsch (Eds),  Advances in Data Analysis, Data Handling, and Business Intelligence, Heidelberg, pp. 701--712, Springer. 
%}



\author{ JE }
\seealso{ 
\code{\link{lpc}}, \code{\link{ms}}   
}

\examples{
data(calspeedflow)
lpc1 <- lpc(calspeedflow[,3:4])
print(lpc1)
lpc2 <- lpc.spline(lpc1)
print(lpc2)
ms1<- ms(calspeedflow[,3:4], plot=FALSE)
print(ms1)
}

%\keyword{ smooth}% at least one, from doc/KEYWORDS
\keyword{multivariate}% __ONLY ONE__ keyword per line
