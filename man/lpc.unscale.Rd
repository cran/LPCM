\name{lpc.unscale}
\Rdversion{1.1}
\alias{lpc.unscale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Unscaling LPC objects.
}
\description{
\code{lpc.unscale} takes a fitted object of type \code{lpc} or \code{lpc.spline}, which had been fitted using option \code{scaled=TRUE}, and transforms the scaled components back to the original data scale. 
}

\usage{
lpc.unscale(object) 
}

\arguments{
  \item{object}{Object of class \code{lpc} or \code{lpc.spline}.}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
A list with items \code{LPC}, \code{data}, \code{starting.points}, \code{knots.coords}, \code{closest.coords}, which gives the unscaled versions of these quantities (some of them may carry the value \code{NULL}, depending on which object was given as the function argument). 
}


\references{
Einbeck, Tutz, and Evers (2005). Local principal curves. Statistics and Computing 15, 301-313.
}

\author{
  JE
}
%\note{
%
%}

%\section{Warning}{
%
% } 


\seealso{
\code{\link{lpc}}, \code{\link{lpc.spline}}
}

\examples{
data(gvessel)
lpc.unscale(lpc(gvessel[,c(2,4,5)],  h=0.11,  x0=c(35, 1870, 6.3)) )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}% __ONLY ONE__ keyword per line
