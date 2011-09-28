\name{unscale}
\Rdversion{1.1}
\alias{unscale}
\alias{unscale.lpc}
\alias{unscale.lpc.spline}
\alias{unscale.ms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Unscaling local principal objects.
}
\description{
\code{unscale} takes an object of type \code{lpc},
\code{lpc.spline}, or \code{ms},  which had been fitted using option
\code{scaled=TRUE}, and transforms the scaled components back to the
original data scale.   
}

\usage{
unscale(x, ...)

\method{unscale}{lpc}(x,...)
\method{unscale}{lpc.spline}(x,...)
\method{unscale}{ms}(x,...)
}

\arguments{
  \item{x}{ an object used to select a method.}
  \item{...}{Further arguments passed to or from other methods (not
    needed yet). }
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{  
A list of relevant items, such as \code{LPC}, \code{start}, \code{cluster.centers},
etc., which gives the unscaled versions of these
quantities (some of them may carry the value \code{NULL}, if the
corresponding information was not available from \code{x}). 
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
\code{\link{lpc}}, \code{\link{lpc.spline}}, \code{\link{ms}}
}

\examples{
data(gvessel)
unscale(lpc(gvessel[,c(2,4,5)],  h=0.11,  x0=c(35, 1870, 6.3)) )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}% __ONLY ONE__ keyword per line
