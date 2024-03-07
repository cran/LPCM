\name{followx}
\Rdversion{1.1}
\alias{followx}

\title{
Fit an individual branch of a local principal curve.
}
\description{
Internal function of package \pkg{LPCM} called by \code{lpc}. Do not use! 
}
\usage{
followx(Xi, x0, h, t0, iter, way, weights, pen = 2, 
    lasteigenvector = 0, rho0 = 0.4, boundary=0.005,
    convergence.at= 0.000001, cross=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xi}{
data matrix
}
  \item{x0}{
branch starting point
}
  \item{h}{
bandwidth
}
  \item{t0}{
step length
}
  \item{iter}{
number of iterations (within the given curve branch)
}
  \item{way}{
possible values "one", "two", "back" (in which directions the curve should proceed)
}
  \item{weights}{
  vector of weights
}
  \item{pen}{
power used for angle penalization
}
\item{lasteigenvector}{
to be passed on from \code{lpc}
}
  \item{rho0}{
constant; see reference [1] in \code{\link{lpc.control}}.
}
  \item{boundary}{
boundary correction, see Einbeck and Zayed (2014)
}
  \item{convergence.at}{
convergence parameter
}
  \item{cross}{
Boolean; are curves allowed to cross?
}

}

\references{
Einbeck, J., & Zayed, M. (2014). Some asymptotics for localized principal components and curves. Communications in Statistics - Theory and Methods, 43(8), 1736-1749. doi: 10.1080/03610926.2012.673676
}

\author{
JE
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{lpc}}
}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{multivariate}