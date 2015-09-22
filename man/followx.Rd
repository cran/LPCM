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
followx(Xi, x0, h, t0, iter, way, weights, pen = 2, phi = 1,
    lasteigenvector = 0, rho0 = 0.4, boundary=0.005,
    convergence.at= 0.000001, cross=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xi}{
%%     ~~Describe \code{Xi} here~~
}
  \item{x0}{
%%     ~~Describe \code{x0} here~~
}
  \item{h}{
%%     ~~Describe \code{h} here~~
}
  \item{t0}{
%%     ~~Describe \code{t0} here~~
}
  \item{iter}{
%%     ~~Describe \code{iter} here~~
}
  \item{way}{
%%     ~~Describe \code{way} here~~
}
  \item{weights}{
}
  \item{pen}{
%%     ~~Describe \code{pen} here~~
}
  \item{phi}{
%%     ~~Describe \code{phi} here~~
}
  \item{lasteigenvector}{
%%     ~~Describe \code{lasteigenvector} here~~
}
  \item{rho0}{
%%     ~~Describe \code{rho0} here~~
}
  \item{boundary}{
%%     ~~Describe \code{plotlpc} here~~
}
  \item{convergence.at}{
%%     ~~Describe \code{plotlpc} here~~
}
  \item{cross}{
%%     ~~Describe \code{plotlpc} here~~
}

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