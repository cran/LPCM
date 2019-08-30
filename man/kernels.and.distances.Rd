\name{kernels.and.distances}
\Rdversion{1.1}
\alias{kern}
\alias{kernd}
\alias{kdex}
\alias{vecdist}
\alias{mindist}
\alias{distancevector}
\alias{enorm}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Auxilary  kernel and distance functions.}
\description{
 Internal  \pkg{LPCM} functions which are normally not to be called by the user.
}
\usage{
kern(y, x = 0, h = 1)
kernd(X, x, h)
kdex(X, x, h)
distancevector(X, y, d = "euclid", na.rm = TRUE)
vecdist(X,Y)
mindist(X,y)
enorm(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a number or vector.}
  \item{y}{a vector.}
  \item{h}{a bandwidth.}
  \item{X}{a matrix.}
  \item{Y}{a matrix.}
\item{d}{type of distance measure (only `euclid').}
\item{na.rm}{...}
}

\details{
\code{kern} specifies the base kernel (by default Gaussian) used in
\code{lpc} ; \code{kernd} is the corresponding multivariate product
kernel.  \code{kdex} is a pointwise multivariate kernel density estimator.
  
\code{distancevector} makes use of function \code{vdisseuclid} from \R package \pkg{hopach} (but that package does not need to be loaded). \code{enorm} is the Euclidean norm.
}

\references{
Pollard, van der Laan, and Wall (2010). Hierarchical Ordered Partitioning and Collapsing Hybrid (HOPACH). \R package \pkg{hopach} version 2.9.1.
}

\author{
JE
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{smooth}
%\keyword{multivariate}% __ONLY ONE__ keyword per line
