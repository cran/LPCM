\name{LPCM-package}
\Rdversion{1.1}
\alias{LPCM-package}
\alias{LPCM}
\docType{package}
\title{
Local principal curve methods
}
\description{
   Fitting multivariate data patterns with local principal curves;
   including simple tools for data compression (projection), bandwidth
   selection, and measuring goodness-of-fit.

   This packages implements the techniques introduced in Einbeck, Tutz
   & Evers (2005), and successive related papers.

   The main functions to be called by the user are
   \itemize{
   \item \code{\link{lpc}}, for the estimation of the local centers of mass
   which make up the principal curve;
   \item \code{\link{lpc.spline}}, which is a smooth and fully parametrized
   cubic spline respresentation of the latter; 
   \item \code{\link{lpc.project}}, which enables to compress data by
   projecting them orthogonally onto the curve; 
   \item \code{\link{lpc.coverage}} and \code{\link{lpc.Rc}} for assessing
   goodness-of-fit;
   \item \code{\link{lpc.self.coverage}} for bandwidth selection;
   \item the generic \code{plot} and \code{print} functions for objects
   of class \code{lpc} and \code{lpc.spline}. 
   }
This package also contains some (rather experimental) code for density
mode detection and mean shift clustering, see \code{\link{ms}}.

A second R package which will implement the extension of local principal
curves to local principal surfaces and manifolds, as proposed in
Einbeck, Evers & Powell (2010), is in preparation.

}
\details{
\tabular{ll}{
Package: \tab LPCM\cr
Type: \tab Package\cr
License: \tab GPL 2.0 or newer.\cr
LazyLoad: \tab yes\cr
}
}

\author{
Jochen Einbeck and Ludger Evers

Maintainer: Jochen Einbeck <jochen.einbeck@durham.ac.uk>

}
\references{
Einbeck, J., Tutz, G., & Evers, L. (2005), Local principal curves, Statistics and Computing 15, 301-313.

Einbeck, J., Evers, L., & Powell, B. (2010): Data compression and regression through local principal curves and surfaces, International Journal of Neural Systems, 
20, 177-192.
}
\keyword{smooth}
\keyword{multivariate}
\seealso{pcurve, princurve
}
%\examples{
%
%}
