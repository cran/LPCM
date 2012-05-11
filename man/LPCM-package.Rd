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

   This package implements the techniques introduced in Einbeck, Tutz
   & Evers (2005), and successive related papers.

   The main functions to be called by the user are
   \itemize{
   \item \code{\link{lpc}}, for the estimation of the local centers of mass
   which make up the principal curve;
   \item \code{\link{lpc.spline}}, which is a smooth and fully parametrized
   cubic spline respresentation of the latter; 
   \item \code{\link{lpc.project}}, which enables to compress data by
   projecting them orthogonally onto the curve; 
   \item \code{\link{lpc.coverage}} and \code{\link{Rc}} for assessing
   goodness-of-fit;
   \item \code{\link{lpc.self.coverage}} for bandwidth selection;
   \item the generic \code{plot} and \code{print} methods for objects
   of class \code{lpc} and \code{lpc.spline}. 
 }

This package also contains some  code for density
mode detection (`local principal points') and  mean shift clustering (as well as bandwidth
selection in this context), which implements the methods presented in
Einbeck (2011). See the help file for \code{\link{ms}}.

A second R package which will implement the extension of local principal
curves to local principal surfaces and manifolds, as proposed in
Einbeck, Evers & Powell (2010), is in preparation.

}
\details{
\tabular{ll}{
Package: \tab LPCM\cr
Type: \tab Package\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}
}

\author{
Jochen Einbeck and Ludger Evers

Maintainer: Jochen Einbeck <jochen.einbeck@durham.ac.uk>

}



\section{Acknowledgements}{
 Contributions (in form of pieces of code, or useful suggestions for
 improvements) by Jo Dwyer, Mohammad Zayed, and
 Ben Oakley are gratefully acknowledged.
}

\references{
Einbeck, J., Tutz, G., & Evers, L. (2005): Local principal curves, Statistics and Computing 15, 301-313.

Einbeck, J., Evers, L., & Powell, B. (2010): Data compression and regression through local principal curves and surfaces, International Journal of Neural Systems, 
20, 177-192.

Einbeck, J. (2011): Bandwidth selection for nonparametric unsupervised
learning techniques -- a unified approach via self-coverage. Journal of
Pattern Recognition Research 6, 175-192. 
}
\keyword{smooth}
\keyword{multivariate}
\seealso{pcurve, princurve
}
%\examples{
%
%}
