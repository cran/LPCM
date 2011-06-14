\name{Rc}
\Rdversion{1.1}
\alias{Rc}
\alias{lpc.Rc}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Goodness-of-fit for principal objects
}
\description{
These functions compute the goodness-of-fit criterion \eqn{R_C}{R_c}
proposed in Einbeck, Tutz, and Evers (2005), and extended beyond the
scope of principal curves in Einbeck (2011).
}
\usage{
Rc(data,  closest.coords, type="curve")
lpc.Rc(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A data matrix.}
  \item{object}{An object of class \code{lpc} or \code{lpc.spline}.}
\item{closest.coords}{A matrix of coordinates of the projected data.}
\item{type}{For principal curves, don't modify. For principal points,
  set "points". }
}

\details{
\code{Rc} computes the value \eqn{R_C}{R_c}, a quantity which
estimates the goodness-of-fit of a fitted principal curve.   This
quantity can be interpreted similar to the coeffient of determination in
regression analysis: Values close to 1 indicate a good fit, while values
close to 0 indicate a `bad' fit (corresponding to linear PCA).

In principle, function \code{Rc} can be used for assessing
goodness-of-fit of any principal curve algorithm provided that
the coordinates (\code{closest.coords}) of the projected data are
available (For instance, for HS principal curves fitted via
\code{princurve}, this information is contained in component \code{$s},
and for LPCs this information is given in component \code{$closest.coords} of the
spline representation). It can also be used for assessing the
goodness-of-fit of algorithms which find "principal points" (such as
iterative mean shift, or k-means); set \code{type="points"} in
this case (see also help file to \code{\link{ms}}).

\code{lpc.Rc} is  a wrappper around \code{Rc}, which takes an object of
type \code{lpc} or \code{lpc.spline}. This function computes all missing information, so computation will take the longer the less informative the given \code{object} is.

If the data were scaled, then the scaled data and results should also be
used as arguments in \code{Rc}.   The function \code{lpc.Rc} looks up the option \code{scaled} in the fitted
object, and cares for this automatically.  Important: If the data
were scaled, then do NOT unscale the results by hand in order to feed
the unscaled version into  \code{Rc}, this will give a wrong result.

In terms of methodology, these functions compute  \eqn{R_c}{R_c} directly through the mean
reduction of absolute residual length, rather than through the
area above the coverage curve.

These functions do currently not account for observation
weights, i.e.  \eqn{R_C}{R_c} is computed through the unweighted mean
reduction in absolute residual length (even if weights have been used for
the curve fitting).  
}

\references{
Einbeck, Tutz, and Evers (2005). Local principal curves. Statistics and
Computing 15, 301-313.

Einbeck (2011).  Bandwidth selection for nonparametric unsupervised
learning techniques -- a unified approach via self-coverage. Journal of
Pattern Recognition Research, to appear.
}

\author{
J. Einbeck and L. Evers.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lpc.spline}}
}

\examples{
data(calspeedflow)
lpc1 <- lpc(calspeedflow[,3:4])
lpc.Rc(lpc1)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
\keyword{multivariate}% __ONLY ONE__ keyword per line
