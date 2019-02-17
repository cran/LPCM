\name{Rc}
\Rdversion{1.1}
\alias{Rc}
\alias{Rc.lpc}
\alias{Rc.lpc.spline}
\alias{Rc.ms}
\alias{base.Rc}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Measuring goodness-of-fit for principal objects.
}
\description{
These functions compute the `coverage coefficient' \eqn{R_C}{R_c}
 for local principal curves, local principal points
 (i.e., kernel density estimates obtained through iterated mean shift), and other principal objects.  
}
\usage{
Rc(x,...)

\method{Rc}{lpc}(x,...)
\method{Rc}{lpc.spline}(x,...)
\method{Rc}{ms}(x,...)

base.Rc(data,  closest.coords, type="curve")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object used to select a method.}
  \item{...}{Further arguments passed to or from other methods (not
    needed yet).}
  \item{data}{A data matrix.}
 \item{closest.coords}{A matrix of coordinates of the projected data.}
 \item{type}{For principal curves, don't modify. For principal points,
  set "points". }
}

\details{
\code{Rc} computes the coverage coefficient \eqn{R_C}{R_c}, a quantity which
estimates the goodness-of-fit of a fitted principal object.   This
quantity can be interpreted similar to the coeffient of determination in
regression analysis: Values close to 1 indicate a good fit, while values
close to 0 indicate a `bad' fit (corresponding to linear PCA).


For objects of type \code{lpc}, \code{lpc.spline}, and \code{ms},  S3 methods are available which use the generic function
\code{Rc}.  This, in turn, calls the base function \code{base.Rc}, which
can also be used manually if the fitted object is of another class.
In principle, function \code{base.Rc} can be used for assessing
goodness-of-fit of any principal object provided that
the coordinates (\code{closest.coords}) of the projected data are
available. For instance, for HS principal curves fitted via
\code{princurve}, this information is contained in component \code{$s},
and for a a k-means object, say \code{fitk}, this information can be
obtained via \code{fitk$centers[fitk$cluster,]}. Set \code{type="points"} in
the latter case.

The function \code{Rc} attempts to compute all missing information, so
computation will take the longer the less informative the given
object \code{x} is. Note also,  \code{Rc} looks up the option \code{scaled} in the fitted
object, and accounts for the scaling automatically. Important: If the data
were scaled, then do NOT unscale the results by hand in order to feed
the unscaled version into  \code{base.Rc}, this will give a wrong result.

In terms of methodology, these functions compute  \eqn{R_C}{R_c} directly through the mean
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
Pattern Recognition Research  6, 175-192. 
}

\author{
J. Einbeck.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\seealso{
\code{\link{lpc.spline}}, \code{\link{ms}}, \code{\link{coverage}}.
}

\examples{
data(calspeedflow)
lpc1 <- lpc.spline(lpc(calspeedflow[,3:4]), project=TRUE)
Rc(lpc1)
# is the same as:
base.Rc(lpc1$lpcobject$data, lpc1$closest.coords)

\dontrun{
ms1 <- ms(calspeedflow[,3:4], plot=FALSE)
Rc(ms1)
# is the same as:
base.Rc(ms1$data, ms1$cluster.center[ms1$closest.label,], type="points")
}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
\keyword{multivariate}% __ONLY ONE__ keyword per line
