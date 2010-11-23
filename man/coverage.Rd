\name{coverage}
\Rdversion{1.1}
\alias{coverage}
\alias{coverage.raw}
\alias{lpc.coverage}
\alias{lpc.self.coverage}
\alias{select.self.coverage}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Coverage and self-coverage plots.
}
\description{
These functions compute coverages (for any principal object),
and self-coverages  (only for local principal curves, these may be
used for bandwidth selection).
}
\usage{
coverage.raw(X, vec, tau, weights=1, plot.type="p", print=FALSE,
      label=NULL,...)

coverage(X, vec, taumin=0.02, taumax, gridsize=25, weights=1,
      plot.type="o", print=TRUE,...)

lpc.coverage(object, taumin=0.02, taumax, gridsize=25, quick=TRUE,
      plot.type="o", print=TRUE, ...)

lpc.self.coverage(X,  taumin=0.02, taumax=0.5,   gridsize=25, x0, mult=1, 
     way = "two", scaled=TRUE,  weights=1, pen=2, depth=1,
     control=lpc.control(boundary=0, cross=FALSE),   quick=TRUE,
     plot.type="o", print=TRUE, ... )

select.self.coverage(self, sens=0.02, from, plot.type="o", auto=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{X}{a \eqn{N \times d}{N x d} data matrix.}
\item{object}{An object of type \code{lpc} or \code{lpc.spline}.}
\item{vec}{A matrix with \eqn{d}{d} columns. The rows contain the points which
  make up the fitted object.}
\item{tau}{tube size.}
\item{taumin}{Minimal tube size.}
\item{taumax}{Maximal tube size.}
\item{weights}{An optional vector of weights.  If weights are specified,
     then the coverage is the weighted mean of the indicator functions
     for falling within the tube. The function \code{lpc.coverage} does not have a
     \code{weights} argument, as it extracts the weights from the
     \code{$weights} component of the fitted \code{object}. }
\item{label}{Experimental option; don't use.}   
\item{gridsize}{The number of different tube sizes to consider.}
\item{quick}{If TRUE, an approximate coverage curve is provided by
  computing distances between data points and the curve through  the
  closest local centers or mass; whereas with FALSE  we use the
  distances of the points when projected orthogonally onto the spline
  representation of the local principal curve. The latter takes considerably
  more computing time. The resulting coverage curves are generally very
  similar, but the quick version may deliver little spurious peaks occasionally.  
}
\item{self}{An object of class \code{self}, or a matrix with two colums
     providing a self-coverage curve.} 
   \item{sens}{Sensitivity parameter: only bandwidths corresponding to
     curvature values of the self-coverage curve, which are  less than the negative value of \code{sens}, are selected.}
 \item{from}{Minimum coverage for bandwidth selection. Default: 1/3 for
   clustering, 2/3 for principal curves.}
 \item{auto}{Unimplemented. If \code{TRUE}, the sensitivity is decreased
   until a bandwidth is selected.}
\item{plot.type}{If set to 0, no plotted output is given. Otherwise,
  an appropriate plot is provided, using the plotting type as specified. }
\item{print}{If TRUE, coverage values are printed on the screen as soon as
  computed. This is quite helpful especially if \code{gridsize} is large.} 
\item{x0, way, scaled, pen, mult, depth, control}{LPC parameters as outlined in
  \code{\link{lpc}} and \code{\link{lpc.control}}.}
\item{...}{Optional graphical parameters passed to the corresponding
  plotting functions.} 
}

\value{
  A list of items, and a plot (unless \code{plot.type=0}).

  For function \code{lpc.self.coverage}, the item \code{$select} recommends suitable
  bandwidths for the use in \code{lpc}. These correspond to
  points of strong negative curvature (implemented via second
  differences) of the self-coverage curve.    

} 

\details{
See Einbeck, Tutz, and Evers (2005). This paper uses `quick' coverage curves.

Functions \code{coverage.raw} and \code{coverage} can be used for any
object fitted by an unsupervised learning technique (for instance, HS principal curves, or even clustering
algorithms), while the functions prefixing with \code{lpc.} can only be
used for local principal curves. The function \code{lpc.coverage} is a wrapper around 
\code{coverage} which takes directly a fitted \code{lpc} object, rather
than a data matrix.

Function \code{select.self.coverage} is called by
\code{lpc.self.coverage}. It extracts suitable bandwidths from the
self-coverage curve, and produces a plot. The function can be called
directly by the user if the sensitivity \code{sens} or the minimum
coverage \code{from} are to be modified. 
}

\references{
Einbeck, J.,  Tutz, G.,  & Evers, L. (2005). Local principal curves. Statistics and Computing 15, 301-313.
}

\author{
J. Einbeck 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lpc}}
}

\examples{
data(gvessel)
gvessel.lpc <- lpc(gvessel[,c(2,4,5)], h=0.11,  x0=c(35, 1870, 6.3))
lpc.coverage(gvessel.lpc, gridsize=10, print=FALSE)
lpc.self.coverage(gvessel[,c(2,4,5)], x0=c(35, 1870, 6.3), print=FALSE)

data(calspeedflow)
fitms <- ms(calspeedflow[,3:4], h=0.1)
coverage(fitms$data, fitms$cluster.center)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
\keyword{multivariate}% __ONLY ONE__ keyword per line
