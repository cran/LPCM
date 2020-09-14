\name{coverage}
\Rdversion{1.1}
\alias{coverage}
\alias{coverage.raw}
\alias{lpc.coverage}
\alias{lpc.self.coverage}
\alias{ms.self.coverage}
\alias{select.self.coverage}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Coverage and self-coverage plots.
}
\description{
These functions compute coverages and self-coverages, and produce corresponding plots, for any principal curve object.  The former may be used as goodness-of-fit measures, and the latter for for bandwidth selection.
}
\usage{
coverage.raw(X, vec, tau, weights=1, plot.type="p", print=FALSE,
      label=NULL,...)

coverage(X, vec, taumin=0.02, taumax, gridsize=25, weights=1,
      plot.type="o", print=FALSE,...)

lpc.coverage(object, taumin=0.02, taumax, gridsize=25, quick=TRUE,
      plot.type="o", print=FALSE, ...)

lpc.self.coverage(X,  taumin=0.02, taumax=0.5,   gridsize=25, x0=1,
     way = "two", scaled=1,  weights=1, pen=2, depth=1,
     control=lpc.control(boundary=0, cross=FALSE),   quick=TRUE,
     plot.type="o", print=FALSE, ... )
 
ms.self.coverage(X, taumin=0.02, taumax=0.5, gridsize=25,
       thr=0.001, scaled=1, cluster=FALSE, plot.type="o", 
       print=FALSE, ...)
       
select.self.coverage(self,  smin, plot.type="o", plot.segments=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{X}{a \eqn{N \times d}{N x d} data matrix.}
\item{object}{An object of type \code{lpc}, \code{lpc.spline} or \code{ms}.}
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
  \item{thr}{adjacent mean shift clusters are merged if their relative
  distance falls below this threshold.}
 \item{cluster}{if \code{TRUE}, distances are always measured to the
    cluster to which an observation is assigned, rather than to the
    nearest cluster.}
\item{self}{An object of class \code{self}, or a matrix with two columns
     providing a self-coverage curve.} 
 \item{smin}{Minimum coverage for bandwidth selection. Default: 1/3 for
   clustering, 2/3 for principal curves.}
\item{plot.type}{If set to 0, no plotted output is given. Otherwise,
  an appropriate plot is provided, using the plotting type as
  specified. }
\item{plot.segments}{A list with default \code{list(lty=c(1,2,3),
    lwd=c(2,1,1),lcol=c(3,3,3))} which specifies how (and how many)
  bandwidth candidates, in order of decreasing negative second derivative of
  self-coverage, are to be highlighted.}
\item{print}{If TRUE, coverage values are printed on the screen as soon as
  computed. This is quite helpful especially if \code{gridsize} is large.} 
\item{x0, way, scaled, pen,  depth, control}{Auxiliary parameters as outlined in
  \code{\link{lpc}}, \code{\link{lpc.control}}, and \code{\link{ms}}.}
\item{...}{Optional graphical parameters passed to the corresponding
  plotting functions.} 
}

\value{
  A list of items, and a plot (unless \code{plot.type=0}).

  The functions \code{lpc.self.coverage} and \code{ms.self.coverage} produce an object of class
  \code{self}.  The component \code{$select} recommends suitable
  bandwidths for the use in \code{lpc}, in the order of strength of
  evidence. These correspond to points of strong negative curvature (implemented via second
  differences) of the self-coverage curve.    

} 

\details{

The function \code{coverage.raw} computes the coverage, i.e. the
proportion of data points lying inside a circle or band with radius
\eqn{\tau}{tau}, for a fixed value \code{tau}.  The whole coverage curve
\eqn{C(\tau)}{C(tau)} is constructed through function \code{coverage}.  
  
Functions \code{coverage.raw} and \code{coverage} can be used for any
object fitted by an unsupervised learning technique (for instance, HS principal curves, or even clustering
algorithms), while the functions prefixing with \code{lpc.} and \code{ms.} can only be
used for the corresponding objects. The functions \code{lpc.coverage} and \code{ms.coverage} are wrappers around 
\code{coverage} which operate directly a fitted object, rather
than a data matrix.

Function \code{select.self.coverage} extracts suitable bandwidths from the
self-coverage curve, and produces a plot. The function is called from
within \code{lpc.self.coverage} or \code{ms.self.coverage}
but can also be called directly by the user (for instance, if the graphical output is to be reproduced, or if
the minimum coverage \code{smin} is to be modified). The component
\code{$select}  contains the selected candidate bandwidths, in the order
of strength of evidence provided by the self-coverage criterion (the
best bandwidth comes first, etc.). A plot is produced as a by-product,
which symbolizes the best bandwidth by a thick solid line, the
second-best by a dashed line, and the third-best by a dotted line. It is
recommended to run the self-coverage functions with fixed starting
points, as in the examples below, and to scale by the range only.

See Einbeck (2011) for details. Note that the original publication by Einbeck, Tutz, and Evers (2005) uses `quick' coverage curves.

}

\references{
Einbeck, J.,  Tutz, G.,  & Evers, L. (2005). Local principal curves. Statistics and Computing 15, 301-313.

Einbeck, J. (2011). Bandwidth selection for mean-shift based
unsupervised learning techniques: a unified approach via
self-coverage. Journal of Pattern Recognition Research 6, 175-192.
}

\author{
J. Einbeck 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lpc}}, \code{\link{ms}}
}

\examples{

data(faithful)
mfit <- ms(faithful)
coverage(mfit$data, mfit$cluster.center, gridsize=18)

\dontrun{
f.self <- ms.self.coverage(faithful,gridsize= 50, taumin=0.1, taumax=0.5, plot.type="o")   
h <- select.self.coverage(f.self)$select
mfit2 <- ms(faithful,h=h[2]) # using `second-best' suggested bandwidth 
}

\dontrun{
data(gvessel)
g.self <-lpc.self.coverage(gvessel[,c(2,4,5)], x0=c(35, 1870, 6.3), print=FALSE, plot.type=0)
h <- select.self.coverage(g.self)$select
g.lfit <- lpc(gvessel[,c(2,4,5)], h=h[1],  x0=c(35, 1870, 6.3))
lpc.coverage(g.lfit, gridsize=10, print=FALSE)
}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
\keyword{multivariate}% __ONLY ONE__ keyword per line
