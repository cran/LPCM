\name{ms}
\Rdversion{1.1}
\alias{ms}
\alias{meanshift}
\alias{ms.rep}
\alias{ms.self.coverage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean shift clustering.
}
\description{
Functions for mean shift, iterative mean shift,  mean shift clustering,
and bandwidth selection for mean shift clustering (based on
self-coverage). The main function is \code{ms} which, for a given
bandwidth, detects the local modes (`local principal points') and performs the clustering.

These  functions implement the techniques presented in Einbeck (2011).
}

\usage{
meanshift(X, x, h)
ms.rep(X, x, h, plotms=1, thresh= 0.00000001, iter=200)
ms(X, h, subset,  thr=0.0001, scaled= TRUE, iter=200, plotms=2,
       or.labels=NULL, ...)
ms.self.coverage(X, taumin=0.02, taumax=0.5, gridsize=25,
       thr=0.0001, scaled=TRUE, cluster=FALSE, plot.type="o", 
       or.labels=NULL, print=FALSE, ...)
}

\arguments{
  \item{X}{data matrix.}
  \item{h}{bandwidth.}
  \item{x}{point from which we wish to shift to the local mean.}
  \item{subset}{vector specifying a subset of 1:n, where n is the sample
  size. This allows to run the iterative mean shift procedure only
  from a subset of points (if unspecified, 1:n is used here,
  i.e. each data point serves as a starting point).} 
\item{scaled}{logical.}
  \item{taumin,taumax,gridsize}{determine the grid of bandwidths to
    investigate.}
  \item{thresh, iter}{mean shift iterations are stopped when the
  mean shift length (relative to the length of \code{x}) falls below
      \code{thresh}, or  after \code{iter} iterations (whatever event
  happens first).} 
  \item{thr}{adjacent mean shift clusters are merged if their relative
  distance falls below this threshold.}
  \item{cluster}{if \code{TRUE}, distances are always measured to the
    cluster to which an observation is assigned, rather than to the
    nearest cluster.}
  \item{plotms, plot.type, or.labels, ...}{graphical parameters.}
  \item{print}{if TRUE, coverage values are printed on the screen as soon as
  computed. This is quite helpful especially if \code{gridsize} is large.} 
}

\details{
  The methods implemented here can be used for density mode estimation,
  clustering, and the selection of starting points for the LPC algorithm.

 Chen (1995) showed that, if the mean shift is computed iteratively, the
 resulting sequence of local means converges to a mode of the estimated
 density function.  By assigning each data point to the mode to which it
 has converged, this turns into a clustering technique.

The concepts of coverage and self-coverage, which were originally
introduced in the principal curve context, adapt straightforwardly to
this setting.

The goodness-of-fit messure \code{Rc} can also be  applied in this context. For
instance, a value of \eqn{R_C=0.8}{R_C=0.8} means that,
after the clustering, the mean absolute residual length has been
reduced by \eqn{80\%}{80\%} (compared to the distances to the overall mean).
}


\value{
 
  The main function \code{ms} produces an object of class \code{ms},
  with components:
  
   \item{cluster.center}{a matrix which gives the coordinates of the
  estimated density modes (i.e., of the mean-shift based cluster centers).}
  \item{cluster.label}{assigns each data point to the cluster center to
  which its mean shift trajectory has converged. }
  \item{closest.label}{assigns each data point to the closest cluster
  center in terms of Euclidean distance.}
\item{data}{the data frame (scaled if \code{scaled=TRUE}).}
\item{scaled}{boolean.}
\item{scaled.by}{the data were scaled by dividing each variable through
  the values provided in this vector.}

 For all other functions, use \code{names()}.
}


\references{
Chen, Y. (1995).  Mean Shift, Mode Seeking, and Clustering. IEEE
Transactions on Pattern Analysis and Machine Intelligence, 17, 790-799.

Einbeck, J. (2011). Bandwidth selection for mean-shift based
unsupervised learning techniques: a unified approach via
self-coverage.  Journal of Pattern Recognition Research 6, 175-192. 
}

\author{
  J. Einbeck. See \code{\link{LPCM-package}} for further
  acknowledgements. 
}
%\note{
%
%}

%\section{Warning}{
%
% } 


\seealso{
\code{\link{Rc}}, \code{\link{lpc.self.coverage}}
}

\examples{
data(faithful)
foo <- ms.self.coverage(faithful,gridsize= 10, taumin=0.1, taumax=0.5,
    plot.type="o")    # need higher gridsizes in practice!
h <- select.self.coverage(foo)$select
fit <- ms(faithful,h=h[1])
coverage(fit$data, fit$cluster.center)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}% __ONLY ONE__ keyword per line
