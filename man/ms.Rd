\name{ms}
\Rdversion{1.2}
\alias{ms}
\alias{meanshift}
\alias{ms.rep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean shift clustering.
}
\description{
Functions for mean shift, iterative mean shift,  and mean shift clustering.
The main function is \code{ms} which, for a given bandwidth, detects the local modes 
(`local principal points') and performs the clustering.

These  functions implement the mean shift techniques presented in Einbeck (2011).
}

\usage{
meanshift(X, x, h)
ms.rep(X, x, h, thresh= 0.000001, iter=200)
ms(X, h, subset,  thr=0.001, scaled= 1, iter=200, plot=TRUE, ...)
}

\arguments{
  \item{X}{data matrix or vector.}
  \item{h}{scalar or vector-valued bandwidth (by default, 5 percent of
  the  data range, or 20 percent of the standard deviation, respectively, in each direction). If set manually and \code{scaled>0}, this
  bandwidth needs to be set on the scaled scale; for instance setting
  \code{scaled=1} and \code{h=0.10} will use a bandwidth of \eqn{10} percent of the data range in either direction. }
  \item{x}{point from which we wish to shift to the local mean.}
  \item{subset}{vector specifying a subset of 1:n, where n is the sample
  size. This allows to run the iterative mean shift procedure only
  from a subset of points (if unspecified, 1:n is used here,
  i.e. each data point serves as a starting point).} 
\item{scaled}{if equal to 1 (default), each variable is divided by its range, and if equal to 2 (or any other positive value other than 1), each variable is divided by its standard deviation. If equal to 0, then no scaling is applied.}
   \item{thresh, iter}{mean shift iterations are stopped when the
  mean shift length (relative to the length of \code{x}) falls below
      \code{thresh}, or  after \code{iter} iterations (whatever event
  happens first).} 
  \item{thr}{adjacent mean shift clusters are merged if their relative
  distance falls below this threshold.}
  \item{plot}{if equal to 0, then no plotted output. For bivariate
  data, \code{plot=1} gives by default a dynamically created color plot showing the mean
  shift trajectories and the resulting clustering.} 
  \item{...}{further graphical parameters.}
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
\item{scaled}{the user-supplied value, could be boolean or numerical.}
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
\note{
All values provided in the output refer to the scaled data, unless \code{scaled=0} or (equivalently) \code{scaled=FALSE}. 

The default option \code{scaled=1} or  \code{scaled=TRUE}  scales the data by dividing each variable through their range (differing from the scaling through the standard deviation as common e.g. for  PCA).  All other settings  \code{scaled>0} will scale the data by their standard deviation.

If \code{scaled=1} or if no scaling is applied, then the default bandwidth setting is 5 percent of the data range in each direction. If the data are scaled through the standard deviation, then the default setting is 20 percent of the standard deviation in each direction. 
}

%\section{Warning}{
%
% } 


\seealso{
\code{\link{Rc}}, \code{\link{plot.ms}}
}

\examples{
data(faithful)
# Mean shift clustering with default bandwidth (5 percent of data range)
ms(faithful)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
\keyword{multivariate}% __ONLY ONE__ keyword per line
