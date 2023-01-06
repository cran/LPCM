\name{ms}
\alias{ms}
\Rdversion{1.2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean shift clustering.
}
\description{
Function for mean shift clustering, which, for a given bandwidth, detects the local modes and performs the clustering.
}

\usage{
ms(X, h, subset,  thr=0.01, scaled= 1, iter=200, plot=TRUE, ...)
}

\arguments{
  \item{X}{data matrix or vector.}
  \item{h}{scalar or vector-valued bandwidth (by default, 5 percent of
  the  data range, or 20 percent of the standard deviation, respectively, in each direction). If set manually and \code{scaled>0}, this
  bandwidth needs to be set on the scaled scale; for instance setting scale; for instance \code{scaled=1} and \code{h=0.10} will use a bandwidth of \eqn{10} percent of the  data range in either direction. }
  \item{subset}{vector specifying a subset of 1:n, where n is the sample
  size. This allows to run the iterative mean shift procedure only
  from a subset of points (if unspecified, 1:n is used here,
  i.e. each data point serves as a starting point).} 
  \item{thr}{adjacent mean shift clusters are merged if their relative
  distance falls below this threshold (see Note section).}
 \item{scaled}{if equal to 1 (default), each variable is divided by its range, and if equal to 2 (or any other positive value other than 1), each variable is divided by its standard deviation. If equal to 0, then no scaling is applied.}
   \item{iter}{maximum mean shift iterations (passed to \code{ms.rep}).}
  \item{plot}{if equal to 0, then no plotted output. For bivariate
  data, \code{plot=1} gives by default a dynamically created color plot showing the mean
  shift trajectories and the resulting clustering.} 
  \item{...}{further graphical parameters.}
}

\details{
The methods implemented here can be used for density mode estimation,
clustering, and the selection of starting points for the LPC algorithm.

It can be shown (Chen, 1995, Comaniciu & Meer, 2002, Li, 2005) that, if the mean shift is computed iteratively, the
resulting sequence of local means converges to a mode of the estimated
density function.  By assigning each data point to the mode to which it
has converged, this turns into a clustering technique.
}


\value{
 
  The function \code{ms} produces an object of class \code{ms},
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
}


\references{
Chen, Y. (1995).  Mean Shift, Mode Seeking, and Clustering. IEEE
Transactions on Pattern Analysis and Machine Intelligence, 17, 790-799.

Comaniciu, D. and Meer,P. (2002). Mean shift: a robust approach toward feature 
space analysis, IEEE Transactions on Pattern Analysis and Machine Intelligence
24, 603-619.

Li, X,  Hu, Z, and Wu, F. (2007). A note on the convergence of the mean shift, Pattern
  Recognition 40, 1756 - 1762.

}

\author{
  J. Einbeck. See \code{\link{LPCM-package}} for further
  acknowledgements. 
}
\note{
All values provided in the output refer to the scaled data, unless \code{scaled=0} or (equivalently) \code{scaled=FALSE}. 

The default option \code{scaled=1} or  \code{scaled=TRUE}  scales the data by dividing each variable through their range (differing from the scaling through the standard deviation as common e.g. for  PCA).  All other settings  \code{scaled>0} will scale the data by their standard deviation.

If \code{scaled=1} or if no scaling is applied, then the default bandwidth setting is 5 percent of the data range in each direction. If the data are scaled through the standard deviation, then the default setting is 20 percent of the standard deviation in each direction. 

The threshold \code{thr} for merging cluster centers works as follows: After identification of a new cluster center, we compute the Euclidean distance of the new center to (each) existing center, relative to the Euclidean distance of the existing center to the overall mean. If this distance falls below \code{thr}, then the new center is deemed identical to the old one.

The goodness-of-fit measure \code{Rc} can also be  applied in this context. For
instance, a value of \eqn{R_C=0.8}{R_C=0.8} means that,
after the clustering, the mean absolute residual length has been
reduced by \eqn{80\%}{80\%} (compared to the distances to the overall mean).
}

%\section{Warning}{
%
% } 


\seealso{
\code{\link{ms.rep}}, \code{\link{Rc}}, \code{\link{plot.ms}}
}

\examples{
data(faithful)
# Mean shift clustering with default bandwidth (5 percent of data range)
ms(faithful)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
\keyword{multivariate}% __ONLY ONE__ keyword per line
