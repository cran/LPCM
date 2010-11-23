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
self-coverage). These are experimental (and not fully documented) functions which implement the
techniques presented in Einbeck (2010).
}

\usage{
meanshift(X, x, h)
ms.rep(X, x, h, plotms=1, thresh= 0.00000001, iter=100)
ms(X, h, subset,  thr=0.001, scaled= TRUE, plotms=2, or.labels=NULL, ...)
ms.self.coverage(X, taumin=0.02, taumax=0.5, gridsize=25,
       thr=0.001, scaled=TRUE, draw=1/3,  cluster=FALSE, plot.type="o", 
       or.labels=NULL, print=FALSE, ...)
}

\arguments{
  \item{X}{Data matrix.}
  \item{h}{Bandwidth.}
  \item{x}{Point from which we wish to shift to the local mean.}
  \item{subset}{Vector specifying a subset of \code{1:n}. This allows to run the iterative mean shift procedure only
  from a subset of points (if unspecified, \code{1:n} is used here,
  i.e. each data point serves as a starting point).} 
\item{scaled}{Logical.}
  \item{taumin,taumax,gridsize}{Determine the grid of bandwidths to
    investigate.}
  \item{draw}{Only cluster centers belonging to this (randomly selected)
  fraction of the original
  data cloud are used for the computation of the self-coverage.}
  \item{thresh, thr, iter}{Parameters controlling convergence behavior.}
  \item{cluster}{if \code{TRUE}, distances are always measured to the
    cluster to which an observation is assigned, rather than to the
    nearest cluster.}
  \item{plotms, plot.type, or.labels, ...}{Graphical parameters.}
  \item{print}{If TRUE, coverage values are printed on the screen as soon as
  computed. This is quite helpful especially if \code{gridsize} is large.} 
}

\details{
The methods implemented here are not directly related to local principal curves, but
have with them the building block "mean shift" in common.

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
 
  For the function \code{ms}:
  
   \item{cluster.center}{a matrix which gives the coordinates of the
  estimated density modes (i.e., of the mean-shift based cluster centers).}
  \item{cluster.label}{assigns each data point to the cluster center to
  which its mean shift trajectory has converged. }
  \item{closest.label}{assigns each data point to the closest cluster
  center in terms of Euclidean distance.}
\item{data}{the data frame (scaled if \code{scaled=TRUE}).}
\item{scaled.by}{the data were scaled by dividing each variable through
  the values provided in this vector.}

 For all other functions, use \code{names()}.
}


\references{
Chen, Y. (1995).  Mean Shift, Mode Seeking, and Clustering. IEEE
Transactions on Pattern Analysis and Machine Intelligence, 17, 790-799.

Einbeck, J. (2010). Bandwidth selection for mean-shift based
unsupervised learning techniques: a unified approach via
self-coverage. Working paper, Durham University. 
}

\author{
  JE, with contributions by M. Zayed.
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
Rc(fit$data, fit$cluster.center[fit$closest.label,], type="points")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}% __ONLY ONE__ keyword per line
