\name{ms.rep}
\Rdversion{1.2}
\alias{ms.rep}
\alias{meanshift}
\alias{ms.rep.min}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean shift procedures.
}
\description{
Functions for mean shift, iterative mean shift, and inverse mean shift. 
}

\usage{
meanshift(X, x, h)
ms.rep(X, x, h, thresh= 0.0001, iter=200)
ms.rep.min(X, x, h, thresh=0.000001, iter=200, adjust.convergence=FALSE, verbose=TRUE)
}

\arguments{
  \item{X}{data matrix or vector.}
  \item{x}{point from which we wish to shift to the local mean.}
  \item{h}{scalar or vector-valued bandwidth; see also description in \code{\link{ms}}.}
  \item{thresh, iter}{mean shift iterations are stopped when the
  mean shift length (relative to the distance of of \code{x} to the overall mean; see Note section) falls below \code{thresh}, or  after \code{iter} iterations (whatever event happens first).} 
\item{adjust.convergence}{Only required for estimation of antimodes via \code{ms.rep.min}. See also Note section below.}  
\item{verbose}{Allows or suppresses text output messages.} 
}

\details{
The function \code{meanshift} computes a single mean shift iteration, and \code{ms.rep} computes an iterative series of mean shift iterations. Both of these functions are rarely used on their own, but are typically called by the overarching function \code{\link{ms}}.

The function \code{ms.rep.min} implements an inverse version of the mean shift procedure which can be used for the computation of antimodes (minima of the density). The methodology has been presented in the univariate setting  in Ameijeiras-Alonso and Einbeck (2023). The function has been tested to a moderate extent in two dimensions, and it can be expected to work here. It has not been tested for higher dimensions. That is, for data of dimension 3 or higher, the function may or may not produce a result which may or may not correspond to an antimode. 

}


\value{
 The function \code{meanshift} delivers a single (vector-valued) value.
 
The functions \code{ms.rep} and \code{ms.rep.min} produce a list with the following items:

\item{Meanshift.points}{(called \code{M} for \code{ms.rep.min}); the trajectory of points found while proceeding from the starting value \code{x} to the mode (or antimode, respectively)}
\item{Threshold.values}{These give the iteration-wise values of the relative length of the mean shift step (explained in the Note section) which are then compared to \code{thresh}.}
\item{start}{The starting value.}
\item{final}{The mode or antimode, respectively.}
  
}


\references{
Ameijeiras-Alonso, J. and Einbeck, J. (2023). A fresh look at mean-shift based modal
clustering. Advances in Data Analysis and Classification, \doi{10.1007/s11634-023-00575-1}.
}

\author{
  J. Einbeck. See \code{\link{LPCM-package}} for further
  acknowledgements. 
}

\note{
The threshold \code{thresh} for stopping mean shift iterations works as follows. At each iteration, we compare
the length of the mean shift, that is the Euclidean distance between the point \code{x} and its local mean \code{m}, to the Euclidean distance between the point \code{x} and the overall data mean. If this distance falls below \code{thresh}, the mean shift procedure is stopped.

When \code{ms.rep} is called by function \code{ms},  the relation of the thresholds \code{thr} and \code{thresh} is \code{thresh = thr^2}.

Convergence of the inverse mean shift algorithm is not mathematically guaranteed. Of course, no antimode will be found if there is none (i.e., if the value \code{x} is not situated between two modes), in which case the method will return a NA.
But the algorithm may also fail to converge if the antimode is associated with a very small density value, which however rarely happens in practice unless the two distributions are fully separated. In this case the inverse mean shift algorithm will oscillate around the antimode. By reducing the step length successively once  that such a situation is identified,  convergence can still be ensured algorithmically. This adjustment is activated if the
 option \code{check.convergence} is set equal to \code{TRUE}.  This functionality is experimental, and details are to be reported elsewhere.
}

%\section{Warning}{
%
% } 


\seealso{
\code{\link{ms}}
}

\examples{
data(stamps, package="multimode")
h0 <- 0.005
hist(stamps, breaks=20)
# Take arbitrary starting value x=0.08, sitting between a mode and antimode
mode <- ms.rep(stamps, 0.08,h0)$final
antimode <- ms.rep.min(stamps, 0.08,h0, verbose=FALSE)$final
segments(mode, 0, mode, 100, col=2, lwd=3)
segments(antimode, 0, antimode,100, col=3, lwd=3)
}
\keyword{multivariate}% __ONLY ONE__ keyword per line
