\name{plot.lpc}
\alias{plot.lpc}
\alias{plot.lpc.spline}
\alias{plot.ms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting local principal curves and mean shift trajectories}
\description{
Takes an object of class \code{lpc}, \code{lpc.spline} or \code{ms}. In the case of principal curves, it plots any subset of the following components of the local principal curve:  Centers of mass; the curve connecting the local centers of mass; the cubic spline representation of the curve; the projections onto the curve; the starting points. For the mean shift procedure, it produces a plot of mean shift trajectories and cluster centers. 
}
\usage{
\method{plot}{lpc}(x, type, unscale = TRUE, lwd = 1, datcol = "grey60", 
    datpch = 21, masscol = NULL, masspch = 15, curvecol = 1, splinecol = 3, 
    projectcol = 4, startcol = NULL,  startpch=NULL,...)  
\method{plot}{lpc.spline}(x, type, unscale = TRUE, lwd = 1, datcol = "grey60", 
    datpch = 21, masscol = NULL, masspch = 15, curvecol = 1, splinecol = 3, 
    projectcol = 4, startcol = NULL,  startpch=NULL,...)
\method{plot}{ms}(x, unscale=FALSE, lwd=1, datcol="grey70",   datpch=21, masscol=NULL, 
    masspch=15, curvecol=NULL, ...)    
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{lpc} or \code{lpc.spline}.   }
  \item{type}{a vector of type \code{c("mass", "spline",...)} with possible entries \code{mass, curve}, \code{spline, project, start.}}
  \item{unscale}{if TRUE, then data (and all fitted componens) are scaled back to their original scale; otherwise the scaled data are plotted (only relevant if \code{scaled=TRUE}  in the fitted object). For \code{ms}, this is currently unimplemented.}  
  \item{lwd}{width of principal curves or trajectories.}
  \item{datcol}{color of data points.}
  \item{datpch}{plotting symbol for data points.}
  \item{masscol}{color of centers of mass (see below) or cluster centers.  }
  \item{masspch}{plotting symbol for centers of mass or cluster centers.}
  \item{curvecol}{color of the curve interpolating the local centers of mass (this is the "local principal curve"!).}
  \item{splinecol}{color of the spline representation of the local principal curve.} 
  \item{projectcol}{color of projections onto the spline representation of the local principal curve.}
  \item{startcol}{color of the plotted starting points.}
  \item{startpch}{plotting symbol for starting points; needs to be either a
  single symbol, or a vector of symbols of the same length as the number
of starting points.}
  \item{\dots}{ further arguments passed to \code{plot} or \code{scatterplot3d}.  }
}
%\details{
% Plots a local principal curve  
%}
\value{
A 2D plot, 3D plot, or a pairs plot (depending on the type of object and the data dimension \eqn{d}.).

The most flexible plotting option is \code{masscol}. Depending on the
length of the specified vector, this will be interpreted differently. If
a scalar is provided, the corresponding color will be given to all centers of
mass (or cluster centers). For LPCs, if the length of the vector is larger than 1, then this option
will assign different colours to  different depths, or different branch
numbers, or to individual data points, depending on the length. The
default setting is assigning colours according to depth, in the order
red, blue, black. 
 
With increasing dimension \eqn{d}, less plotting options tend to be supported. The nicest plots are obtained  for \eqn{d=2} and \eqn{d=3}.
 }


\section{Warning}{
This function computes all missing information (if posssible), so computation will take the longer the less informative the given object is, and the more advanced aspects are asked to plot!
} 


\references{
Einbeck, J., Tutz, G.,  and Evers, L. (2005). Local principal curves. Statistics and Computing 15, 301-313.

Einbeck, J., Evers, L. & Hinchliff, K. (2010):  Data compression and regression based on local principal curves. In A. Fink, B. Lausen, W. Seidel, and A. Ultsch (Eds),  Advances in Data Analysis, Data Handling, and Business Intelligence, Heidelberg, pp. 701--712, Springer. 
}



\author{ JE }
\seealso{ 
\code{\link{lpc}}, \code{\link{lpc.spline}}  , \code{\link{ms}}  
}

\examples{
data(calspeedflow)
lpc1 <- lpc(calspeedflow[,3:4])
plot(lpc1, type=c("spline","project"), lwd=2)
ms1<- ms(calspeedflow[,3:4], plot=FALSE)
plot(ms1, masscol=1)
plot(ms1, curvecol="grey30")

data(mussels, package="dr")
ms2 <- ms(mussels[,-3], scaled=1, h=0.1, plot=FALSE)
plot(ms2, datpch=20, masspch=24)
}
%\keyword{ smooth}% at least one, from doc/KEYWORDS
\keyword{multivariate}% __ONLY ONE__ keyword per line
