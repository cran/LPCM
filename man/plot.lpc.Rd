\name{plot.lpc}
\alias{plot.lpc}
\alias{plot.lpc.spline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plotting local principal curves}
\description{
  Takes an object of class \code{lpc} or \code{lpc.spline} and plots any subset of the following components of the local principal curve:  Centers of mass; the curve connecting the local centers of mass; the cubic spline representation of the curve; the projections onto the curve; the starting points.
}
\usage{
\method{plot}{lpc}(x, type, unscale = TRUE, lwd = 1, datcol = "grey60", 
    datpch = 21, masscol = 2, masspch = 15, curvecol = 1, splinecol = 3, 
    projectcol = 4, startcol = 1:3,  startpch=NULL,...)  
\method{plot}{lpc.spline}(x, type, unscale = TRUE, lwd = 1, datcol = "grey60", 
    datpch = 21, masscol = 2, masspch = 15, curvecol = 1, splinecol = 3, 
    projectcol = 4, startcol = 1:3,  startpch=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{lpc} or \code{lpc.spline}.   }
  \item{type}{a vector of type \code{c("mass", "spline",...)} with possible entries \code{mass, curve, spline, project, start.}}
  \item{unscale}{if TRUE, then data (and all fitted componens) are scaled back to their original scale; otherwise the scaled data are plotted (only relevant if \code{scaled=TRUE}  in the fitted object). }
  \item{lwd}{width of curves.}
  \item{datcol}{color of data points.}
  \item{datpch}{plotting symbol for data points.}
  \item{masscol}{color of centers of mass.}
  \item{masspch}{plotting symbol for centers of mass.}
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
 A 2D plot, 3D plot, or a pairs plot (depending on the data dimension \eqn{d}.).

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
\code{\link{lpc}}, \code{\link{lpc.spline}}  
}

\examples{
data(calspeedflow)
lpc1 <- lpc(calspeedflow[,3:4])
plot(lpc1, type=c("spline","project"), lwd=2)
}

%\keyword{ smooth}% at least one, from doc/KEYWORDS
\keyword{multivariate}% __ONLY ONE__ keyword per line
