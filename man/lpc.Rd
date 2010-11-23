\name{lpc}
\Rdversion{1.1}
\alias{lpc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Local principal curves
}
\description{
This is the main function which computes the actual local principal curve, i.e. a sequence of local centers of mass.
}
\usage{
lpc(X, h, t0 = mean(h),  x0,  mult = 1, way = "two",  scaled = TRUE,
      weights=1, pen = 2, depth = 1, control=lpc.control())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
data matrix with \eqn{N} lines (observations) and \eqn{d} columns (variables). Note that \eqn{d>2} restricts \code{depth} to 1.
}
  \item{h}{
bandwidth. May be either specified as a single number, then the same bandwidth is used in
    all dimensions, or as a d-dimensional bandwidth vector. The default
    setting is 10 percent of the  range in each direction. If \code{scaled =TRUE} then the bandwidth has to be
    specified in fractions of the data range, e.g. h= c(0.2,0.1), rather than
    absolute values. 
}
  \item{t0}{scalar step length. Default setting is \code{t0=h}, if \code{h} is a scalar, and \code{t0=mean(h)}, if \code{h} is a vector.
}

  \item{x0}{
    optionally, one can set one or more  starting points here. This can
    be done in form of a matrix, where each row corresponds to a
    starting point, or in form of a vector, where starting points are
    read in consecutive order from the entries of the vector. 
    The starting point has always to be specified on the original data scale, even if \code{scaled=TRUE}. 
   }
    \item{mult}{
numerical value which specifies the number of starting points.  If the
    number  given here is larger than the number of starting points
    provided at \code{x0},  then the missing points will be set at
    random. For example, if \eqn{d=2}, \code{mult=3}, and
    \code{x0=c(58.5, 17.8, 80,20)},  then one gets the starting points (58.5, 17.8), (80,20), and a randomly
     chosen  third one.  Another example for such a situation  is \code{mult=1}, with
    \code{x0} being unspecified.
     }
  \item{way}{"one": go only in direction of the first local eigenvector,
     "back": go only in opposite direction,
 	 "two": go from starting point in both directions.
       }
\item{scaled}{
if TRUE, scales each variable by dividing through its range (see also
the Notes section below).
}
 \item{weights}{a vector of observation weights (can also be used to exclude
 individual observations from the computation by setting their weight to
 zero.)} 
  \item{pen}{
  power used for angle penalization (see [1]). If set to 0, the angle
 penalization is switched off.
}
  \item{depth}{
maximum depth of branches (phimax in [2]), restricted to the values 1,2 or 3, where 2 and 3 are only allowed for \eqn{d=2}.
   }
\item{control}{
Additional parameters steering partcularly the boundary and convergence
 behavior of the fitted curve. See \code{\link{lpc.control}}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A list of items:
  \item{LPC}{
        The coordinates of the local centers of mass of the fitted
	principal curve.
	}
\item{Parametrization}{Curve parameters and branch labels for
 each local center of mass.} 
\item{h}{
   The bandwidth used for the curve estimation.
 }
 \item{to}{
   The constant \eqn{t_0}{t0} used for the curve estimation.
 }
 \item{starting.points}{
   The coordinates of the starting point(s) used.
 }  
 \item{data}{
    The data frame used for curve estimation.
  }
 \item{scaled}{Logical.} 
 \item{weights}{The vector of weights used for curve estimation.}
 \item{control}{The settings used in \code{lpc.control()}} 
 \item{Misc}{ Miscellannea.}
 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
[1] Einbeck, J., Tutz, G., & Evers, L. (2005). Local principal curves. Statistics and Computing 15, 301-313.

[2] Einbeck, J., Tutz, G., & Evers, L. (2005): Exploring Multivariate Data Structures with Local Principal Curves. In: Weihs, C. and Gaul, W. (Eds.): Classification - The Ubiquitous Challenge. Springer, Heidelberg, pages 256-263. 
}

\author{
J. Einbeck and L. Evers, with contributions by J.
Dwyer and M. Zayed.
}
\note{
All values provided in the output refer to the scaled data, if
\code{scaled=TRUE}.  Use \code{\link{lpc.unscale}} to convert the
results back to the original data scale. 
  
The option \code{scaled=TRUE} scales the data by dividing each variable through their
range. This differs from the usual scaling through
the standard deviation as common for  PCA, but we found the
algorithm and the default bandwidth selection to work more reliably
this way. If you wish to scale by the standard deviation, please do that
by feeding the scaled data directly into the \code{lpc}
function, i.e.

\code{lpc(sweep(data, 2,sd(data), "/"),h, t0, ..., scaled=FALSE, ...)}.

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

data(calspeedflow)
lpc1 <- lpc(calspeedflow[,3:4])
plot(lpc1)

data(mussels, package="dr")
lpc2 <- lpc(mussels[,-3], x0=as.numeric(mussels[49,-3]), scaled=FALSE)
pairs(lpc2$LPC, panel="lines", labels= c("H","L","S","W") )

data(gaia)
s <- sample(nrow(gaia),200)
gaia.pc <-  princomp(gaia[s,5:20])
lpc3 <- lpc(gaia.pc$scores[,c(2,1,3)],scaled=FALSE)
plot(lpc3,curvecol=2)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smooth}
\keyword{multivariate}% __ONLY ONE__ keyword per line
