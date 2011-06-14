\name{lpc.control}
\Rdversion{1.1}
\alias{lpc.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Auxiliary parameters for controlling local principal curves.
}
\description{
This function bundles parameters controlling mainly the convergence-, boundary-,
and stopping-behaviour of the local principal curve. It will be used
only inside the \code{lpc()} function argument. 
}
\usage{
lpc.control(iter =100, boundary = 0.005, convergence.at =0.00001,
     pruning.thresh=0.0, rho0=0.4, cross=TRUE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iter}{
Maximum number of iterations on either side of the starting point within each branch.
}
  \item{boundary}{
This boundary correction [2] reduces the bandwidth adaptively once the
relative difference of parameter values between two centers of mass
falls below the given threshold. This measure delays convergence and
enables the curve to proceed further into the end points. If set to 0,
this boundary correction is switched off. 
}
\item{convergence.at}{
This forces the curve to stop if  the
relative difference of parameter values between two centers of mass
falls below the given threshold.  If set to 0, then the curve will
always stop after exactly \code{iter} iterations.
}
  \item{pruning.thresh}{
Prunes branches corresponding to higher-depth starting points if
their density estimate falls below this threshold.  Typically, a value between 0.0
and 1.0. The setting 0.0 means no pruning.
}
  \item{rho0}{
A numerical value which  steers the birth process of higher-depth starting
points. Usually, between 0.3 and 0.4 (see reference [1]). 
}
  \item{cross}{Logical parameter. If TRUE, curves are stopped when they
    come too close to an existing branch. Used in the self-coverage function. 
       }
     }


     
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list of the six specified parameters, which can be read by the
\code{control} argument of the \code{lpc} function. 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
[1] Einbeck, J., Tutz, G. & Evers, L. (2005): Exploring Multivariate Data Structures with Local Principal Curves. In: Weihs, C. and Gaul, W. (Eds.): Classification - The Ubiquitous Challenge. Springer, Heidelberg, pages 256-263. 

[2] Einbeck, J. and Zayed, M. (2011). Some asymptotics for localized
principal components and curves. Working paper, Durham University. Unpublished.
}

\author{
JE
}


\examples{
data(calspeedflow)
fit1 <- lpc(calspeedflow[,c(3,4)], x0=c(50,60),scaled=TRUE,
   control=lpc.control(iter=20, boundary=0))
plot(fit1, type=c("curve","start","mass"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{smooth}
%\keyword{multivariate}% __ONLY ONE__ keyword per line
