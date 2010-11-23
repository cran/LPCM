\name{gvessel}
\alias{gvessel}
\docType{data}
\title{
 North Atlantic Water Temperature Data.
}
\description{
These are observations taken over nine days in May 2000 by the German vessel Gauss in the North Atlantic. 
}
\usage{data(gvessel)}
\format{
  A data frame with 643 observations on the following 7 variables. 
  \describe{
    \item{\code{day2g}}{an integer for the day at which the measurement was taken.}
    \item{\code{salg}}{a numeric vector with measurements of salinity according to the PSS (Practical Salinity Scale).}
    \item{\code{tempg}}{a numeric vector with measurements of water temperature in degrees Celsius.}
    \item{\code{depthg}}{a numeric vector with the water depths (in meters) at which the measurements were taken.}
    \item{\code{oxyg}}{a numeric vector with measurements of oxygen content (mm per litre of water)}
    \item{\code{longg}}{longitude}
    \item{\code{latg}}{latitude}
  }
}
%\details{
%}
\source{
Retrieved by B. Powell from the World Ocean Database, \url{http://www.nodc.noaa.gov/OC5/WOD/pr_wod.html}.
}
\references{
Einbeck, J., Evers, L., and Powell, B. (2010): Data compression and regression through local principal curves and surfaces, International Journal of Neural Systems, 20, 177-192.
}
\examples{
data(gvessel)
pairs(gvessel[,c(3,2,4,5)])
tcol <- (gvessel$tempg- min(gvessel$tempg))/(max(gvessel$tempg)- min(gvessel$tempg))
require(scatterplot3d)
scatterplot3d(gvessel[,2],gvessel[,4],gvessel[,5], color=rgb(tcol,0,1-tcol))
}
\keyword{datasets}
