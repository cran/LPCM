\name{gaia}
\alias{gaia}
\docType{data}
\title{
Gaia data
}
\description{
(Simulated) spectral decomposition of stellar objects, generated in the framework of the Gaia project.
}
\usage{data(gaia)}
\format{
  A data frame with 8286 observations on the following 22 variables.
  \describe{
    \item{\code{ID}}{ID of the object}
    \item{\code{metallicity}}{metallicity (abundance); that is proportion of matter  other than hydrogen and helium relative to that of the sun.}
    \item{\code{gravity}}{the surface gravity; that is acceleration due to gravity at the surface of the star. }
    \item{\code{temperature}}{the `effective' temperature (K); that is the temperature of the observable part of the stellar atmosphere.}
    \item{\code{band1}}{photon counts in band 1}
    \item{\code{band2}}{photon counts in band 2}
    \item{\code{band3}}{photon counts in band 3}
    \item{\code{band4}}{photon counts in band 4}
    \item{\code{band5}}{photon counts in band 5}
    \item{\code{band6}}{photon counts in band 6}
    \item{\code{band7}}{photon counts in band 7}
    \item{\code{band8}}{photon counts in band 8}
    \item{\code{band9}}{photon counts in band 9}
    \item{\code{band10}}{photon counts in band 10}
    \item{\code{band11}}{photon counts in band 11}
    \item{\code{band12}}{photon counts in band 12}
    \item{\code{band13}}{photon counts in band 13}
    \item{\code{band14}}{photon counts in band 14}
    \item{\code{band15}}{photon counts in band 15}
    \item{\code{band16}}{photon counts in band 16}
  }
}
\details{
Gaia is an astrophysics mission of the European Space Agency (ESA) which will
undertake a detailed survey of over 10^9 stars in our Galaxy and
extragalactic objects.  An important part of the scientific analysis of these
data is the classification of all the objects as well as the estimation of
stellar astrophysical parameters (effective stellar temperature, surface gravity, metallicity). This will be done on the basis of
high-dimensional spectroscopic and astrometric data such as those ones given here. 

More precisely, the  spectral data come in form of photon counts ("fluxes") observed in (originally) 96 wavelength intervals ("bands"), see Bailer-Jones (2010) for more details. The data given here are a 16-dimensional subset created by binning/selecting from the 96 bands.  The counts given here are standardized, i.e. they are divided by the total number of incoming photons over all filters  (in other words, they add up to 1).  Note that these data are simulated using computer models. The satellite which will collect the actual data will be launched in 2012.

The 16-d spectral data  have been used in Einbeck, Evers and Bailer-Jones (2008) as well as Einbeck, Evers and Powell (2010) in order to predict the stellar temperature.

}
\source{
Coryn Bailer-Jones (MPIA Heidelberg).
}
\references{
Bailer-Jones, C.A.L. (2010). The ILIUM forward modelling algorithm for multivariate parameter estimation and its application to derive stellar parameters from Gaia spectrophotometry, Monthly Notices of the Royal Astronomical Society, vol. 403, pp. 96-116.

Einbeck, J., Evers, L., and Bailer-Jones, C.A.L. (2008).  Representing complex data using localized principal components with application to astronomical data. In: Gorban, A, Kegl, B, Wunsch, D, & Zinovyev, A:  Principal Manifolds for Data Visualization and Dimension Reduction; Lecture Notes in Computational Science and Engineering  58, 180-204, ISSN/ISBN: 978-3-540-73749-0. 

Einbeck, J., Evers, L., and Powell, B. (2010): Data compression and regression through local principal curves and surfaces, International Journal of Neural Systems, 20, 177-192.
}

\examples{
data(gaia)
s <- sample(nrow(gaia),200)
library(lattice)
splom(gaia[s,5:20], cex=0.3,pscales=0)

gaia.pc <-  princomp(gaia[s,5:20])
temp <- gaia$temperature
tempcol     <- (temp[s]- min(temp[s]))/max(temp[s]- min(temp[s]))
library(scatterplot3d)
scatterplot3d(gaia.pc$scores[,c(2,1,3)], pch="+",
     color=rgb(sqrt(tempcol),0,1-sqrt(tempcol)))
     # This is a 3D scatterplot of the first three principal component scores;
     # with higher stellar temperatures shaded in red colour.
}

\keyword{datasets}
