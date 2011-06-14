\name{calspeedflow}
\Rdversion{1.1}
\alias{calspeedflow}
\docType{data}
\title{
Speed-flow data from Calfornia. 
}
\description{
A `fundamental diagram' with observations of speed and flow recorded from 9th of July 2007, 9am, to 10th of July 2007, 10pm, on Line 5 of the Californian Freeway SR57-N, VDS number 1202263. The data were originally measured in intervals of thirty seconds, and then aggregated over intervals of 5 minutes length. 
}
\usage{data(calspeedflow)}
\format{
  A data frame with 444 observations on the following 4 variables.
  \describe{
    \item{\code{Date}}{a factor with levels \code{07/09/2007}... \code{07/10/2007}.}
    \item{\code{Timestamp}}{a factor with a timestamps in intervals of five minutes.}
    \item{\code{Lane5Flow}}{a numeric vector of vehicle flow in vehicles per 5 minutes.}
    \item{\code{Lane5Speed}}{a numeric vector of vehicle speed in miles per hour.}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
Retrieved from PeMS.
}
\references{
Einbeck, J.,  and Dwyer, J. (2011). Using principal curves to analyze traffic patterns on freeways. Transportmetrica 7, 229-246.
}
\examples{
data(calspeedflow)
plot(calspeedflow[,3:4])
}
\keyword{datasets}
