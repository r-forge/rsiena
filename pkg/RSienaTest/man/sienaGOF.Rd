\name{sienaGOF}
\alias{sienaGOF}
\alias{evaluateAuxiliaryStatistics}
\alias{evaluateAuxiliaryStatistics.simulated}
\alias{evaluateAuxiliaryStatistics.observed}
\alias{print.sienaGOF}
\alias{plot.sienaGOF}
\title{Functions to assess goodness of fit for SAOMs}
\description{
 A series of functions which start with a \code{sienaFit} object, a
 \code{siena} data object, and a function for calculating auxiliary
 statistics on graphs. 
 This function uses the \code{snopgof} package to evaluate goodness
 of fit from the resulting vectors of statistics.
 }

\usage{
sienaGOF(obsList, simList)
evaluateAuxiliaryStatistics.observed(sienaDataObject, varName, wave, auxiliaryFunction, verbose=FALSE)
evaluateAuxiliaryStatistics.simulated(sienaFitObject, groupName, varName, wave, auxiliaryFunction, verbose=FALSE)
# Alias for the above two functions
evaluateAuxiliaryStatistics(sienaObject, ...)
}
\arguments{
  \item{obsList}{The results from a call to \code{evaluateAuxiliaryStatistics.observed}.}
  \item{simList}{The results from a call to \code{evaluateAuxiliaryStatistics.simulated}.}
  \item{sienaObject}{Either a \code{sienaDataObject} or a \code{sienaFitObject}; alias will determine which and call
  the appropriate function.}
  \item{sienaDataObject}{Results from a call to \code{sienaDataCreate}.}
  \item{sienaFitObject}{Results from a call to \code{siena07}.}
  \item{groupName}{The name of the sienaGroup to be used.}
  \item{varName}{The value of the variable to be used, whether a network, behavior, etc.}
  \item{wave}{Wave(s) to be used for calculating the statistics.}
  \item{auxiliaryFunction}{Function to be used to calculate statistics, e.g. from the \code{sna} package.}
  \item{verbose}{Whether to print intermediate results. Calculations can take some time, and user feedback may be useful.}
  \item{...}{Other arguments to be used by the appropriate subfunction which will be called by the alias.}
}
\details{
 This function is currently under construction. Contact the author for working papers on how it works.
 In general, it tests goodness of fit of statistics that should not be explicitly fit. If the p value is
 below the desired false positive rate, there is sufficient evidence to reject the goodness of fit of the model
 with respect to the implied statistic.
}
\value{
  \code{sienaGOF} Returns a list:
  \item{ObservedStatistics }{ The values of the calculated statistics for the observations. }
  \item{SimulatedStatistics }{ A preprocess object (see \code{snopgof} package) calculated on the statistics for the simulations. }
  \item{Results }{ The p values of the \code{snopgof} test. }
}

\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}
  for general information on RSiena.

Lospinoso, J.A. and Snijders, T.A.B, "A Non-parametric Test for
Goodness of Fit." Working Paper
}
\author{Josh Lospinoso}
\seealso{\code{\link{siena07}}, \code{\link{sienaTimeTest}},
  \code{snopgof}}
\examples{
\dontrun{
# Fifty iterations being used, but we recommend using many more (perhaps 1000 or 1500)
mymodel <- sienaModelCreate(fn=simstats0c, nsub=4, n3=50)
mynet1 <- sienaNet(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTrip, balance)
ans <- siena07(mymodel, data=mydata, effects=myeff, returnDeps=T, batch=TRUE)

library(sna)
# install.packages("sna")
library(network)
# install.packages("network")
triadCensus <- function (x) {
	triad.census(network(x))[-1]
}

sim <- evaluateAuxiliaryStatistics(ans, "Data1", "mynet1", 1:2, triadCensus, verbose=T)
obs <- evaluateAuxiliaryStatistics(mydata, "mynet1", 1:2, triadCensus, verbose=T)
(res <- sienaGOF(obs, sim))

# And plots if desired
plot(res, wave=1)
plot(res, wave=2)
}
}
\keyword{models}
