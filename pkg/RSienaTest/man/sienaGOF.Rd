\name{sienaGOF}
\alias{sienaGOF}
\alias{evaluateAuxiliaryStatistics}
\alias{evaluateAuxiliaryStatistics.simulated}
\alias{evaluateAuxiliaryStatistics.observed}
\alias{print.sienaGOF}
\alias{plot.sienaGOF}
\alias{sienaGOFPreprocess}
\alias{sienaGOFProcess}

\title{Functions to assess goodness of fit for SAOMs}
\description{
 A series of functions which start with a \code{sienaFit} object, a
 \code{siena} data object, and a function for calculating auxiliary
 statistics on graphs. A Monte Carlo test based on the Mahalanobis
 distance based test is used to calculate frequentist p-values and
 plotting functions can be used to diagnose bad fit.
 }

\usage{
sienaGOF(obsList, simList)
evaluateAuxiliaryStatistics.observed(sienaDataObject, varName, auxiliaryFunction, wave=NULL, 
verbose=FALSE, join=TRUE)
evaluateAuxiliaryStatistics.simulated(sienaFitObject, groupName, varName, 
auxiliaryFunction, wave=NULL, verbose=FALSE, join=TRUE)
# Alias for the above two functions, for user convenience
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
  \item{join}{Boolean: should sienaGOF do tests on all of the waves individually, or sum across waves?}
  \item{wave}{Wave(s) to be used for calculating the statistics. Ignored if join=TRUE}
  \item{auxiliaryFunction}{Function to be used to calculate statistics, e.g. from the \code{sna} package.}
  \item{verbose}{Whether to print intermediate results. Calculations can take some time, and user feedback may be useful.}
  \item{...}{Other arguments to be used by the appropriate subfunction which will be called by the alias.}
}
\details{
 This function is used to assess the goodness of fit of a stochastic actor oriented model
 for an arbitrarily defined auxiliary statistic. These statistics should be chosen to
 represent features of the network that are not explicitly fit but can be considered
 important properties that the model at hand should represent well. Some examples are:
 
 -Geodesic distances
 
 -Triad census
 
 -Outdegree distribution
 
 -Indegree distribution
 
 -Behavioral distribution
 
 -Edgewise homophily counts
 
 -Edgewise shared partners

The function is written so that the user can easily define a function to capture some
other relevant aspects of the network, behaviors, etc.

We recommend the following heuristic approach to fitting your model:

1. Check convergence of your estimation.

2. Assess time heterogeneity and either modify the base effects or include time dummy terms.

3. Assess goodness of fit (i.e. using join=TRUE) on auxiliary statistics, and refine the model.
}
\value{
  \code{sienaGOF} Returns a list:
  \item{ObservedStatistics }{ The values of the calculated statistics for the observations. }
  \item{SimulatedStatistics }{ A preprocess object calculated on the statistics for the simulations. }
  \item{Results }{ A list including, among other things, p values of the Monte Carlo Mahalanobis distance test. }
}

\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}
  for general information on RSiena and \url{http://www.stats.ox.ac.uk/~lospinos/}
 for information on the following presentation and working paper:
 
Lospinoso, J.A. and Snijders, T.A.B, "Goodness of fit for
Stochastic Actor Oriented Models." Presentation given at Sunbelt XXXI, St. Pete's Beach, Fl. 2011.

Lospinoso, J.A. and Snijders, T.A.B, "Goodness of fit for
Stochastic Actor Oriented Models." Working Paper.
}
\author{Josh Lospinoso}
\seealso{\code{\link{siena07}}, \code{\link{sienaTimeTest}} }
\examples{
# Fifty iterations being used, but we recommend using many more (perhaps 1000 or 1500)
mymodel <- sienaModelCreate(fn=simstats0c, nsub=4, n3=50)
mynet1 <- sienaNet(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTrip, balance)
ans <- siena07(mymodel, data=mydata, effects=myeff, returnDeps=TRUE, batch=TRUE)

require(sna)
require(network)
triadCensus <- function (x) {
	triad.census(network(x))[-1]
}

sim <- evaluateAuxiliaryStatistics(ans, "Data1", "mynet1", triadCensus)
obs <- evaluateAuxiliaryStatistics(mydata, "mynet1", triadCensus)
(res <- sienaGOF(obs, sim))

# And plots if desired
triadCensus.names <- c("012","102","021D","021U","021C","111D","111U","030T","030C","201","120D","120U","120C","210","300")
plot(res, key=triadCensus.names)
}
\keyword{models}
