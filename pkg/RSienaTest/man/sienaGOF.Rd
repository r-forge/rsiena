\name{sienaGOF}
\alias{sienaGOF}
\alias{print.sienaGOF}
\alias{plot.sienaGOF}

\title{Functions to assess goodness of fit for SAOMs}
\description{
 Goodness of fit functionality requiring a \code{sienaFit} object, a
 \code{siena} data object, and a function for calculating auxiliary
 statistics on graphs. A Monte Carlo test based on the Mahalanobis
 distance based test is used to calculate frequentist p-values and
 plotting functions can be used to diagnose bad fit. Alternately, a
 Kolmogorov-Smirnov test statistic can be used for cumulative
 auxiliary statistics.
 }

\usage{
sienaGOF(sienaDataObject,
		sienaFitObject, groupName, varName,  auxiliaryFunction, wave=NULL,
		verbose=FALSE, join=TRUE, expectationFunction=mean,
		twoTailed=FALSE, cumulative=FALSE, \dots)
\method{plot}{sienaGOF}(x,  standardize=NA, violin=TRUE, 
			ylim=NULL, xlim=NULL, 
			xlab=NULL, ylab=NULL, key=NULL, 
			perc=.05, wave=1, main=NULL, 
			image=FALSE, imageThreshold=.3, \dots)
\method{print}{sienaGOF}(x, \dots)
}
\arguments{
  \item{sienaDataObject}{Results from a call to \code{sienaDataCreate}.}
  \item{sienaFitObject}{Results from a call to \code{siena07}.}
  \item{groupName}{The name of the sienaGroup to be used.}
  \item{varName}{The value of the variable to be used, whether a 
  network, behavior, etc.}
  \item{join}{Boolean: should sienaGOF do tests on all of the waves 
  individually, or sum across waves?}
  \item{auxiliaryFunction}{Function to be used to calculate statistics, 
  e.g. from the \code{sna} package.
  	The function must accept two arguments of the order (dependentVariable, 
  	missingData), and return a
  	numeric vector of constant length.}
  \item{verbose}{Whether to print intermediate results. Calculations can 
  take some time, and user feedback may be useful.}
  \item{expectationFunction}{The function to be used for centering the 
  Mahalanobis distances. 
  Not used if \code{cumulative=TRUE}.}
  \item{cumulative}{If the auxiliary function returns cumulative values, 
  this option will use a Kolmogorov-Smirnov
  test statistic instead of the Mahalanobis distance.}
  \item{twoTailed}{Whether to use two tails for calculating p values on the
   Monte Carlo test. Recommended for 
  advanced users only, as it is probably only applicable in rare cases.}
  \item{x}{ Object from a call to sienaGOF. }
  \item{standardize}{ One of three levels of standardization to be performed
  on the auxiliary statistics (for layout):
  3 -- center by median and scale;
  2 -- map to [0,1] interval;
  1 -- center by mean;
  0 -- perform no manipulations.
  Default for non-cumulative is 3, default for cumulative is 0.}
  \item{violin}{ Use violin plots (vs. box plots)? }
  \item{ylim}{ Limits on x axis. For plotting sections of the test. }
  \item{xlim}{ Limits on y axis. For plotting sections of the test. }
  \item{xlab}{ Label for the x axis. }
  \item{ylab}{ Label for the y axis. }
  \item{key}{ Keys for the auxiliary statistic levels. }
  \item{perc}{ P value for the confidence bands. }
  \item{main}{ Main title. }
  \item{image}{ Plot the image matrix of dyadic errors only.
  This option will negate most of the others. }
  \item{imageThreshold}{ Threshold for the ratio of errors required
  to plot a false dyad in the image.}
  \item{wave}{Wave(s) to be used.}
  \item{\dots}{Other arguments.}
}
\details{
 This function is used to assess the goodness of fit of a stochastic actor 
 oriented model
 for an arbitrarily defined auxiliary statistic. These statistics should be
  chosen to
 represent features of the network that are not explicitly fit but can be 
 considered
 important properties that the model at hand should represent well. Some 
 examples are:
 
 -Geodesic distances
 
 -Triad census
 
 -Outdegree distribution
 
 -Indegree distribution
 
 -Behavioral distribution
 
 -Edgewise homophily counts
 
 -Edgewise shared partners

The function is written so that the user can easily define a function to 
capture some
other relevant aspects of the network, behaviors, etc.

We recommend the following heuristic approach to fitting your model:

1. Check convergence of your estimation.

2. Assess time heterogeneity and either modify the base effects or include
 time dummy terms.

3. Assess goodness of fit (i.e. using join=TRUE) on auxiliary statistics,
 and refine the model.
}
\value{
  \code{sienaGOF} Returns a list of class \code{sienaGofTest}, which
   contains the following:
  \item{SimulatedMhd}{ Mahalanobis distances for the simulations. }
  \item{ObservedMhd}{ Mahalanobis distances for the observation. }
  \item{TwoTailed}{ Whether the p value corresponds to a one or two 
  tailed Monte Carlo test. }
  \item{Simulations}{ A copy of the statistics from the auxiliaryFunction 
  for simulations. }
  \item{Observations}{ A copy of the statistics from the auxiliaryFunction 
  for observations. }
  \item{Rank}{ Rank of the covariance matrix of the simulated auxiliary 
  statistics. }
}

\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}
  for general information on RSiena and \url{http://www.stats.ox.ac.uk/~lospinos/}
 for information on the following presentation and working paper:
 
Lospinoso, J.A. and Snijders, T.A.B, "Goodness of fit for
Stochastic Actor Oriented Models." Presentation given at Sunbelt XXXI,
 St. Pete's Beach, Fl. 2011.

Lospinoso, J.A. and Snijders, T.A.B, "Goodness of fit for
Stochastic Actor Oriented Models." Working Paper.
}
\author{Josh Lospinoso}
\seealso{\code{\link{siena07}}, \code{\link{sienaTimeTest}} }
\examples{
# Fifty iterations being used, but we recommend using
# many more (perhaps 1000 or 1500)
mymodel <- sienaModelCreate(fn=simstats0c, nsub=4, n3=50)
mynet1 <- sienaNet(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, gwespFB, balance)
# The 'parm' column only accepts integers, so the convention is
# to specify alpha * 100, i.e. in hundredths:
myeff[myeff$shortName=='gwespFB', 'parm'] <- 30
ans <- siena07(mymodel, data=mydata, effects=myeff, returnDeps=TRUE, batch=TRUE)

require(sna)
require(network)

missingsToZero <- function(x, missingData) {
	x[is.na(x)] <- 0
	1*( (x==1 & !missingData) )
}

triadCensus <- function (x, missingData) {
	# Zero out the missing data
	zeroed <- missingsToZero(x, missingData)
	triad.census(as.matrix(zeroed))
}

res <- sienaGOF(mydata, ans,"Data1", "mynet1", triadCensus, verbose=TRUE,
 join=TRUE)

# And plots if desired
plot(res, key=c("003", "012","102","021D","021U","021C",
	"111D","111U","030T","030C","201",
	"120D","120U","120C","210","300"))

# We can also do cumulative tests:
outdegreeDistribution <- function(x, missingData) {
	zeroed <- missingsToZero(x, missingData)
	outdegrees <- apply(zeroed,1,sum)
	# Returns the cumulative distribution of outdegrees from levels 0 to 7
	cumulativeOutdegrees = sapply(0:7, function(i)sum(outdegrees<=i))
}
res <- sienaGOF(mydata, ans,"Data1", "mynet1", outdegreeDistribution, 
join=TRUE, cumulative=TRUE)
plot(res,key=paste("OD<=",0:7,sep=""))

# Finally, see the image matrix for the simulations
plot(res, image=TRUE)
}
\keyword{models}
