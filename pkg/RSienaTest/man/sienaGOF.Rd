\name{sienaGOF}
\alias{sienaGOF}
\alias{plot.sienaGOF}
\title{Functions to assess goodness of fit for SAOMs}
\description{
 Goodness of fit functionality requiring a \code{sienaFit} object, a
 \code{siena} data object, and a function for calculating auxiliary
 statistics on graphs. A Monte Carlo test based on the Mahalanobis
 distance based test is used to calculate frequentist p-values and
 plotting functions can be used to diagnose bad fit. There are basic
 functions for calculating auxiliary statistics available out of the
 box, and the user is also permitted to create custom functions.
 }

\usage{
sienaGOF(sienaFitObject, auxiliaryFunction, groupName=NULL, 
		varName=NULL, wave=NULL, verbose=FALSE, join=TRUE, twoTailed=FALSE, 
		cluster=NULL, robust=FALSE, \dots)
\method{plot}{sienaGOF}(x,  center=FALSE, scale=FALSE, violin=TRUE, 
		key=NULL, perc=.05, wave=1, main=main, ylab=ylab, \dots)
}
\arguments{
  \item{sienaFitObject}{Results from a call to \code{siena07}.}
  \item{groupName}{The name of the sienaGroup to be used.}
  \item{varName}{The value of the variable to be used, whether a 
  network, behavior, etc.}
  \item{join}{Boolean: should sienaGOF do tests on all of the waves 
  individually, or sum across waves?}
  \item{auxiliaryFunction}{Function to be used to calculate statistics, 
  e.g. one implementing the \code{sna} or \code{igraph} packages.
  	There functions \code{IndegreeDistribution}, \code{OutdegreeDistribution},
  	\code{GeodesicDistribution}, \code{TriadCensus}, and \code{KnnDistribution}
  	are already available with \code{RSiena}. The user can write custom functions
  	fairly easily, however. See \code{Examples} for more information.}
  \item{verbose}{Whether to print intermediate results. Calculations can 
  take some time, and user feedback may be useful.}
  \item{twoTailed}{Whether to use two tails for calculating p values on the
   Monte Carlo test. Recommended for 
  advanced users only, as it is probably only applicable in rare cases.}
  \item{x}{ Object from a call to sienaGOF. }
  \item{ylab}{ The y-axis label for the plot. }
  \item{robust}{ Whether to use robust estimation of the covariance matrix. }
  \item{center}{ Whether to center the statistics by median during plotting.}
  \item{scale}{ Whether to scale the statistics by range during plotting.}
  \item{violin}{ Use violin plots (vs. box plots only)? }
  \item{key}{ Keys for the auxiliary statistic levels in the plot. }
  \item{perc}{ P value for the confidence bands (two sided). }
  \item{main}{ Main title. }
  \item{wave}{Wave(s) to be used.}
  \item{cluster}{Optionally, a \code{snow} cluster to execute the auxiliary
  function calculations on.}
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
 
 The \code{print} function will display some useful information to help with
 model selection if some effects are set to FIX and TEST on the effects
 object. A rough estimator for the Mahalanobis distance at each proposed
 specification is given in the output. There are three headings for various
 estimators of this specification:
 
 --MM(Full) gives the one step estimator
 which will correspond with the estimator available on the \code{siena07}
 object
 
 --MM(Par.) only adds the candidate effect coordinate to the current estimate
 
 --GMM, for "Generalized Method of moments", gives the approximate global
 minimum Mahalanobis distance that could be achieved under the specification.
 
 For these specifications, Mahalanobis distances as well as parameter estimates
 are given to help guide model selection. See the manual, or the citations for
 more information.

These functions are pre-fabricated for ease of use, and can be passed in as the
\code{auxiliaryFunction} with no extra effort:
 
\code{ GeodesicDistribution(\dots)}

\code{ IndegreeDistribution(\dots)}

\code{ KnnDistribution(\dots)}

\code{ OutdegreeDistribution(\dots)}

\code{ TriadCensus(\dots)}

These functions are pre-fabricated for those users wanting to create their own
custom auxiliary functions. They will help with extracting the appropriate
ingredients from the observations and simulations:

\code{igraphEdgelistExtraction(\dots)}

\code{snaEdgelistExtraction(\dots)}

\code{snaSociomatrixExtraction(\dots)}

\code{sparseMatrixExtraction(\dots)}

See the manual for instructions on how to implement these extractor functions.

}
\value{
  \code{sienaGOF} Returns a list of class \code{sienaGofTest}, which
   includes the following:
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
\dontrun{
	mymodel <- sienaModelCreate(fn=simstats0c, nsub=4, n3=1000)
	mynet1 <- sienaNet(array(c(tmp3,tmp4), dim=c(32, 32, 2)))
	mydata <- sienaDataCreate(mynet1)
	myeff <- getEffects(mydata)
	myeff <- includeEffects(myeff, transTrip, inAct, cycle3, balance)
	myeff[myeff$shortName=='inAct', c('fix', 'test')] <-TRUE
	myeff[myeff$shortName=='cycle3', c('fix', 'test')] <-TRUE
	myeff[myeff$shortName=='balance', c('fix', 'test')] <-TRUE
	ans <- siena07(mymodel, data=mydata, effects=myeff, returnDeps=TRUE, batch=TRUE)
	( res <- sienaGOF(ans, KnnDistribution) )
	plot(res)
	
	## We can also make a custom function:
	#
	CustomTriadCensus <- function(...) {
		require(sna)
		x <- snaEdgelistExtraction(...)
		triad.census(x)[6:16]
	}
	# Execute the test:
	( res <- sienaGOF(ans, CustomTriadCensus) )	
}
}
\keyword{models}
