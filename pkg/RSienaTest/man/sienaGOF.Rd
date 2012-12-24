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
sienaGOF(sienaFitObject, auxiliaryFunction,
		wave=NULL, verbose=FALSE, join=TRUE, twoTailed=FALSE,
		cluster=NULL, robust=FALSE, \dots)
\method{plot}{sienaGOF}(x,  center=FALSE, scale=FALSE, violin=TRUE,
		key=NULL, perc=.05, wave=1, main=main, ylab=ylab, \dots)
}
\arguments{
  \item{sienaFitObject}{Results from a call to \code{siena07}.}
  \item{join}{Boolean: should sienaGOF do tests on all of the waves
  individually, or sum across waves?}
  \item{auxiliaryFunction}{Function to be used to calculate statistics,
  e.g. one implementing the \code{sna} or \code{igraph} packages.
  	See \code{Examples} for more information on the method signature. The method signature
    is \code{function(index, data, sims, wave)}, where the \code{index} is the index
    of the simulated network, or \code{NULL} if the observed variable is needed. \code{data}
    is the observed data object that you must extract the relevant variables from. \code{sims} is the
    list of simulations returned from \code{siena07}. \code{wave} is the index of the wave.}
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
 specification is given in the output. This estimator is called the
 modified model distance (MMD). See Lospinoso (2012) for more information.

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
  for general information on RSiena.:

Lospinoso, J.A.  and Snijders, T.A.B., "Goodness of fit for
Stochastic Actor Oriented Models." Presentation given at Sunbelt XXXI,
 St. Pete's Beach, Fl. 2011.

Lospinoso, J.A., "Statistical Models for Social
Network Dynamics." Ph.D. Thesis. University of Oxford: U.K.
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

  # An example function:
  sparseMatrixExtraction <- function (i, data, sims, wave, groupName="Data1", varName="mynet1") {
    require(Matrix)
    dimsOfDepVar<-
      attr(data[[groupName]]$depvars[[varName]],
           "netdims")
    missing <- Matrix(is.na(data[[groupName]]$depvars[[varName]][,,wave+1])*1)
    if (is.null(i)) {
      # sienaGOF wants the observation:
      returnValue <- Matrix(data[[groupName]]$depvars[[varName]][,,wave+1])
      returnValue[is.na(returnValue)] <- 0
    }
    else
    {
      #sienaGOF wants the i-th simulation:
      returnValue <- sparseMatrix(
        sims[[i]][[groupName]][[varName]][[wave]][,1],
        sims[[i]][[groupName]][[varName]][[wave]][,2],
        x=sims[[i]][[groupName]][[varName]][[wave]][,3],
        dims=dimsOfDepVar[1:2] )
    }
    ## Zero missings:
    1*((returnValue - missing) > 0)
  }

  IndegreeDistribution <- function (i, data, sims, wave,
                                    levels=0:6, extractor=sparseMatrixExtraction) {
    x <- extractor(i, data, sims, wave)
    a <- apply(x, 2, sum)
    sapply(levels, function(i){ sum(a<=i) })
  }

  ( res <- sienaGOF(ans, IndegreeDistribution) )
	plot(res)

  # Here are some useful functions for building your own auxiliary statistics:
  snaEdgelistExtraction <- function (i, data, sims, groupName, varName, wave) {
    require(sna)
    returnValue <- snaSociomatrixExtraction(i, data, sims, groupName, varName, wave)
    as.edgelist.sna(returnValue)
  }

  igraphEdgelistExtraction <- function (i, data, sims, groupName, varName, wave) {
    require(igraph)
    returnValue <- snaSociomatrixExtraction(i, data, sims, groupName, varName,
                                            wave)
    graph.adjacency(returnValue)
  }

  OutdegreeDistribution <- function (i, data, sims, groupName, varName, wave,
                                     levels=0:8, extractor=snaSociomatrixExtraction) {
    x <- extractor(i, data, sims, groupName, varName, wave)
    a <- apply(x, 1, sum)
    sapply(levels, function(i){ sum(a<=i) })
  }

  IndegreeDistribution <- function (i, data, sims, groupName, varName, wave,
                                    levels=0:8, extractor=snaSociomatrixExtraction) {
    x <- extractor(i, data, sims, groupName, varName, wave)
    a <- apply(x, 2, sum)
    sapply(levels, function(i){ sum(a<=i) })
  }

  GeodesicDistribution <- function (i, data, sims, groupName, varName, wave,
                                    extractor=snaEdgelistExtraction, levels=1:8) {
    require(sna)
    x <- extractor(i, data, sims, groupName, varName, wave)
    a <- geodist(x)$gdist
    sapply(levels, function(i){ sum(a<=i) })
  }

  TriadCensus <- function (i, data, sims, groupName, varName, wave,
                           extractor=snaEdgelistExtraction) {
    require(sna)
    x <- extractor(i, data, sims, groupName, varName, wave)
    triad.census(x)
  }

  KnnDistribution <- function (i, data, sims, groupName, varName, wave,
                               extractor=igraphEdgelistExtraction, levels=0:25) {
    require(igraph)
    x <- extractor(i, data, sims, groupName, varName, wave)
    a <- graph.knn(x)$knn
    a[is.nan(a)] <- 0
    sapply(levels, function(i){ sum(a<=i) })
  }
}
}
\keyword{models}
