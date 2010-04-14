\name{sienaTimeTest}
\alias{sienaTimeTest}
\alias{plot.sienaTimeTest}
\title{Functions to assess and account for time heterogeneity of parameters}
\description{
 Takes a \code{sienaFit} object from a \code{\link{siena07}} estimation
 and tests the addition of dummy parameters at waves \code{m=2...(M-1)}
 through the score type test of Schweinberger (2007b).
 Tests for joint significance,
 parameter-wise significance, individual significance, and one step
 estimation of the unrestricted model parameters are returned in a list.

 Plotting facilities are also included which permit visual
 representation of the results.

 After assessing time heterogeneity, effects objects can be modified via
 the \code{timeDummy} column. Simply type the periods which you would
 like to have time dummied separated by commas into this column. Using
 the \code{timeDummy} column within the effects object, you may specify
 which time dummy interacted parameters are to be estimated.

 Note that this function is only appropriate for sienaFit objects which
 use the score function method of derivative estimation of Schwienberger
 and Snijders (2007a).
 }

\usage{
timetest <- sienaTimeTest(sienaFit)
\method{plot}{sienaTimeTest}(x, pairwise=FALSE, effects=1:2,
	dims=c(2, 1), scale=.2, plevels=c(.1, .05, .025), multiplot=FALSE, ...)
sienaEffects$timeDummy <- "2,4,6" #Include time dummies for periods
                                  # 2, 4, and 6
sienaEffects$timeDummy <- "all"	  #Include all time dummies
}
\arguments{
  \item{sienaFit}{ A sienaFit object returned by \code{siena07}. }
  \item{x}{ A sienaTimeTest object returned by \code{sienaTimeTest}. }
  \item{pairwise}{ A Boolean value corresponding to whether the user
  would like a pairwise plot of the simulated statistics to assess
  correlation among the effects (\code{pairwise=TRUE}), or a plot of the
  estimates across waves in order to assess graphically the results of
  the score type test. }
\item{effects}{ A vector of integers corresponding to the indices given
  in the \code{sienaTimeTest} output for effects which are to be
  plotted.}
\item{dims}{ If multiple effects are to be plotted, a vector of the
  form (rows, columns) corresponding to how the plot is to lay out the
  effects. }
\item{scale}{ A positive number corresponding to the number of standard
  deviations on one step estimates to use for computing the maximum and
  minimum of the plotting range. We recommend experimenting with this
  number when the y-axes of the plots are not satisfactory. Smaller
  numbers shrink the axes.}
\item{plevels}{ A list of three decimals indicating the gradients at
  which to draw the confidence interval bars.}
\item{multiplot}{Boolean: are the plots to be one per page (TRUE) or not}
 \item{\dots}{For extra arguments (none used at present)}
}
\details{

  This test follows the score type test of Schweinburger (2007b) as
 implemented by Lospinoso et. al. (2010) by using statistics already
 calculated at each wave to populate vectors of partitioned moment
 functions corresponding to a restricted model (the model that has been
 fit by \code{sienaFitObject}) and an unrestricted model (which contains
 dummies for \code{m=2...(M-1)}. A covariance matrix of these statistics
 is calculated, and a delta matrix is constructed through the score
 functions' outer products with these statistics. Through an
 orthogonalization and a Delta method of approximation, the
 variance-covariance structure of the dummy statistics is calculated,
 and appropriate statistical tests can be used.

 The \code{pairwise=TRUE} plot may be used to assess whether effects are
 highly correlated. This information may be important when considering
 forward-model selection, since highly correlated effects may have
 highly correlated one-step estimates, particularly since the individual
 score type tests are not orthogonalized against the scores and
 deviations of yet-unestimated dummies.  For example, reciprocity and
 outdegree may have highly correlated statistics as indicated by a
 strong, positive correlation coefficient.  When considering whether to
 include dummy terms, it may be a good idea to include, e.g., outdegree,
 estimate the parameter, and see whether reciprocity dummies remain
 significant after method of moments estimation of the updated model--as
 opposed to including both outdegree and reciprocity.

 The \code{pairwise=FALSE} plot displays the most of the information
 garnered from \code{sienaTimeTest} in a graphical fashion.  For a each
 effect, the method of moments parameter estimate for the base period
 (i.e. wave 1) is given as a blue, horizontal reference line. One step
 estimates are given for all of the parameters by dots at each wave. The
 dots are colored black if the parameter has been included in the model
 already (i.e. has been estimated via method of moments), or red if they
 have not been included. Confidence intervals are given based on pivots
 given at \code{pvalues}. Evidence of time heterogeneity is suggested by
 points with confidence intervals not overlapping with the base period.

 If it is determined that a time heterogeneity occurs for any number of
 time periods or effects (or any combination therein), the
 \code{timeDummy} column provides facilities within a sienaFit object to
 quickly adjust the model and re-estimate in the usual way.

}
\value{
  \code{sienaTimeTest} Returns a list containing many items,
  including the following:
  \item{JointTest }{A \code{chi^2} test for joint significance of the dummies.}
  \item{ParameterTest }{ A \code{chi^2} test for joint significance
    across dummies for each separate effect.}
  \item{IndividualTest }{ A matrix displaying initial estimates, one
	step estimates, and a p-value for \code{H0:} the unrestricted
	parameters are equal to zero.}
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}
  for general information on RSiena.

Lospinoso, J.A., Schweinberger, M., Snijders, T.A.B, and Ripley,
R.M. "Assessing and Accounting for Time Heterogeneity in Stochastic
Actor Oriented Models". Advances in Data Analysis and Computation.
Special Issue on Social Networks. Submitted. Available from
\url{http://www.stats.ox.ac.uk/~lospinos/}.

M. Schweinberger and T.A.B. Snijders. Markov models for digraph panel
data: Monte carlo-based derivative estimation. Comput. Stat. Data Anal.,
51(9):4465-4483, 2007.  ISSN 0167-9473.

M. Schweinberger. Chapter 4: Statistical Modeling of Network Panel Data:
Goodness of Fit. PhD thesis, University of Groningen, 2007.}
\author{Josh Lospinoso}
\seealso{\code{\link{siena07}}, \code{\link{sienaTimeTest}}}
\examples{
## Estimate a restricted model
mymodel <- sienaModelCreate(fn=simstats0c, nsub=2, n3=100)
mynet1 <- sienaNet(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
myeff$include[myeff$shortName=='transTrip' & myeff$type=='eval'] <- TRUE
myeff$include[myeff$shortName=='balance' & myeff$type=='eval'] <- TRUE
ans <- siena07(mymodel, data=mydata, effects=myeff, batch=TRUE)

## Conduct the score type test to assess whether heterogeneity is present.
tt <- sienaTimeTest(ans)

## Suppose that we wish to include two time dummies.
## Add them in the following way:
myeff$timeDummy[myeff$shortName=='recip' & myeff$type=='eval'] = "2"
myeff$timeDummy[myeff$shortName=='balance' & myeff$type=='eval'] = "2"
ans2 <- siena07(mymodel, data=mydata, effects=myeff, batch=TRUE)

## Re-assess the time heterogeneity
tt2 <- sienaTimeTest(ans2)

## And so on..

## A demonstration of the plotting facilities, on a larger dataset:
\dontrun{
mymodel <- sienaModelCreate(fn=simstats0c, nsub=4, n3=500)
mynet1 <- sienaNet(array(c(s501, s502, s503, s501, s503, s502), dim=c(50, 50, 6)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
myeff$include[myeff$shortName=='transTrip' & myeff$type=='eval'] <- TRUE
myeff$include[myeff$shortName=='balance' & myeff$type=='eval'] <- TRUE
myeff$timeDummy[myeff$shortName=='recip' & myeff$type=='eval'] = "2,3,5"
myeff$timeDummy[myeff$shortName=='balance' & myeff$type=='eval'] = "2,5"
myeff$timeDummy[myeff$shortName=='recip' & myeff$type=='eval'] = "4"
myeff$timeDummy[myeff$shortName=='density' & myeff$type=='eval'] = "all"
ansp <- siena07(mymodel, data=mydata, effects=myeff, batch=FALSE)
ttp <- sienaTimeTest(ansp)

## Pairwise plots show
plot(ttp, pairwise=TRUE)

## Time test plots show
plot(ttp, effects=1:4, dims=c(2,2))
}
}
\keyword{models}
