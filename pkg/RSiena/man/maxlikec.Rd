\name{maxlikec}
\alias{maxlikec}
\title{A version of FRAN}
\description{
 The function to be called as 'FRAN' by \code{\link{siena07}}
 when using estimation by Maximum likelihood. Calls compiled C++.
}
\usage{
maxlikec(z, x, INIT = FALSE, TERM = FALSE, initC=FALSE,
           data=NULL, effects=NULL, profileData=FALSE,
           prevAns=NULL, returnChains=FALSE, byGroup=FALSE,
           returnDataFrame=FALSE, byWave=FALSE)
}
\arguments{
  \item{z}{Control object, passed in automatically by
	\code{\link{siena07}} or \code{\link{bayes}}.}
  \item{x}{Model object, passed in automatically by
	\code{\link{siena07}} or \code{\link{bayes}}.}
  \item{INIT}{If TRUE, do initial processing. May be required to set up z.}
  \item{TERM}{If TRUE, do end processing. }
  \item{initC}{If TRUE, call is to setup the data and model in C++. For
	use with multiple processes only.}
  \item{data}{A sienaData object as returned by \code{\link{sienaDataCreate}}.}
  \item{effects}{A sienaEffects object as
	returned by \code{\link{getEffects}}.}
  \item{profileData}{Boolean to force dumping of the data for profiling
	with \command{sienaProfile.exe}.}
  \item{prevAns}{An object of class "sienaFit" as returned by
   \code{\link{siena07}}, from which scaling
	information (derivative matrix and standard deviation of the
	deviations) will be extracted along with the latest version of the
	parameters which will be used as the initial values, unless the
	model requests the use of standard initial values. If the previous
	model is exactly the same as the current one, Phase 1 will be omitted. If
	not, any parameter estimates for effects which are included in
	the new model will be used as initial
	values, but phase 1 will still be carried out.
	If the results used as \code{prevAns} are a reasonable starting point,
	this will increase the efficiency of the algorithm. }
  \item{returnChains}{Boolean. Whether to return the chains.}
  \item{byGroup}{Boolean. For internal use: allows different thetas for each
	group to be used in \code{\link{bayes}}.}
  \item{returnDataFrame}{Boolean. Whether to return the chains as a list
	or a data frame. For debugging use only: will cause multiple
	processes to fail.}
  \item{byWave}{Boolean: whether to return the Hessian for each wave
	separately. Necessary for using \code{\link{sienaTimeTest}}.}
}
\details{
  This should be used for the element FRAN of the model object
  estimation by Method of Moments.
  The arguments with no defaults must be passed in on the call to
  \code{\link{siena07}} or \code{\link{bayes}}.
}
\value{Depends on the call. If \code{INIT} or \code{initC} or
  \code{TERM}
  are true, returns \code{z}, the control object. Otherwise,
  returns a list containing:
  \item{fra}{Simulated scores.}
  \item{dff}{Simulated Hessians: stored as lower triangular matrices}
  \item{ntim}{NULL, compatibility only}
  \item{feasible}{Currently set to TRUE.}
  \item{OK}{Could be set to FALSE if serious error has occurred.}
  \item{dff}{Simulated Hessian}
  \item{sims}{NULL, for compatibility only}
  \item{chain}{A list of simulated chains, one for each iteration. Each
	list consists of a list for each data object, each of which consists of
	a list for each network, each of which consists of a list for each period,
	each component of which is a list or a data frame depending on the
	value of \code{returnDataFrame}.
	Only if \code{returnDeps} is TRUE.}
  \item{accepts}{Number of accepted MH steps by dependent variable
	(permute steps are counted under first dependent variable)}
  \item{rejects}{Number of rejected MH steps by dependent variable
 	(permute steps are counted under first dependent variable)}
 \item{aborts}{Number of aborted MH steps by dependent variable
	(permute steps are counted under first dependent variable)}
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley }
\seealso{\code{\link{siena07}} , \code{\link{simstats0c}}, \code{\link{bayes}}}
\examples{
\donttest{
mynet1 <- sienaNet(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTrip)
mymodel <- sienaModelCreate(maxlike=TRUE, nsub=2, n3=100)
ans< - siena07(mymodel, data=mydata, effects=myeff,  batch=TRUE)
}
}
\keyword{models}
