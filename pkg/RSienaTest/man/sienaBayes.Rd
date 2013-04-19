\name{sienaBayes}
\alias{sienaBayes}
\alias{bayes}
\title{A function for fitting Bayesian models}
\description{A function to fit Bayesian models to Siena Data
  objects, including fitting hierarchical random effects models
  to sienaGroup data objects. Uses the function \code{\link{maxlikec}} for the MCMC part,
  the Bayesian part is performed in R. Still under development.
}
\usage{
sienaBayes(data, effects, model,  saveFreq=100,
    initgainGlobal=0.1, initgainGroupwise = 0.02, initfgain=0.01,
    priorMu=NULL, priorSigma=NULL, priorDf=NULL, priorKappa=NULL,
    frequentist=FALSE, incidentalBasicRates=FALSE,
	gamma=0.5, delta=1e-04,
	nwarm=50, nmain=250, nrunMHBatches=20,
	lengthPhase1=round(nmain/5), lengthPhase3=round(nmain/5),
    storeAll=FALSE, prevAns=NULL,
	plotit=FALSE, nbrNodes=1, clusterType=c("PSOCK", "FORK"),
	getDocumentation=FALSE)
}
\arguments{
  \item{data}{A sienaData object as returned by \code{\link{sienaDataCreate}},
    or a sienaGroup object as returned by \code{\link{sienaGroupCreate}}.}
  \item{effects}{List of sienaEffects objects as
    returned by \code{\link{getEffects}}.}
  \item{model}{Algorithm object, as created by
    \code{\link{sienaAlgorithmCreate}}. Should contain all options required
    for the MCMC scheme, and a random seed if required.}
  \item{saveFreq}{Integer: if this is larger than 1, the provisional results
    are saved after each multiple of saveFreq iterations in the main phase,
    in a file with name PartialBayesResult.RData (if a file with this name
    exists, it will be overwritten).
    This is to guard against crashes or power failure.}
  \item{initgainGlobal}{Step sizes in initial searches for good parameter values
  across the groups.}
  \item{initgainGroupwise}{Step sizes in initial searches for good parameter
  values by group; can be up to 0.1 for larger networks,
  0 for very small networks.}
  \item{initfgain}{Positive number, used only for frequentist estimation:
  the gain factor in the Robbins Monro algorithm is
  initfgain * ((iteration number)^(-gamma)).}
  \item{priorMu}{Prior mean of mu (global population mean); default:
  0, but data-dependent for basic rate parameters.}
  \item{priorSigma}{Prior global population covariance matrix Sigma; default:
  identity matrix.}
  \item{priorDf}{Prior degrees of freedom for Sigma (global population
    covariance matrix); default: number of parameters + 2.}
  \item{priorKappa}{Proportionality constant between prior covariance matrix
    and covariance matrix of prior distribution for mu; default: 1.}
  \item{frequentist}{Boolean: chooses between frequentist or Bayesian
  estimation of the global parameters.
  Frequentist estimation is possible only for at least 2 groups.}
  \item{incidentalBasicRates}{Boolean: if this is TRUE, the basic rate
  parameters are defined specifically for each group, and estimated
  using a Robbins Monro algorithm; if FALSE, they are
  common for all groups.}
  \item{gamma}{Positive number, used only for frequentist estimation:
  the gain factor in the Robbins Monro algorithm is
  initfgain * ((iteration number)^(-gamma)).}
  \item{delta}{When the global population covariance matrix becomes
  non-positive definite (i.e., has one or more negative correlations)
  during iterations, it is changed so that
  all correlations are at least delta.}
  \item{nwarm}{Number of iterations in the warm up phase.}
  \item{nmain}{Total number of iterations.
  Should be strictly larger than nwarm + 5.}
  \item{nrunMHBatches}{Thinning ratio in MCMC process.}
  \item{lengthPhase1}{Only used for frequentist estimation:
  length of the first phase of the Robbins Monro algorithm.
  lengthPhase1 + lengthPhase3 should be strictly less than nmain.}
  \item{lengthPhase3}{Only used for frequentist estimation:
  length of the third phase of the Robbins Monro algorithm.
  lengthPhase1 + lengthPhase3 should be strictly less than nmain.}
  \item{storeAll}{Boolean: whether to store parameters for all MCMC iterations,
  i.e., before thinning. \code{storeAll=TRUE} may lead to producing
  very large objects and is not recommended for usual operation.}
  \item{prevAns}{An object of class "sienaFit" as returned by
    \code{\link{siena07}}, from which scaling
    information (derivative matrix and standard deviation of the
    deviations) will be extracted along with the latest version of the
    parameters which will be used as the initial values, unless the
    model requests the use of standard initial values.
    Influences only the start of the initialization.
    If the results used as \code{prevAns} are a reasonable starting point,
    this will increase the efficiency of the algorithm. }
  \item{plotit}{Boolean: whether to plot parameters during the run}
  \item{nbrNodes}{Number of processes to be used. Cannot be more than
    the number of waves summed over number of groups.}
   \item{clusterType}{If using multiple processes, whether to use
     forking processes or not. (Only "PSOCK" can be used on Windows.)}
   \item{getDocumentation}{Flag to allow documentation of internal
    functions, not for use by users.}
}
\details{
  This function is for Bayesian estimation of one group or
  of multiple groups all having the same number of waves.
  It wraps Bayesian sampling of parameters around calls to
  \code{\link{maxlikec}}.\cr
  Initial parameter values and the proposal covariance matrix
  for Metropolis-Hastings steps for groupwise parameters are obtained
  from, first, Method of Moments estimation of a parameter assumed to be
  the same across the groups, followed by one subphase of the
  Robbins-Monro algorithm for Method of Moments estimates for the groups
  separately, with step size \code{initgain}.
  The proposal covariance matrices then are scaled to achieve about
  40 out of 100 acceptances of Bayes proposals after single MH steps.\cr
  For the groupwise parameters normal distributions are assumed.
  The prior distribution for the basic rate parameters is determined
  in a data dependent way.\cr
  After initialization and scaling of the proposal covariance matrices,
  a warming phase is done of \code{nwarm} Bayesian proposals
  each with a number of MH steps.
  Finally \code{nmain - nwarm} repeats (of \code{nrunMHBatches} of a
  number of MH steps plus Bayesian proposals) are performed.
  In the warming as well as the final phase, the number of MH steps is
  determined by parameter \code{mult} ('multiplication factor')
  in the call of \code{sienaAlgorithmCreate} that created the algorithm object.
  If \code{plotit} is TRUE,
  plots are produced at intervals during the run showing progress. They
  can be memory hogs, but the code may be useful for use on the returned
  values.
}
\value{
  Returns an object of class "sienaBayesFit".
  This is a list containing, among other things:
  \item{priorMu}{prior global population mean (not quite the same
  as corresponding input parameter)}
  \item{priorSigma}{prior global population covariance matrix
  (not quite the same as corresponding input parameter)}
  \item{priorKappa}{proportionality constant between prior covariance matrix
                and covariance matrix of prior distribution for the mean}
  \item{priorDf}{prior degrees of freedom for covariance matrix}
  \item{effectName}{array of names of effects included in the model}
  \item{f$groupNames}{array of names of groups included in the model}
  \item{initialResults}{sienaFit object: result of abbreviated MoM estimation
  under the assumption of same parameters across groups}
  \item{ThinParameters}{array of dimensions (nmain iterations by parameters
  by groups): sampled groupwise parameters}
  \item{ThinPosteriorMu}{array of dimensions (nmain iterations by parameters):
  sampled global mean parameters}
  \item{ThinPosteriorSigma}{array of dimensions (nmain iterations by parameters
  by parameters): sampled global covariance matrix parameters}
  \item{candidates}{if \code{storeAll=TRUE}: array of proposed parameters,
  by group}
  \item{acceptances}{if \code{storeAll=TRUE}: matrix of booleans:
   whether the corresponding change to the parameters was accepted, by group}
  \item{MHacceptances}{if \code{storeAll=TRUE}:
  array of acceptances of the MH steps, by step    type and group
  but summed over dependent variables}
  \item{MHrejections}{if \code{storeAll=TRUE}:
  array of rejections of the ML steps}
  \item{MHproportions}{if \code{storeAll=TRUE}:
  array of proportions of the MH steps accepted}
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}

  Koskinen, J. H. and T. A. B. Snijders (2007).
  Bayesian inference for dynamic social network data.
  \emph{Journal of Statistical Planning and Inference}, 13, 3930-3938.
}
\author{Ruth Ripley, Johan Koskinen, Tom Snijders }
\seealso{\code{\link{siena07}}, \code{\link{sienaGroupCreate}} \cr
  There is a print method for sienaBayesFit objects,
  \code{\link{print.sienaBayesFit}}.}
\examples{
\dontrun{
  Group1 <- sienaDependent(array(c(N3401, HN3401), dim=c(45, 45, 2)))
  Group3 <- sienaDependent(array(c(N3403, HN3403), dim=c(37, 37, 2)))
  Group4 <- sienaDependent(array(c(N3404, HN3404), dim=c(33, 33, 2)))
  Group6 <- sienaDependent(array(c(N3406, HN3406), dim=c(36, 36, 2)))
  dataset.1 <- sienaDataCreate(Friends = Group1)
  dataset.3 <- sienaDataCreate(Friends = Group3)
  dataset.4 <- sienaDataCreate(Friends = Group4)
  dataset.6 <- sienaDataCreate(Friends = Group6)
  FourGroups <- sienaGroupCreate(
        list(dataset.1, dataset.3, dataset.4, dataset.6))
  FourEffects <- getEffects(FourGroups)
  FourModel <- sienaModelCreate(projname = 'FourGroups', maxlike=TRUE)
  bayes.model <- sienaBayes(FourModel, data = FourGroups,
        effects = FourEffects, nwarm=10, nmain=25, nrunMHBatches=10)
}
}
\keyword{models}
