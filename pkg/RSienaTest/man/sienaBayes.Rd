\name{sienaBayes}
\alias{sienaBayes}
\alias{sienaBayesFit}
\alias{bayes}
\alias{glueBayes}
\title{A function for fitting Bayesian models}
\description{A function to fit hierarchical Bayesian models random effects
  to sienaGroup data objects. Uses the function \code{\link{maxlikec}} for the
  MCMC part, the Bayesian part is performed in R.
}
\usage{
sienaBayes(data, effects, algo,  saveFreq=100,
    initgainGlobal=0.1, initgainGroupwise = 0.02, initfgain=0.2, initML=FALSE,
    priorMu=NULL, priorSigma=NULL, priorDf=NULL, priorKappa=NULL,
    priorRatesFromData=2,
    frequentist=FALSE, incidentalBasicRates=FALSE,
    reductionFactor=1, gamma=0.5, delta=1e-04,
    nwarm=50, nmain=250, nrunMHBatches=20,
    nSampVarying=1, nSampConst=1, nSampRates=0, nImproveMH=100,
    lengthPhase1=round(nmain/5), lengthPhase3=round(nmain/5),
    storeAll=FALSE, prevAns=NULL, prevBayes=NULL, silentstart = TRUE,
    nbrNodes=1, clusterType=c("PSOCK", "FORK"),
    getDocumentation=FALSE)

glueBayes(z1, z2, nwarm2=0)
}
\arguments{
  \item{data}{A sienaGroup object as returned by \code{\link{sienaGroupCreate}}.
    It is planned to also allow a \code{\link{siena}} data object as returned by
    \code{\link{sienaDataCreate}.}}
  \item{effects}{sienaEffects object as returned by \code{\link{getEffects}(data)}.
    The effects indicated by its column \code{randomEffects} get
    effects that are varying across the groups.}
  \item{algo}{Algorithm object, as created by
    \code{\link{sienaAlgorithmCreate}}. Should contain all options required
    for the MCMC scheme, and a random seed if required.}
  \item{saveFreq}{Integer. If this is larger than 1, the provisional results
    are saved after each multiple of saveFreq iterations in the main phase,
    in a file with name PartialBayesResult.RData (if a file with this name
    exists, it will be overwritten).
    This file contains an object z of class "sienaBayesFit",
    with the provisional results.
    This is to guard against crashes or power failure.
    It can be used as value for \code{prevBayes} as indicated below.}
  \item{initgainGlobal}{Step sizes in initial searches for good parameter values
    across the groups.}
  \item{initgainGroupwise}{Step sizes in initial searches for good parameter
    values by group; can be up to 0.1 for larger networks,
    0 for very small networks.}
  \item{initfgain}{Positive number, used only for frequentist estimation:
    the gain factor in the Robbins Monro algorithm is
    initfgain * ((iteration number)^(-gamma)).}
  \item{initML}{Boolean, whether to use maximum likelihood estimation for the
     groupwise initial estimation.}
  \item{priorMu}{Vector of length equal to the number of varying parameters.
     Prior mean of mu (global population mean for varying parameters); default: 0.
     For the basic rate parameters the defaults are data-dependent,
     but if ((\code{priorRatesFromData >= 1}) and \code{priorMu=NULL}),
     then for these parameters the prior mean will be 2.}
  \item{priorSigma}{Square matrix of dimension equal to length of priorMu.
     Prior global population covariance matrix Sigma for the varying parameters;
     default: identity matrix.}
  \item{priorDf}{Prior degrees of freedom for Sigma (global population
    covariance matrix); default: number of varying parameters + 2.}
  \item{priorKappa}{Proportionality constant between prior covariance matrix
    and covariance matrix of prior distribution for mu; default: 1.}
  \item{priorRatesFromData}{0, 1, or 2. Determines the prior distribution for the
     rate parameters. 0: prior is taken from \code{priorMu} and \code{priorSigma};
     1: prior is defined by mean and covariance matrix of estimated rate parameters
     from initialization phase;
     2: prior is defined by robust estimates for location and multivariate scale
    of estimated rate parameters from initialization phase.}
  \item{frequentist}{Currently only frequentist=FALSE works.
    Boolean: chooses between frequentist or Bayesian
    estimation of the global parameters.
    Frequentist estimation is possible only for at least 2 groups.}
  \item{incidentalBasicRates}{Boolean. If this is TRUE, the basic rate
    parameters are defined specifically for each group, and estimated
    using a Robbins Monro algorithm; if FALSE, they have a
    common prior for all groups.}
  \item{reductionFactor}{Positive number.
   If \code{priorRatesFromData} = 1 or 2, the prior distribution for the rate
   parameters is the estimated covariance matrix of the estimated rate parameters
   in the initialization phase multiplied by \code{reductionFactor}.
   For many data sets the default will be OK. If the rate parameters diverge,
   as evidenced by the traceplots of the posterior values,
   then a smaller value is advisable, such as 0.01.
   Has no effect if \code{priorRatesFromData} = 0.}
  \item{gamma}{Positive number, used only for frequentist estimation:
    the gain factor in the Robbins Monro algorithm is
    initfgain * ((iteration number)^(-gamma)).}
  \item{delta}{When the global population covariance matrix becomes
    non-positive definite (i.e., has one or more negative correlations)
    during iterations, it is changed so that
   all correlations are at least delta.}
  \item{nwarm}{Number of iterations in the warm up phase.
    Used only if \code{prevBayes} is \code{NULL}.
    Then it should be at least 5.}
  \item{nmain}{Number of iterations in the main phase.
     Should be at least 10.}
  \item{nrunMHBatches}{Thinning ratio in MCMC process.}
  \item{nSampVarying}{Number of samples of varying parameters
        for each chain sample.}
  \item{nSampConst}{Number of samples of constant parameters ("eta")
        for each chain sample.}
  \item{nSampRates}{Number of extra samples of basic rate parameters
        for each chain sample.}
  \item{nImproveMH}{Number of iterations per improveMH step.}
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
      \code{\link{siena07}} for the same data set and effects object.
      This then is the result of a multi-group estimation for these data,
      from which scaling information (derivative matrix and standard deviation of
      the deviations) will be extracted along with the parameters estimates
      which will be used as the initial values, unless
      \code{algo} requests the use of standard initial values.
      If \code{prevAns=NULL}, then a multi-group estimation will be performed
      as part of the initialization. }
  \item{prevBayes}{An object of class "sienaBayes" as returned by
      \code{sienaBayes}, on which the current function will continue.
      For example, the object z contained in the intermediate saved result
      PartialBayesResult.RData as mentioned above.
      Initialization and warming phases are skipped. If this is given,
      the values of \code{prevAns}, \code{nwarm}, and the prior
      distributions are disregarded.
      A new set of \code{nmain} iterations in the main phase are done.
      All further relevant parameters of the function call should be identical.
      This is not completely checked, so errors
      may occur if this is not the case.}
  \item{silentstart}{Boolean: whether to suppress most information
      to the console during the calculation of initial values.}
  \item{nbrNodes}{Number of processes to be used. Cannot be more than
      the number of periods summed over number of groups.}
   \item{clusterType}{If using multiple processes, whether to use
      forking processes or not. (Only "PSOCK" can be used on Windows.)}
   \item{getDocumentation}{Flag to allow documentation of internal
      functions, not for use by users.}
   \item{z1}{\code{sienaBayes} object.}
   \item{z2}{\code{sienaBayes} object with the same data and model specification
      as \code{z1}.}
   \item{nwarm2}{Number of warming iterations in \code{z2}; the first
      \code{nwarm2} iterations of \code{z2} will be left out of the
      combined object.}
}
\details{
  The function \code{sienaBayes} is for Bayesian estimation of one group or
  of multiple groups all having the same number of waves and the same
  model specification.
  It wraps Bayesian sampling of parameters around calls to
  \code{\link{maxlikec}}. The RSiena manual has a lot of information. \cr
  Effects can be randomly varying between groups, or global (i.e., constant
  across groups). This is indicated by the keyword \code{random} in the
  call of \code{\link{setEffect}}.\cr
  For the groupwise parameters normal distributions are assumed with conjugate
  priors. The prior distribution for the basic rate parameters is determined
  in a data-dependent way. For the non-varying parameters,
  a flat prior is assumed.\cr
  The frequentist option currently is not supported (may be broken).\cr
  The procedure consists of three parts: initialization, warming,
  main phase.\cr
  The required dimensions of the prior parameters \code{priorMu} and
  \code{priorSigma} are given in the output for
  \code{\link{print.sienaEffects}}.\cr
  In the initialization phase, initial parameter values and the proposal
  covariance matrix for Metropolis-Hastings steps for groupwise parameters
  are obtained from, first, Method of Moments estimation of a parameter vector
  assumed to be the same across the groups (in a multi-group estimation),
  followed by one subphase of the Robbins-Monro algorithm for Method of
  Moments estimates for the groups separately, with step size \code{initgain}.
  The proposal covariance matrices then are scaled, in the
  function 'improveMH', to achieve about 25 out of 100 acceptances of
  Bayes proposals after single MH steps.\cr
  After initialization and scaling of the proposal covariance matrices,
  a warming phase is done of \code{nwarm} Bayesian proposals
  each with a number of MH steps, followed again by the function 'improveMH'.\cr
  Finally \code{nmain} repeats (of \code{nrunMHBatches} of a
  number of MH steps sampling chains, plus \code{nSampVarying} MH steps
  sampling the varying parameters ('theta_j') plus \code{nSampConst} MH steps
  sampling the non-varying parameters ('eta') plus one Gibbs
  step sampling the global mean and covariance matrix of the varying parameters
  ('mu' and 'Sigma') are performed.
  In the warming as well as the final phase, the number of MH steps is
  determined by parameter \code{mult} ("multiplication factor")
  in the call of \code{sienaAlgorithmCreate}
  that created the algorithm object.\cr
  The function is time-consuming. When starting to use it, it is advisable
  to start with low values of \code{nmain} to explore computing time.\cr
  When the procedure seems to diverge, and for very small groups, it is
  advisable to use smaller values of the parameters \code{initgainGlobal}
  and \code{initgainGroupwise}.\cr
  \code{glueBayes(z1,z2)} combines two existing \code{sienaBayesFit} objects
  with the same data and model specifications (this is not
  checked) into one. This is intended to be used when \code{z2} was produced
  by calling \code{sienaBayes} with \code{prevAns=z1}.\cr
  For the function \code{\link{print.sienaEffects}}, the option
  \code{dropRates=TRUE} is meant to be useful especially for effects objects
  used for \code{sienaBayes}.
}
\value{
  \code{sienaBayes} as well as \code{glueBayes} return an object of class
  \code{sienaBayesFit}. This is a list containing, among other things:
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
  \item{ThinParameters}{array of dimensions (\code{nwarm+nmain} iterations
  by parameters by groups): sampled groupwise parameters}
  \item{ThinPosteriorMu}{array of dimensions (\code{nwarm+nmain} iterations
  by parameters): sampled global means of varying parameters}
  \item{ThinPosteriorEta}{array of dimensions (\code{nwarm+nmain} iterations
  by fixed parameters): sampled fixed parameters}
  \item{ThinPosteriorSigma}{array of dimensions (\code{nwarm+nmain} iterations
  by parameters by parameters): sampled global covariance matrix
  of varying parameters}
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
\seealso{\code{\link{siena07}}, \code{\link{sienaGroupCreate}}, \code{\link{bayesTest}},
   \code{\link{extract.sienaBayes}}\cr
  There are print and summary methods for sienaBayesFit objects,
  \code{\link{print.sienaBayesFit}}.}
\examples{
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
  FourEffects <- includeEffects(FourEffects, transTrip)
  FourEffects <- setEffect(FourEffects, density, random=TRUE)
  FourEffects <- setEffect(FourEffects, recip, random=TRUE)
  print(FourEffects, includeRandoms=TRUE)
  # Note this also shows the "randomEffects" column.
  FourAlgo <- sienaAlgorithmCreate(projname = "FourGroups", maxlike=TRUE,
                                   seed=123)
\dontrun{
  bayes.model <- sienaBayes(FourAlgo, data = FourGroups,
        effects = FourEffects, nwarm=10, nmain=25, nrunMHBatches=10)
  bayes.model
  summary(bayes.model)
  bayes.nextModel <- sienaBayes(FourAlgo, data = FourGroups,
        effects = FourEffects, nmain=15, nrunMHBatches=10,
        prevBayes = bayes.model)
  bayes.combinedModel <- glueBayes(bayes.model, bayes.nextModel)
  summary(bayes.combinedModel)
}
}
\keyword{models}
