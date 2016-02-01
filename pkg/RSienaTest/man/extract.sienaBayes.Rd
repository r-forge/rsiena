\name{extract.sienaBayes}
\alias{extract.sienaBayes}
\alias{extract.Bayes}
\title{
  Extract posterior samples from sienaBayes results
}
\description{
  This function extracts posterior samples from a \code{\link{sienaBayesFit}} object
  to be used, e.g., for assessing convergence.
}
\usage{
extract.sienaBayes(zlist, nfirst=zlist[[1]]$nwarm+1, extracted,
                   sdLog=TRUE)
}
\arguments{
  \item{zlist}{A list of \code{\link{sienaBayesFit}} objects,
    further called 'chains',
    resulting from calls to \code{\link{sienaBayes}} with a common
    data set and model specification.}
  \item{nfirst}{Integer: the first element for the first MCMC chain
      used for calculating properties of the chain.}
  \item{extracted}{The parameters for which posterior samples are to be
     extracted: \cr
     \code{"all"}: all parameters;\cr
     \code{"rates"}: all groupwise rate parameters;\cr
     \code{"varying"}: all varying non-rate parameters:
                      global means and standard deviations;\cr
     \code{"non-varying"}: all estimated non-varying
                   (and therefore, non-rate) parameters;\cr
     \code{"objective"}: all non-rate parameters.
     }
  \item{sdLog}{Logical: should the logarithms be taken of the posterior
      standard deviations (applies only to varying parameters).}
}
\details{
 This function produces a 3-dimensional array, iterations by chains
 by parameters, suitable for use, e.g., in function \code{monitor()}
 of package \code{stan}. The 'iterations' are draws from the
 posterior distribution (provided that convergence was achieved),
 after the thinning in \code{sienaBayes} implied by its parameters
 \code{nrunMHBatches}, \code{nSampVarying}, \code{nSampConst} and
 \code{nSampRate}.\cr
 The number of iterations extracted is determined by the first chain and
 \code{nfirst}. The same number of iterations will be extracted
 from the other chains; for each, the last part of the chain
 will be used. If one of the chains is not long enough,
 this will lead to an error message.\cr
 It is not checked that the specifications of the chains are identical.
 For the parameter names, the names in the first chain are used.
}
\value{A 3-dimensional array, iterations by chains by parameters,
  the third dimension having the names of the parameters.
  }
\references{See the manual and \url{http://www.stats.ox.ac.uk/~snijders/siena/}
}
\author{Tom Snijders}
\seealso{\code{\link{sienaBayes}}}
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
  FourAlgo <- sienaAlgorithmCreate(projname = "FourGroups", maxlike=TRUE)
\dontrun{
  bayes.models <- replicate(4, sienaBayes(FourAlgo, data = FourGroups,
        effects = FourEffects, nwarm=10, nmain=25, nrunMHBatches=10),
        simplify=FALSE)
  bayes.extracted <- extract.sienaBayes(bayes.models, extracted="all")
}
}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{ tests }