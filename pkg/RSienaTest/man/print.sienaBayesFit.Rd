\name{print.sienaBayesFit}
\alias{print.sienaBayesFit}
\alias{summary.sienaBayesFit}
\alias{print.summary.sienaBayesFit}
\title{Methods for processing sienaBayes objects}
\usage{
\method{print}{sienaBayesFit}(x, nfirst=NULL, ...)
\method{summary}{sienaBayesFit}(object, nfirst=NULL, ...)
\method{print}{summary.sienaBayesFit}(x, nfirst=NULL, ...)
}
\description{
  \code{print} and \code{summary} methods for \code{\link{sienaBayesFit}} objects.
}
\arguments{
  \item{object}{An object of class \code{\link{sienaBayesFit}}.}
  \item{x}{An object of class  \code{\link{sienaBayesFit}} or
    \code{summary.sienaBayesFit} as appropriate.}
  \item{nfirst}{The first element of the MCMC chain used for calculating
      properties of the chain; i.e., the first element for which
      it is assumed that convergence has occurred. If \code{NULL},
      \code{x$nwarm} will be used.}
  \item{\dots}{For extra arguments (none used at present).}
}
\value{
  The function \code{print.sienaBayesFit} prints results (estimated posterior
  and posterior standard deviations) of a Bayesian analysis, as resulting
  from function \code{\link{sienaBayes}}.\cr
  The function \code{summary.sienaBayesFit} prints more extensive results,
  including groupwise results.\cr
  The columns "cred. from" and "cred. to" give the left and right end points
  of the estimated 95 percent credibility interval for the parameters, and the
  column "p" gives the estimated posterior probability that the parameter is
  greater than 0.

}
\seealso{\code{\link{sienaBayes}}}
\references{
  See also \url{http://www.stats.ox.ac.uk/~snijders/siena/}
}
\author{Tom Snijders}
\keyword{ methods }
