\name{siena07}
\alias{siena07}
\title{Function to estimate parameters in a Siena model}
\description{
 Estimates parameters in a Siena model using method of moments, based on direct
 simulation, conditional or otherwise; or using Maximum Likelihood by MCMC
 simulation. Estimation is done using a Robbins-Monro algorithm. Note that
 the data and particular model to be used
 must be passed in using named arguments as the \code{...},
 and the specification for the algorithm must be passed on as x, which is
 a \code{\link{sienaAlgorithm}} object as produced by
 \code{\link{sienaAlgorithmCreate}} (see examples).}
\usage{
siena07(x, batch=FALSE, verbose=FALSE, silent=FALSE,
        useCluster=FALSE, nbrNodes=2, initC=TRUE,
        clusterString=rep("localhost", nbrNodes), tt=NULL,
        parallelTesting=FALSE, clusterIter=!x$maxlike,
        clusterType=c("PSOCK", "FORK"),
        logLevelConsole='WARNING', logLevelFile='INFO',
        logIncludeLocation=F, ...)
      }
\arguments{
  \item{x}{A control object, of class \code{\link{sienaAlgorithm} }}
  \item{batch}{ Desired interface: FALSE gives a gui (graphical user interface
     implemented as a tcl/tk screen),
     TRUE gives a small (if \code{verbose=FALSE}) amount of printout
     to the console.}
  \item{verbose}{Produces various output to the console if TRUE.}
  \item{silent}{Produces no output to the console if TRUE, even if batch mode.}
  \item{useCluster}{Boolean: whether to use a cluster of processes (useful if
  multiple processors are available).}
  \item{nbrNodes}{Number of processes to use if useCluster is TRUE.}
  \item{initC}{Boolean: set to TRUE if the simulation will use C
   routines (currently always needed). Only for use if using
   multiple processors, to ensure all copies are initialised
   correctly. Ignored otherwise, so is set to TRUE by default.}
  \item{clusterString}{Definitions of clusters. Default set up to use
   the local machine only.}
  \item{tt}{A \code{tcltk} toplevel window. Used if called from the
   model options screen.}
  \item{parallelTesting}{Boolean. If TRUE, sets up random numbers to
   parallel those in Siena 3.}
  \item{clusterIter}{Boolean. If TRUE, multiple processes execute
   complete iterations at each call.
   If FALSE, multiple processes execute a single wave at each call.}
  \item{clusterType}{Either "PSOCK" or "FORK". On Windows, must be
  "PSOCK". On a single non-Windows machine may be "FORK", and
  subprocesses will be formed by forking. If "PSOCK", subprocesses are
  formed using \R scripts.}
  \item{logLevelConsole}{
    Minimal level of logging output to the console. One of 'DEBUG', 'VERBOSE',
    'INFO', 'WARNING', 'ERROR', 'FATAL'. Where 'DEBUG' prints all messages,
    'FATAL' only the messages with highest priority.}
  \item{logLevelFile}{
    Like \code{logLevelConsole} but for logging to files.}
  \item{logIncludeLocation}{
    When \code{TRUE} also log the code location (file, method, line) where the
    logging message originated.}
  \item{\dots}{Arguments for the simulation function, see
    \code{\link{simstats0c}}:
    usually, \code{data} and \code{effects}, as in the examples below;\cr
    possibly also \code{prevAns} if a previous reasonable
    provisional estimate was obtained for a similar model;\cr
    possibly also \code{returnDeps} if the simulated dependent variables
    (networks, behaviour) should be returned;\cr
    possibly also \code{returnChains} if the simulated sequences (chains)
    of ministeps should be returned; this will produce a very big file.}
}
\details{
 \code{siena07} runs a Robbins-Monro algorithm for parameter estimation
 according to the Method of Moments using the three-phase
 implementation in Snijders (2001) and
 Snijders, Steglich and Schweinberger (2007),
 with (if \code{findiff=FALSE} in the \code{\link{sienaAlgorithm}} object)
 derivative estimation as in Schweinberger and Snijders (2007).
 Phase 1 does a few iterations to
 estimate the derivative matrix of the targets with respect to the
 parameter vector. Phase 2 does the estimation. Phase 3 runs a
 simulation to estimate standard errors and check convergence of the model. The
 simulation function is called once for each iteration in these phases
 and also once to initialise the model fitting and once to complete it.
 Unless in batch mode, displays a tcl/tk screen to allow interruption
 and to show progress.

 If  \code{maxlike=TRUE} in the \code{\link{sienaAlgorithm}} object,
 estimation is done by Maximum Likelihood implemented as in
 Snijders, Koskinen and Schweinberger (2010) also using the three-phase
 Robbins-Monro algorithm.

 It is necessary to check that convergence has been achieved.
 The rule of thumb is that the all t-ratios for convergence should be
 in absolute value less than 0.1 and
 the overall maximum convergence ratio should be less than 0.25.
 If this was not achieved, the result can be used to start
 another estimation run from the estimate obtained, using
 the parameter \code{prevAns} as illustrated in the example below.\cr
 (This parameter is inserted through '\dots' into the function
 \code{\link{initializeFRAN}}.)

 An alternative implementation is in \code{\link{sienacpp}},
 where the algorithm is coded entirely in C++;
 in \code{siena07} the 'front end' is coded in R.
}
\value{
  Returns an object of class \code{\link{sienaFit}}, some parts of which are:
  \item{OK}{Boolean indicating successful termination}
  \item{termination }{Character string, values: "OK", "Error", or
    "UserInterrupt". "UserInterrupt" indicates that the user asked for early
  termination before phase 3. }
  \item{f}{Various characteristics of the data and model definition.}
  \item{theta}{Estimated value of theta.}
  \item{covtheta}{Estimated covariance matrix of theta; this is not available
    if the sienaAlgorithm object \code{x} was produced with \code{simOnly=TRUE}.}
  \item{se}{Vector of standard errors of estimated theta.}
  \item{dfra}{Matrix of estimated derivatives.}
  \item{sf}{Matrix of deviations from targets in phase 3.}
  \item{sf2}{Array of periodwise deviations from simulations in phase 3.}
  \item{tconv}{t-statistics for convergence.}
  \item{tmax}{maximum absolute t-statistic for convergence for non-fixed parameters.}
  \item{tconv.max}{overall maximum convergence ratio.}
  \item{targets}{Observed statistics; for ML, zero vector.}
  \item{targets2}{Observed statistics by wave,
    starting with second wave; for ML, zero matrix.}
  \item{ssc}{Score function contributions for each wave for each
    simulation in phase 3. Zero if finite difference method is used.}
 \item{sims}{If \code{returnDeps} is \code{TRUE}:
    list of simulated dependent variables (networks, behaviour).
    Networks are given as a list of edgelists, one for each period.}
 \item{chain}{If \code{returnChains = TRUE}:
    list, or data frame, of simulated chains of ministeps.
    The chain has the structure
    \code{chain[[run]][[depvar]][[period]][[ministep]]}.
    }
   \item{Phase3nits}{Number of iterations actually performed in phase 3.}
  Writes text output to the file named "projname.out", where projname is defined
  in the sienaAlgorithm object \code{x}.
}
\references{\itemize{
   \item Schweinberger, Michael, and Snijders,
    Tom A.B. (2007). Markov models for digraph panel data:
    Monte Carlo-based derivative estimation.
    \emph{Computational Statistics and Data Analysis} 51, 4465-4483.
    \item Snijders, Tom A.B. (2001).
    The statistical evaluation of social network dynamics.
    \emph{Sociological Methodology}, 31, 361-395.
   \item Snijders, Tom A. B., Koskinen, Johan, and Schweinberger, Michael (2010).
    Maximum likelihood estimation for social network dynamics.
    \emph{Annals of Applied Statistics}, 4, 567-588.
   \item Snijders, Tom A.B., Steglich, Christian E.G., and Schweinberger,
    Michael (2007). Modeling the co-evolution of networks and behavior.
    Pp. 41-71 in \emph{Longitudinal models in the behavioral
    and related sciences},
    edited by Kees van Montfort, Han Oud and Albert Satorra; Lawrence Erlbaum.
   \item Steglich, Christian E. G., Snijders, Tom A. B., and
    Pearson, Michael A. (2010). Dynamic networks and behavior:
    Separating selection from influence. Sociological Methodology, 40, 329-393.
   \item Information about the implementation of the algorithm is in
     \url{http://www.stats.ox.ac.uk/~snijders/siena/Siena_algorithms.pdf}.
   \item Further see \url{http://www.stats.ox.ac.uk/~snijders/siena/} .
}
}
\author{Ruth Ripley, Tom Snijders}

\seealso{\code{\link{sienaAlgorithmCreate}}, \code{\link{sienacpp}}.\cr
  There are print, summary and xtable methods for \code{\link{sienaFit}}
  objects: \code{\link{xtable}}, \code{\link{print.sienaFit}}.}
\examples{
myalgorithm <- sienaAlgorithmCreate(nsub=2, n3=100)
# nsub=2 and n3=100 is used here for having a brief computation, not for practice.
mynet1 <- sienaDependent(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
ans <- siena07(myalgorithm, data=mydata, effects=myeff, batch=TRUE)

# or for non-conditional estimation
\dontrun{
model <- sienaModelCreate(nsub=2, n3=100, cond=FALSE)
ans <- siena07(myalgorithm, data=mydata, effects=myeff, batch=TRUE)
        }

# or if a previous "on track" result ans was obtained
\dontrun{
ans1 <- siena07(myalgorithm, data=mydata, effects=myeff, prevAns=ans)
         }
}
\keyword{models}
