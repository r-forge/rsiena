\name{sienaRI}
\alias{sienaRI}
\alias{plot.sienaRI}
\alias{print.sienaRI}
\alias{sienaRI.methods}
\title{Functions to assess the relative importance of effects at observation moments}
\description{
 The function \code{sienaRI} returns the relative importance of effects of a SAOM according to the measure of relative importance described in
 Section 3.2 of Indlekofer and Brandes (2013). 
 The measure is based on the influence of effects on potential tie change or behavior change
 decisions of individual actors in network at the given observation moments. It takes the data as well as the complete model specification into account. Therefore, necessary arguments are the analysed data given as a \code{siena} data object as well as the complete model specification represented either by an estimated \code{sienaFit} object or by the triple consisting of a suitable parameter vector \code{theta} and the corresponding \code{sienaAlgorithm} and \code{sienaEffects} objects.
 }
\usage{
sienaRI(data, ans=NULL, theta=NULL, algorithm=NULL, effects=NULL)
\method{plot}{sienaRI}(x, file = NULL, col = NULL, addPieChart = FALSE,
              radius = 1, width = NULL, height = NULL, legend = TRUE, 
              legendColumns = NULL, legendHeight = NULL, 
              cex.legend = NULL, cex.names = NULL,\dots)
}
\arguments{
  \item{data}{\code{siena} data object representing the analyzed data and resulting from a call to \code{\link{sienaDataCreate}}.}
  \item{ans}{\code{\link{sienaFit}} object resulting from a call to \code{\link{siena07}}. The \code{sienaFit} object contains all necessary information on the model specification, in particular, the vector of parameter values \code{ans$theta}, the used algorithm for estimation \code{ans$x}, and information on included model effects \code{ans$effects}. If \code{ans} is a valid \code{sienaFit} object, the calculations of relative importances are based on \code{ans$theta}, \code{ans$x}, and \code{ans$effects}. Alternatively, the necessary information can be given directly as a suitable parameter vector \code{theta}, a \code{sienaAlgorithm} object, and a \code{sienaEffects} object. In this case, \code{ans} has to be unspecified (i.e., \code{ans=NULL}).}
  \item{theta}{Vector of parameter values of effects included in the model. Length of \code{theta} has to be equal to the number of included effects.}
  \item{algorithm}{\code{sienaAlgorithm} object as resulting from a call to \code{\link{sienaAlgorithmCreate}}. Works only for estimation by Method of Moments (i.e., if  \code{maxlike = FALSE}).}
  \item{effects}{\code{sienaEffects} object specifying which effects are included the model. Note that \code{sienaRI} does not yet work for endowment or creation effect, i.e., included effects have to be of type \code{eval} (or \code{rate}).}
  \item{x}{\code{sienaRI} object resulting from a call to \code{sienaRI}.}
  \item{file}{Name (path) of a pdf-file indicating whether the plot should be stored in a file. If \code{file=NULL}, the plot will be opened in a separate window. Otherwise a pdf with the given name will be created in the working directory or at the given path.}
  \item{col}{Colors used in the plot. If \code{col=NULL} a default color scheme is used.}
  \item{addPieChart}{Boolean: If \code{TRUE}, pie charts of aggregated relative importances for the complete network will be added to the plot.}
  \item{radius}{Radius of pie charts. Only effective if \code{addPieCharts = TRUE}.}
  \item{width}{Width of the plot. If \code{width=NULL} a default value is used.}
  \item{height}{Height of the plot. If \code{height=NULL} a default value is used.}
  \item{legend}{Boolean: if \code{TRUE} a legend is added to the plot. \code{x$effectNames} are used as labels.}
  \item{legendColumns}{Number of columns in legend. If \code{legendColumns=NULL} a default value is used. Only effective if \code{legend=TRUE}.}
  \item{legendHeight}{Height of legend. If \code{legendHeight=NULL} a default value is used. Only effective if \code{legend=TRUE}.}
  \item{cex.legend}{Specifies the relative font size of legend labels.}
  \item{cex.names}{Specifies the relative font size of bar graph labels.}
  \item{\dots}{Other arguments.}
}
\details{
\code{sienaRI} takes the data as well as the complete model specification into account. Therefore, necessary arguments are the analyzed data given as a \code{siena} data object as well as the complete model specification represented either by an estimated \code{sienaFit} object or by the triple consisting of a suitable parameter vector \code{theta} and the corresponding \code{sienaAlgorithm} and \code{sienaEffects} objects.

A \code{sienaFit} object contains all necessary information on the model specification, in particular, the vector of parameter values \code{ans$theta}, the used algorithm for estimation \code{ans$x}, and information on included model effects \code{ans$effects}.

If \code{ans} is a valid \code{sienaFit} object, the calculations of relative importances are based on \code{ans$theta}, \code{ans$x}, and \code{ans$effects}.\\ Alternatively, the necessary information can be given directly as a suitable parameter vector \code{theta}, a \code{sienaAlgorithm} object, and a \code{sienaEffects} object. In this case, \code{ans} has to be unspecified, i.e., \code{ans=NULL}.

Note that \code{sienaRI} works only with Method of Moments (i.e., for \code{sienaAlgorithm} objects with \code{maxlike = FALSE}) and that it does not yet work for endowment or creation effect (i.e., included effects have to be of type \code{eval}).
}
\value{
 If the model contains only one dependent variable, \code{sienaRI} returns an object of class \code{sienaRI}. Otherwise, it returns a list of objects of class \code{sienaRI}, each corresponding to one dependent variable.
 
 A returned \code{sienaRI} object stores the expected relative importances of effects of one dependent variable at observation moments as defined in Section 3.2 of Indlekofer and Brandes (2013).

A \code{sienaRI} object consists of following components:
    \describe{
    \item{- \code{dependentVariable}}{the name of the corresponding dependent variable.}
    \item{- \code{effectNames}}{the names of considered effects.}
    \item{- \code{RIActors}}{is a list that contains the expected relative importances of effects for each potential actor decision at observation moments.}
    \item{- \code{expectedRI}}{is a list that contains the expected relative importances of effects aggregated over all actors in each network observations. More precisely, it contains the averages of the actor related values in \code{RIActors}.}
    \item{- \code{absoluteSumActors}}{is a list that contains the sum of the (unstandardized) L1-differences calculated for each potential actor decision at observation moments. (Cf. Section 3.2 of Indlekofer and Brandes (2013))}
    \item{- \code{entropyActors}}{is a list that contains the entropy of each potential actor decision at observation moments.}
    }
}

\references{
\itemize{
Indlekofer, N. and Brandes, U., 
 \dQuote{Relative Importance of Effects in Stochastic Actor-oriented Models.}
 Network Science, 1 (3), 2013.
}
}
\author{Natalie Indlekofer}
\seealso{\code{sienaRIDynamics}}
\examples{
myalgorithm <- sienaAlgorithmCreate(nsub=2, n3=100)
mynet1 <- sienaDependent(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, density, recip, transTies, nbrDist2)
ans <- siena07(myalgorithm, data=mydata, effects=myeff, batch=TRUE)

RI <- sienaRI(mydata, ans)
RI
\dontrun{
plot(RI)
}

myalgorithm <- sienaAlgorithmCreate(nsub=2, n3=100)
mynet2 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mybeh <- sienaDependent(s50a, type="behavior")
mydata2 <- sienaDataCreate(mynet2, mybeh)
myeff2 <- getEffects(mydata2)
myeff2 <- includeEffects(myeff2, density, recip, transTies, nbrDist2)
ans2 <- siena07(myalgorithm, data=mydata2, effects=myeff2, batch=TRUE)
# Use only the parameters for the evaluation function:
theta.eval <- ans2$theta[myeff2$type[myeff2$include]=="eval"]
RI <- sienaRI(mydata2, theta=theta.eval, algorithm=myalgorithm,
              effects = myeff2)
RI
\dontrun{
plot(RI[[2]], col = c("red", "green"), legend=FALSE)
plot(RI[[1]], addPieChart = TRUE)
}
}
\keyword{models}