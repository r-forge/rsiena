\name{sienaCompositionChange}
\Rdversion{1.1}
\alias{sienaCompositionChange}
\alias{sienaCompositionChangeFromFile}
\title{
Functions to create a Siena composition change object
}
\description{Used to create a list of events describing the changes over
  time of a Siena actor set
}
\usage{
sienaCompositionChange(changelist, nodeSet = "Actors", option = 1)
sienaCompositionChangeFromFile(filename, nodeSet = "Actors",
    fileobj=NULL, option = 1)
}
\arguments{
  \item{changelist}{
    A list with an entry for each actor in the node set. Each entry a
    vector of numbers (may be as characters) indicating intervals during
    which the corresponding actor was present.
}
\item{filename}{
  Name of file containing change information. One line per actor, each
  line a series
  of space delimited numbers indicating intervals.
}
\item{fileobj}{
  The result of \command{readLines} on filename.
}
\item{nodeSet}{
  Character string containing the name of a Siena node set.
  If the entire data set contains more than one node set,
  then the node sets must be specified in all data objects.
}
\item{option}{
  Integer controlling the processing of the tie variables for the
  actors not currently present. Values (default is 1)
  \tabular{rl}{
1 \tab  0 before entry, final value carried forward after leaving\cr
2 \tab  0 before entry, missing after (final value carried forward, but
          treated as missing) \cr
3 \tab  missing whenever not in the network. Previous values will be used
        where available, \cr
  \tab  but always treated as missing values.\cr
4 \tab  Convert to structural zeros (not available at present).
}
}
}
\details{
  If there is a composition change object for the first node set in the
  data object, then this will be used in estimation by the Method of Moments
  to make actors active
  (able to send and receive ties) only for the time intervals
  when this is indicated in the composition change object.
  This is done according to the procedure of Huisman and Snijders (2003).
  See the manual for further details.\cr
  For bipartite networks, composition change objects for the second node set
  have no effect and will lead to an error message.\cr
  Intervals are treated as closed at each end.\cr
  For data sets including a composition change object, estimation
  by Method of Moments is forced to be unconditional, overriding
  the specification in the \code{\link{sienaAlgorithm}} object.
}
\value{
  An object of class "compositionChange", a list of numeric vectors,
  with attributes:
  \item{NodeSet}{Name of node set}
  \item{Option}{Option}
}
\references{\itemize{
   \item Huisman, M. E. and Snijders, T. A. B. (2003).
   Statistical analysis of longitudinal network data with changing composition.
   \emph{Sociological Methods & Research}, 32, 253-287.
   \item The manual: \url{http://www.stats.ox.ac.uk/~snijders/siena/RSiena_Manual.pdf}
   \item Further see \url{http://www.stats.ox.ac.uk/~snijders/siena/}
   }}
\author{Ruth Ripley }
\seealso{
 \code{\link{sienaNodeSet}},
 \code{\link{sienaDataCreate}}
}
\examples{
  clist <- list(c(1, 3), c(1.4, 2.5))
  #or
  clist <- list(c('1', '3'), c('1.4', '2.5'))

  compChange <- sienaCompositionChange(clist)

  \dontrun{
  filedata <- c("1 3", "1.4 2.5")
  write.table(filedata, "cc.dat",row.names=FALSE, col.names=FALSE,
          quote=FALSE)
  ## file will be
  ## 1 3
  ## 1.4 2.5
  compChange <- sienaCompositionChangeFromFile("cc.dat")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ classes }
