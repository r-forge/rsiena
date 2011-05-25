\name{sienaGOFsupplements}
\alias{GeodesicDistribution} 
\alias{IndegreeDistribution}
\alias{KnnDistribution}
\alias{OutdegreeDistribution}
\alias{TriadCensus}
\alias{igraphEdgelistExtraction}
\alias{snaEdgelistExtraction}
\alias{snaSociomatrixExtraction}
\alias{sparseMatrixExtraction}
\title{Supplemental functions to sienaGOF}
\description{
 These are supplemental functions for sienaGOF.
 }

\usage{
sparseMatrixExtraction(i, data, sims, groupName, varName, wave)

snaEdgelistExtraction(i, data, sims, groupName, varName, wave)

snaSociomatrixExtraction(i, data, sims, groupName, varName, wave)

igraphEdgelistExtraction(i, data, sims, groupName, varName, wave)

OutdegreeDistribution(i, data, sims, groupName, varName, wave,
		levels=0:8, extractor=snaSociomatrixExtraction)
		
IndegreeDistribution(i, data, sims, groupName, varName, wave,
		levels=0:8, extractor=snaSociomatrixExtraction) 
		
GeodesicDistribution(i, data, sims, groupName, varName, wave,
		extractor=snaEdgelistExtraction, levels=1:8)

TriadCensus(i, data, sims, groupName, varName, wave,
		extractor=snaEdgelistExtraction) 

KnnDistribution(i, data, sims, groupName, varName, wave,
		extractor=igraphEdgelistExtraction, levels=0:25) 
}
\arguments{
  \item{i}{ Iteration number. 0 if observation. }
  \item{data}{ The f object returned on the \code{siena07} object. }
  \item{sims}{ The sims object returned on the \code{siena07} object. }
  \item{groupName}{ The group name. }
  \item{varName}{ The name of the dependent variable. }
  \item{wave}{ The wave index. }
  \item{extractor}{ Which extractor function to use. }
  \item{levels}{ Which levels of auxiliary statistics. }
}
\details{
 See the manual for details on writing your own functions. Otherwise, you
 can use the following functions:
 
\code{ GeodesicDistribution(\dots)}

\code{ IndegreeDistribution(\dots)}

\code{ KnnDistribution(\dots)}

\code{ OutdegreeDistribution(\dots)}

\code{ TriadCensus(\dots)}
 
as the auxiliary functions for \code{sienaGOF}.

}

\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}
  for general information on RSiena and \url{http://www.stats.ox.ac.uk/~lospinos/}
 for information on the following presentation and working paper:
 
Lospinoso, J.A. and Snijders, T.A.B, "Goodness of fit for
Stochastic Actor Oriented Models." Presentation given at Sunbelt XXXI,
 St. Pete's Beach, Fl. 2011.

Lospinoso, J.A. and Snijders, T.A.B, "Goodness of fit for
Stochastic Actor Oriented Models." Working Paper.
}
\author{Josh Lospinoso}

\seealso{\code{\link{siena07}}, \code{\link{sienaGOF}} }

\keyword{models}
