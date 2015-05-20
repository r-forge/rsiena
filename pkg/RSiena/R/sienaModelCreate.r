#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snijders/siena
# *
# * File: sienaModelCreate.r
# *
# * Description: This module contains the function for creating model objects.
# *
# *****************************************************************************/

ModelTypeStrings <- c("Standard actor-oriented model",
                      "Forcing model",
                      "Initiative model",
                      "Pairwise forcing model",
                      "Pairwise mutual model",
                      "Pairwise joint model")

##@sienaModelCreate DataCreate
sienaModelCreate <-
    function(fn,
             projname="Siena", MaxDegree=0, useStdInits=FALSE,
n3=1000, nsub=4, n2start = NULL, dolby=TRUE,
             maxlike=FALSE, diagonalize=1.0*!maxlike,
             condvarno=0, condname='',
             firstg=0.2, reduceg=0.5, cond=NA, findiff=FALSE,  seed=NULL,
             pridg=0.05, prcdg=0.05, prper=0.2, pripr=0.3, prdpr=0.3,
             prirms=0.05, prdrms=0.05, maximumPermutationLength=40,
             minimumPermutationLength=2, initialPermutationLength=20,
modelType=1, mult=5, simOnly=FALSE, localML=FALSE,
truncation=5, doubleAveraging=nsub, standardizeVar=(diagonalize<1))
{
    model <- NULL
    model$projname <- projname
    model$useStdInits <- useStdInits
    model$checktime <- TRUE
    model$n3 <- n3
    model$firstg <- firstg
    model$reduceg <- reduceg
    model$maxrat <- 1.0
    model$maxlike <-  maxlike
	model$simOnly <- simOnly
	model$localML <- localML
    model$FRANname <- deparse(substitute(fn))
    if (maxlike)
    {
        if (missing(fn))
        {
            model$FRANname <- "maxlikec"
        }
        if (is.na(cond))
        {
            cond <- FALSE
        }
        if (cond)
        {
            stop("Conditional estimation is not possible with",
                  "maximum likelihood estimation")
        }
        if (findiff)
        {
            stop("Finite differences estimation of derivatives",
                 "is not possible with maximum likelihood estimation")
        }
    }
    else
    {
        if (missing(fn))
        {
            model$FRANname <- "simstats0c"
        }
    }
    model$cconditional <- cond
    if (!is.na(cond) && cond && condvarno == 0 && condname == "")
    {
        model$condvarno <-  1
        model$condname <- ""
    }
    else
    {
        model$condvarno <-  condvarno
        model$condname <- condname
    }
    model$FinDiff.method <-  findiff
    model$nsub <- nsub
model$n2start <- n2start
	model$dolby <- (dolby && (!maxlike))
	if (diagonalize < 0) {diagonalize <- 0}
	if (diagonalize > 1) {diagonalize <- 1}
    model$diagg <- (diagonalize >= 0.9999)
	model$diagonalize <- diagonalize
    model$modelType <- modelType
    model$MaxDegree <- MaxDegree
    model$randomSeed <- seed
    model$pridg <- pridg
    model$prcdg <- prcdg
    model$prper <- prper
    model$pripr <- pripr
    model$prdpr <- prdpr
    model$prirms <- prirms
    model$prdrms <- prdrms
    model$maximumPermutationLength <- maximumPermutationLength
    model$minimumPermutationLength <- minimumPermutationLength
    model$initialPermutationLength <- initialPermutationLength
    model$mult <- mult
model$truncation <- truncation
model$doubleAveraging <- doubleAveraging
model$standardizeWithTruncation <- standardizeVar
model$standardizeVar <- standardizeVar
# The difference between these two is a hidden, non-documented option,
# perhaps for being tried out
# by later modification of the sienaAlgorithm object.
model$noAggregation <- FALSE
# This also is a hidden, non-documented option, perhaps for being tried out.
#  \item{noAggregation}{Logical:
#   do not replace current parameter value after subphase 1
#   by the mean over subphase 1, if some quasi-autocorrelation
#   then is larger than .5. May be helpful if initial value was very far away.
# The two options model$noAggregation and model$standardizeWithTruncation
# are used only in phase2.r.
model$browse1 <- FALSE # non-documented options for browsing in phase 2.
model$browse2 <- FALSE
model$browse3 <- FALSE
    class(model) <- "sienaAlgorithm"
    model
}

model.create <- sienaModelCreate


##@sienaAlgorithmCreate AlgoritmCreate
sienaAlgorithmCreate <- sienaModelCreate