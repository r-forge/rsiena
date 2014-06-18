##/*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: http://www.stats.ox.ac.uk/~snijders/siena
## *
## * File: bayesTest.r
## *
## * Description: This file contains the code to test results of sienaBayes,
## * and print results of the test.
## *
## ****************************************************************************/


##@simpleBayesTest Tests single parameters of sienaBayesFit objects
simpleBayesTest <- function(z, nfirst=z$nwarm+1, tested0=0,
							probs = c(0.025,0.975), ndigits=4){
	if (length(tested0) > 1)
	{
		stop('tested0 should be a number')
	}
	if (length(probs) != 2)
	{
		stop('probs should be a vector of length 2')
	}
    theEffects <- z$effects
	efNames <- format(paste(ifelse(theEffects$type == "creation",
									"creat", theEffects$type),
							theEffects$effectName)[!z$basicRate])
	if (all(theEffects$type[!z$basicRate] == 'eval'))
	{
		efNames <- format(theEffects$effectName[!z$basicRate])
	}
	else
	{
		efNames <- format(paste(ifelse(theEffects$type == "creation",
									"creat", theEffects$type),
							theEffects$effectName)[!z$basicRate])
	}
	credVal <- credValues(z, tested = tested0, theProbs = probs, nfirst=nfirst)
	mydf <- data.frame(matrix(NA, sum(!z$basicRate), 5))
	names(mydf) <- c(' ', 'varying', 'cred.from', '  cred.to', '  p  ')
	mydf[,1] <- efNames
	mydf[,2] <- ifelse(z$set2[!z$basicRate], '   -   ', '   +   ')
	mydf[,3] <- format(round(credVal[!z$basicRate, 1], digits=ndigits))
	mydf[,4] <- format(round(credVal[!z$basicRate, 2], digits=ndigits))
	mydf[,5] <- format(round(credVal[!z$basicRate, 3], digits=4))
	mydf
}


##@multipleBayesTest Tests single parameters of sienaBayesFit objects
multipleBayesTest <- function(z, ..., nfirst=z$nwarm, tested0=0, ndigits=4){
    theEffects <- z$effects
	efNames <- format(paste(ifelse(theEffects$type == "creation",
									"creat", theEffects$type),
							theEffects$effectName)[!z$basicRate])
	if (all(theEffects$type[!z$basicRate] == 'eval'))
	{
		efNames <- format(theEffects$effectName[!z$basicRate])
	}
	else
	{
		efNames <- format(paste(ifelse(theEffects$type == "creation",
									"creat", theEffects$type),
							theEffects$effectName)[!z$basicRate])
	}

	# 7 lines borrowed from sienaBayes code:
	vec1 <- 4 - theEffects$randomEffects[theEffects$include]
	vec1[z$basicRate] <- 2
	vec1[z$ratePositions[[1]]] <- 1
	vec1[z$effects$fix[z$effects$include]] <- 5
	# now 1=rates group 1; 2=other rates; 3=randomly varying;
	# 4 = estimated non-varying (eta); 5 = non-estimated, fixed.
	# set1 = (1,2,3) set2 = (4)
	ntot <- sum(!is.na(z$ThinPosteriorMu[,1]))
	if (nfirst >= ntot-1)
	{
		stop('Warm sample too short')
	}
	nmax <- ntot-nfirst+1
	z$ThinObjective <- matrix(NA, nmax,
		dim(z$ThinPosteriorEta)[2] +
			dim(z$ThinPosteriorMu[,z$objectiveInVarying, drop=FALSE])[2])
	vio <- vec1[vec1 %in% c(3,4)] == 3 # could be called z$varyingInObjective
	z$ThinObjective[1:nmax, vio] <-
			z$ThinPosteriorMu[nfirst:ntot, z$objectiveInVarying, drop=FALSE]
	z$ThinObjective[1:nmax, !vio] <- z$ThinPosteriorEta[nfirst:ntot,]
	p <- dim(z$ThinObjective)[2]
	testedSet <- c(...)
	if ((min(testedSet) < 1) | (max(testedSet) > p))
	{
		stop(paste('... should be a set of numbers between 1 and', p))
	}
	if (!(length(tested0) %in% c(1,length(testedSet))))
	{
		stop(paste('tested0 must be one number or a vector of',
					length(testedSet), 'numbers'))
	}
	if (length(tested0) == 1)
	{
		testedvec <- rep(tested0, length(testedSet))
	}
	else
	{
		testedvec <- tested0
	}
	efNames <- efNames[testedSet]
	posteriorSample <- z$ThinObjective[,testedSet, drop=FALSE]
	invCovParameters <- chol2inv(chol(cov(posteriorSample)))
	meanParameters <- mean(posteriorSample)
	standLength <- function(x){
		(x - meanParameters) %*% invCovParameters %*% (x - meanParameters)
	}
	posteriorStandLengths <- apply(posteriorSample, 1, standLength)
	boundary <- standLength(testedvec)
	postProb <- mean((posteriorStandLengths - boundary) > 0)
	result <- list(prob=postProb, chisquared=boundary,
				postDistances=posteriorStandLengths,
				nullValue=testedvec, effectNames=efNames)
	class(result) <- 'multipleBayesTest'
	result
}

print.multipleBayesTest <- function(x, ...){
	if (!inherits(x, 'multipleBayesTest'))
	{
        stop("not a legitimate multipleBayesTest object")
	}
	p <- length(x$effectNames)
	mydf <- as.data.frame(matrix(NA,p,2))
	mydf[,1] <- paste(x$effectNames,"= ")
	mydf[,2] <- x$nullValue
	names(mydf) <- c(' ', ' ')
	cat("Tested hypothesis:\n")
	print(mydf)
	cat("\nposterior p-value:", round(x$prob, 4), "\n")
	invisible(x)
}

stretch2 <- function(x){
	x[2] <- x[2] + 0.3*(x[2]-x[1])
	x
}

plot.multipleBayesTest <- function(x, xlim=NULL, ylim=NULL,
            main=NULL, ...){
	if (!inherits(x, 'multipleBayesTest'))
	{
        stop("not a legitimate multipleBayesTest object")
	}
# Makes a density plot of the posterior sample of distances.
	post <- x$postDistances
	par(oma=c(0,1,0,0), mar=c(5.1,5.1,4.1,2.1)) # to accommodate larger font for ylab
# density by group
	if (is.null(xlim)){xlim <- stretch2(range(post))}
	d1 <- density(post)
	if (is.null(ylim)){ylim <- stretch2(range(d1$y))}
	if (is.null(main)){main <- "posterior distances"}
	plot(d1, xlim=xlim, ylim=ylim, xlab='distance', ylab="density", col=4,
			cex=4, cex.lab=2, cex.main=2, main=main, lwd=2, ...)
}
