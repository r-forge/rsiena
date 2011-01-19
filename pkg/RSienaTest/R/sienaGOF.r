##/*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: http://www.stats.ox.ac.uk/~snidjers/siena
## *
## * File: gof.r
## *
## * Description: This file contains the code to assess goodness of fit
## *
## ****************************************************************************/

## Convenience function for evaluating statistics
evaluateAuxiliaryStatistics <- function(sienaObject, ...) {
	if (class(sienaObject) == "sienaFit") {
		evaluateAuxiliaryStatistics.simulated(sienaObject, ...)
	} else if (class(sienaObject) == "siena") {
		evaluateAuxiliaryStatistics.observed(sienaObject, ...)
	}
}

## Sienafit object, returns a list of auxiliary statistics from the given function
evaluateAuxiliaryStatistics.simulated <- function(sienaFitObject, 
		groupName, varName, wave, auxiliaryFunction, verbose=FALSE) {
	if (! sienaFitObject$returnDeps) {
		stop("You must instruct siena07 to return the simulated networks, i.e. returnDeps=TRUE")
	}
	iterations = length(sienaFitObject$sims)
	if (iterations < 1) {
		stop("You need at least one iteration.")
	}
	groups = length(sienaFitObject$sims[[1]])
	if (verbose) cat("Detected", iterations, "iterations and", groups, "groups.\n")
	groupNumber = which(names(sienaFitObject$sims[[1]]) == groupName)
	if (is.na(groupNumber)) {
		stop("Invalid group name.")
	}
	if (verbose) cat("Group", groupName, "corresponds to index",groupNumber,".\n")
	varNumber = which(names(sienaFitObject$sims[[1]][[groupNumber]]) == varName)
	if (varNumber < 1 | varNumber > length(sienaFitObject$sims[[1]][[groupNumber]])) {
		stop("Invalid variable number -- out of bounds.")
	}
	if (verbose) cat("Variable", varName, "corresponds to index",varNumber,".\n")
	if (min(wave) < 1 | max(wave) > length(sienaFitObject$sims[[1]][[groupNumber]][[varNumber]])){
		stop("Invalid wave index -- out of bounds")
	}
	if (verbose) cat("Calculating auxiliary statistics for waves",wave,".\n")
	ret <- lapply(wave, function (j) {
			tmp <- sapply(1:iterations, function (i) { auxiliaryFunction(sienaFitObject$sims[[i]][[groupNumber]][[varNumber]][[j]])})
			dimnames(tmp)[[2]] <-  1:iterations
			t(tmp)
		})
	names(ret) <- paste("Wave",wave)
	class(ret) <- "simulatedAuxiliaryStatistics"
	attr(ret,"auxiliaryStatisticName") <- deparse(substitute(auxiliaryFunction))
	ret
}

## Sienafit object, returns a list of auxiliary statistics from the given function
evaluateAuxiliaryStatistics.observed <- function(sienaDataObject, 
		varName, wave, auxiliaryFunction, verbose=FALSE) {
	varNumber = which(names(sienaDataObject$depvars) == varName)
	if (varNumber < 1 | varNumber > length(sienaDataObject$depvars)) {
		stop("Invalid variable number -- out of bounds.")
	}
	if (verbose) cat("Variable", varName, "corresponds to index",varNumber,".\n")
	if (min(wave) < 1 | max(wave) > dim(sienaDataObject$depvars[[varNumber]])[3] - 1){
		stop("Invalid wave index -- out of bounds")
	}
	if (verbose) cat("Calculating auxiliary statistics for waves",wave,".\n")
	ret <- lapply(wave, function (j) {
				auxiliaryFunction(sienaDataObject$depvars[[varNumber]][,,j+1])
			})
	names(ret) <- paste("Wave",wave)
	class(ret) <- "observedAuxiliaryStatistics"
	attr(ret,"auxiliaryStatisticName") <- deparse(substitute(auxiliaryFunction))
	ret
}

sienaGOF <- function(obsList, simList) {
	require(snopgof)
	if ((class(obsList) != "observedAuxiliaryStatistics") |
			(class(simList) != "simulatedAuxiliaryStatistics") |
			(length(obsList) != length(simList))
			| (attr(obsList,"auxiliaryStatisticName") != attr(simList,"auxiliaryStatisticName"))) {
		stop("Arguments invalid")
	}
	pre <- lapply(1:length(simList), function (i) gof.preprocess(simList[[i]]))
	# We could try inverse covariance for weighting or something here. Just identity weighting for now.
	#weight <- lapply(1:length(simList), function (i) solve(cov(simList[[i]])))
	#res <- lapply(1:length(simList), function (i) gof(obsList[[i]], pre[[i]], weight[[i]]))
	res <- lapply(1:length(simList), function (i) gof(obsList[[i]], pre[[i]]) )
	ret <- list(ObservedStatistics=obsList, SimulatedStatistics=simList, Results=res)
	class(ret) <- "sienaGOF"
	attr(ret,"auxiliaryStatisticName") <- attr(obsList,"auxiliaryStatisticName")
	ret
}

print.sienaGOF <- function (x, ...) {
	cat("Siena Goodness of Fit (", attr(x,"auxiliaryStatisticName") ,")\n > P-values by wave:\n")
	pVals = sapply(1:length(x$Results), function(i) x$Results[[i]]$p)
	for (i in 1:length(pVals)) {
		cat(" > Wave ", i, ": ", pVals[i], "\n")
	}
}

plot.sienaGOF <- function (x, wave=1, ...) {
	plot(x$Results[[wave]], main=paste("Goodness of Fit of ", attr(x, "auxiliaryStatisticName")),
			xlab="Variate Index", ylab="Variate Value", ...)
}