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

## Sienafit object, returns a list of auxiliary 
## statistics from the given function
evaluateAuxiliaryStatistics.simulated <- function(sienaFitObject, 
		groupName, varName, auxiliaryFunction, wave=NULL,
		verbose=FALSE, join=TRUE) {
	if (! sienaFitObject$returnDeps) {
		stop("You must instruct siena07 to return the simulated networks")
	}
	iterations = length(sienaFitObject$sims)
	if (iterations < 1) {
		stop("You need at least one iteration.")
	}
	groups = length(sienaFitObject$sims[[1]])
	if (verbose) cat("Detected", iterations, "iterations and", groups,
				"groups.\n")
	groupNumber = which(names(sienaFitObject$sims[[1]]) == groupName)
	if (is.na(groupNumber)) {
		stop("Invalid group name.")
	}
	if (verbose) cat("Group", groupName, "corresponds to index",
				groupNumber,".\n")
	varNumber = which(names(sienaFitObject$sims[[1]][[groupNumber]]) == 
					varName)
	if (varNumber < 1 | varNumber > 
			length(sienaFitObject$sims[[1]][[groupNumber]])) {
		stop("Invalid variable number -- out of bounds.")
	}
	if (verbose) cat("Variable", varName, "corresponds to index",
				varNumber,".\n")
	if (is.null(wave)) {
		wave = 1:length(sienaFitObject$sims[[1]][[groupNumber]][[varNumber]])
	}
	if (min(wave) < 1 | max(wave) > 
			length(sienaFitObject$sims[[1]][[groupNumber]][[varNumber]])){
		stop("Invalid wave index -- out of bounds")
	}
	if (verbose) cat("Calculating auxiliary statistics for waves",wave,".\n")
	if (join) {
		ttc <- system.time(tmp <- lapply(wave, function (j) {
					tmp <- sapply(1:iterations, function (i) 
					{ auxiliaryFunction(sienaFitObject$
					sims[[i]][[groupNumber]][[varNumber]][[j]])})
					dimnames(tmp)[[2]] <-  1:iterations
					t(tmp)
				}))
		ret <- tmp[[1]]
		if (length(tmp)>1) {
			for (i in 2:length(tmp)) {
				ret = ret + tmp[[i]]
			}
		}
		ret <- list(Joint=ret)
	} else {
		ttc <- system.time(ret <- lapply(wave, function (j) {
					tmp <- sapply(1:iterations, function (i) {
					auxiliaryFunction(sienaFitObject$
					sims[[i]][[groupNumber]][[varNumber]][[j]])})
					dimnames(tmp)[[2]] <-  1:iterations
					t(tmp)
				}))
		names(ret) <- paste("Wave",wave)
	}
	class(ret) <- "simulatedAuxiliaryStatistics"
	attr(ret,"auxiliaryStatisticName") <- deparse(substitute(auxiliaryFunction))
	attr(ret,"joint") = join
	attr(ret,"time") = ttc
	ret
}

## Sienafit object, returns a list of auxiliary statistics from the given function
evaluateAuxiliaryStatistics.observed <- function(sienaDataObject, 
		varName, auxiliaryFunction, wave=NULL, verbose=FALSE, join=TRUE) {
	varNumber = which(names(sienaDataObject$depvars) == varName)
	if (varNumber < 1 | varNumber > length(sienaDataObject$depvars)) {
		stop("Invalid variable number -- out of bounds.")
	}
	if (verbose) cat("Variable", varName, 
				"corresponds to index",varNumber,".\n")
	if (is.null(wave) ) {
		wave = 1:(dim(sienaDataObject$depvars[[varNumber]])[3] - 1)
	}
	if (min(wave) < 1 | max(wave) > 
			dim(sienaDataObject$depvars[[varNumber]])[3] - 1){
		stop("Invalid wave index -- out of bounds")
	}
	if (verbose) cat("Calculating auxiliary statistics for waves",wave,".\n")
	if (join) {
		tmp <- lapply(wave, function (j) {
					auxiliaryFunction(sienaDataObject$
									depvars[[varNumber]][,,j+1])
				})
		ret <- tmp[[1]]
		if (length(tmp)>1) {
			for (i in 2:length(tmp)) {
				ret = ret + tmp[[i]]
			}
		}
		ret <- list(Joint=ret)
	} else {
		ret <- lapply(wave, function (j) {
					auxiliaryFunction(sienaDataObject$
									depvars[[varNumber]][,,j+1])
				})
		names(ret) <- paste("Wave",wave)
	}
	class(ret) <- "observedAuxiliaryStatistics"
	attr(ret,"auxiliaryStatisticName") <-
			deparse(substitute(auxiliaryFunction))
	attr(ret,"joint") = join
	ret
}


sienaGOFPreprocess <- function(simulated, expectationFunction=mean) {
	require(MASS)
	if (class(simulated) != "matrix") {
		stop("Invalid input.")
	}
	a <- cov(simulated)
	 # How do we make this error get thrown
	 #up the call stack? It will look much
	 # nicer if we can get that to work
	 # tryCatch(
	 # ainv <- solve(a)	
	 #, error = 
	 # function(e)simpleError(paste("Your testing functions are producing",
     # "collinear statistics which cannot be handled by the",
 	 # "Mahalanobis distance. Try reducing the number of statistics",
 	 # "returned by your testing function (e.g. Geodesic distances 1:10",
	 # "to 1:5 or so).")))
	# Using generalized inverse instead of regular inverse so
	# we don't have to worry
	# about collinearity
	ainv <- ginv(a)
	simulations=nrow(simulated)
	variates=ncol(simulated)
	expectation = apply(simulated, 2, expectationFunction);
	centeredSimulations <- t(sapply(1:simulations, 
			function(i) simulated[i,] - expectation))
	if (variates==1) {
		centeredSimulations <- t(centeredSimulations)
	}
	ttc <- system.time(mahalanobisDistances <- 
					sapply(1:simulations, function(i) 
					  centeredSimulations[i,] %*% 
					  ainv %*% centeredSimulations[i,]))
	ret <- list(
			Iterations=dim(simulated)[1], 
			Variates=dim(simulated)[2], 
			Original=simulated, 
			Covariance=a, InverseCovariance=ainv,
			Expectation=expectation,
			CenteredSimulations=centeredSimulations, 
			MHDistances = sort(mahalanobisDistances),
			ExpectationFunction=expectationFunction, 
			ComputeTime=ttc)
	class(ret) <- "sienaGOFPreprocess"
	ret
}

sienaGOFProcess <- function(observed, simulated, twoTailed=FALSE) {
	if (class(observed) == "numeric") {
		observed <- matrix(observed,nrow=1)
	}
	if (class(observed) != "matrix" |
			class(simulated) != "sienaGOFPreprocess") {
		stop("simulated must be a sienaGOFPreprocess object.")
	}
	if (ncol(observed) != simulated$Variates) {
		stop("Dimensionality of function parameters do not match.")
	}
	observations = nrow(observed)
	mhd <- sapply(1:observations, function(i) 
			  t(observed[i,] - simulated$Expectation) %*%
			  simulated$InverseCovariance %*% 
			  (observed[i,] - simulated$Expectation) )
	if (twoTailed) {
		p.mhd <- sapply(1:observations, function (i) 
			1 - abs(1 - 2 * sum(mhd[i] <= 
			simulated$MHDistances)/simulated$Iterations) )
	} else {
		p.mhd <- sapply(1:observations, function (i) 
			sum(mhd[i] <=  simulated$MHDistances)
			/simulated$Iterations)
	}
	return = list(
			Observation=observed, 
			Preprocess=simulated, 
			MHDistance=mhd, 
			p.MHD = p.mhd,
			TwoTailed = twoTailed)
	class(return) <- "sienaGOFProcess"
	return
}

sienaGOF <- function(obsList, simList) {
	if ((class(obsList) != "observedAuxiliaryStatistics") |
			(class(simList) != "simulatedAuxiliaryStatistics") |
			(length(obsList) != length(simList))
			| (attr(obsList,"auxiliaryStatisticName")
				!= attr(simList,"auxiliaryStatisticName"))) {
		stop("Arguments invalid")
	}

	pre <- lapply(1:length(simList), function (i)
				sienaGOFPreprocess(simList[[i]]))
	res <- lapply(1:length(simList), function (i) 
				sienaGOFProcess(obsList[[i]], pre[[i]]) )
	
	ret <- list(ObservedStatistics=obsList, SimulatedStatistics=simList, 
			Results=res)
	class(ret) <- "sienaGOF"
	attr(ret,"auxiliaryStatisticName") <-
			attr(obsList,"auxiliaryStatisticName")
	ret
}

print.sienaGOF <- function (x, ...) {
	pVals = sapply(1:length(x$Results), function(i) 
				x$Results[[i]]$p.MHD)
	ctimes = sapply(1:length(x$Results), function(i)
				as.numeric(x$Results[[i]]$Preprocess$ComputeTime["elapsed"]))
	cat("Siena Goodness of Fit (", 
			attr(x,"auxiliaryStatisticName") ,")\n=====\n")
	if (! attr(x$SimulatedStatistics,"joint")) {
		cat(" > Monte Carlo Mahalanobis distance test P-values:\n")
		for (i in 1:length(pVals)) {
			cat(" > Wave ", i, ": ", pVals[i], "\n")
		}
	} else {
		cat("Monte Carlo Mahalanobis distance test P-value: ",pVals[1], "\n")
	}
	if (x$Results[[1]]$TwoTailed) {
		cat("-----\nTwo tailed test used.")	
	} else {
		cat("-----\nOne tailed test used ",
		"(i.e. area under curve for greater distance than observation).")	
	}
	cat("\nBased on ", x$Results[[1]]$Preprocess$Iterations,
			" simulated draws from the null distribution.")
	cat("\nComputation time for auxiliary statistic calculations: ",
			attr(x$SimulatedStatistics,"time")["elapsed"] , "seconds.")
	cat("\nComputation time for distance/test statistic calculations: ",
			sum(ctimes), "seconds.\n")
}

plot.sienaGOF <- function (x, standardize=3, violin=TRUE, 
			ylim=NULL, xlim=NULL, 
			xlab=NULL, ylab=NULL, key=NULL, perc=.05, wave=1, main=NULL, ...) {
		if (is.null(main)) {
			main=paste("Goodness of Fit of ", 
					attr(x$Observed,"auxiliaryStatisticName"))	
		}
		if (!attr(x$Observed,"joint")){
			main = paste(main, "Wave", wave)
		}
		x <- x$Results[[wave]]
		itns <- x$Preprocess$Iterations
		vars <- x$Preprocess$Variates
		sims <- x$Preprocess$Original
		obs <- x$Observation
		n.obs <- nrow(obs)
		
		if (standardize==3) {
			sims.median <- apply(sims, 2, median)
			sims.min <- apply(sims, 2, min)
			sims.max <- apply(sims, 2, max)
			sims <- sapply(1:ncol(sims), function(i)
						(sims[,i] - sims.median[i])/(sims.max[i] - 
									sims.min[i] ) )
			obs <- matrix(sapply(1:ncol(sims), function(i) 
						(obs[,i] - sims.median[i])/(sims.max[i] -
									sims.min[i] ) ), nrow=n.obs )
		} else if (standardize==2) {
			sims.min <- apply(sims, 2, min)
			sims.max <- apply(sims, 2, max)
			sims <- sapply(1:ncol(sims), function(i) (sims[,i] -
						sims.min[i])/(sims.max[i] - sims.min[i] ) )
			obs <- matrix(sapply(1:ncol(sims), function(i) (obs[,i] -
						sims.min[i])/(sims.max[i] - sims.min[i] )
						), nrow=n.obs )
		} else if (standardize==1) {
			sims.mean <- apply(sims, 2, mean)
			sims <- sapply(1:ncol(sims), function(i)
						(sims[,i] - sims.mean[i]) )
			obs <- matrix(sapply(1:ncol(sims), function(i)
						(obs[,i] - sims.mean[i]) ), nrow=n.obs )
		}
		
		if (is.null(ylim)) {
			ylim = c(min(obs, sims), max(obs, sims))
		}
		if (is.null(xlim)) {
			xlim = c(0, ncol(obs)+1)
		}
		if (is.null(xlab)) {
			xlab= paste( paste("p:", round(x$p, 3), 
						collapse = " "), collapse = "\n")
		}
		if (is.null(ylab)) {
			ylab = "Statistic Values"
		}
		xAxis <- (1:vars)
		
		plot(obs[1,]~xAxis, col="white", type="p",
				ylim=ylim, xlim=xlim, main=main,
				xlab=xlab, ylab=ylab, axes=FALSE, ...)
		if (!is.null(key)) {
			if (length(key) != ncol(obs)) {
				stop("Key length does not match the number of variates.")
			}
			axis(1, at=xAxis, lab=key)
		} else {
			axis(1, at=xAxis, lab=paste("v", xAxis, sep=""))
		}
		
		ind.lower = round(itns * perc/2)
		ind.upper = round(itns * (1-perc/2))
		yperc.lower = sapply(1:ncol(sims), function(i)
					sort(sims[,i])[ind.lower]  )
		yperc.upper = sapply(1:ncol(sims), function(i)
					sort(sims[,i])[ind.upper]  )
		lines(yperc.lower~xAxis, lty=3, col = "gray", lwd=3)
		lines(yperc.upper~xAxis, lty=3, col = "gray", lwd=3)
		
		if (violin) {
			require(vioplot)
			for (i in 1:ncol(sims)) {
				vioplot(sims[,i], at=xAxis[i],
						add=TRUE, col="gray", wex=.75, ...)
			}
		} else {
			boxplot(as.numeric(sim)~rep(1:vars, each=itns), add=TRUE, ...)
		}
		for(i in 1:nrow(obs)) {
			lines(obs[i,]~xAxis, col="red", type="l", lwd=1, ...)
			lines(obs[i,]~xAxis, col="red", type="p", lwd=3, pch=19, ...)
			text(xAxis, obs[i,], labels=round(x$Observation[i,],3), pos=4)
		}
}