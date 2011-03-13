## /*****************************************************************************
##  * SIENA: Simulation Investigation for Empirical Network Analysis
##  *
##  * Web: http://www.stats.ox.ac.uk/~snidjers/siena
##  *
##  * File: gof.r
##  *
##  * Description: This file contains the code to assess goodness of fit
##  *
##  ****************************************************************************/

sienaGOF <- function(sienaDataObject,
		sienaFitObject, groupName, varName,  auxiliaryFunction, wave=NULL,
		verbose=FALSE, join=TRUE, expectationFunction=mean,
		twoTailed=FALSE, cumulative=FALSE,  cluster=NULL, robust=FALSE, ...) 
	{

	require(MASS)
	##  require(Matrix)
	##  Check input
	if (! sienaFitObject$returnDeps) 
	{
		stop("You must instruct siena07 to return the simulated networks")
	}
	iterations = length(sienaFitObject$sims)
	if (iterations < 1) 
	{
		stop("You need at least one iteration.")
	}
	groups = length(sienaFitObject$sims[[1]])
	if (verbose) cat("Detected", iterations, "iterations and", groups,
				"groups.\n")
	groupNumber = which(names(sienaFitObject$sims[[1]]) == groupName)
	if (is.na(groupNumber))
	{
		stop("Invalid group name.")
	}
	if (verbose) cat("Group", groupName, "corresponds to index",
				groupNumber,".\n")
	varNumber = which(names(sienaFitObject$sims[[1]][[groupNumber]]) == 
					varName)
	if (varNumber < 1 | varNumber > length(sienaDataObject$depvars))
	{
		stop("Invalid variable number -- out of bounds.")
	}
	if (verbose)
	{
		cat("Variable", varName, 
				"corresponds to index",varNumber,".\n")
	}
	if (is.null(wave) )
	{
		wave = 1:(dim(sienaDataObject$depvars[[varNumber]])[3] - 1)
	}
	if (min(wave) < 1 | max(wave) > 
			dim(sienaDataObject$depvars[[varNumber]])[3] - 1 | 
			max(wave) > 
			length(sienaFitObject$sims[[1]][[groupNumber]][[varNumber]])
		)
	{
		stop("Invalid wave index -- out of bounds")
	}
	if (varNumber < 1 | varNumber > 
			length(sienaFitObject$sims[[1]][[groupNumber]]))
	{
		stop("Invalid variable number -- out of bounds.")
	}
	if (is.null(wave))
	{
		wave = 1:length(sienaFitObject$sims[[1]][[groupNumber]][[varNumber]])
	}
	
	##  Need to cast all of the observations and statistics into sparse matrices
	observedSp <- lapply(wave, function(j) { 
				Matrix(sienaDataObject$depvars[[varNumber]][,, j+1],
						sparse=TRUE)})
	dimsOfDepVar <- dim(observedSp[[1]])
	missingSp <- lapply(wave, function(j) { 
				Matrix(is.na(sienaDataObject$depvars[[varNumber]][,, j+1]),
						sparse=TRUE)})
	simsSp <- lapply(1:iterations, function(i){ lapply(wave, function(j) 
						{ 
				sparseMatrix(sienaFitObject$
					sims[[i]][[groupNumber]][[varNumber]][[j]][,1],
					sienaFitObject$
					sims[[i]][[groupNumber]][[varNumber]][[j]][,2],
					x=1,dims=dimsOfDepVar ) })})

	##  Now observation auxiliary statistics
	if (join) 
	{
		tmp <- lapply(wave, function (j) 
				{
					auxiliaryFunction(observedSp[[j]], missingSp[[j]]) })
		obsStats <- tmp[[1]]
		if (length(tmp)>1) {
			for (i in 2:length(tmp)) {
				obsStats = obsStats + tmp[[i]]
			}
		}
		ttcObservation <- system.time(obsStats <- list(Joint=obsStats))
	} else {
		ttcObservation <- system.time( obsStats <- lapply(wave, function (j) {
			auxiliaryFunction(observedSp[[j]], missingSp[[j]]) }))
		screen = lapply(wave, function (j) {
					is.na(sienaDataObject$depvars[[varNumber]][,,j+1])})
		names(obsStats) <- paste("Wave",wave)
	}
	class(obsStats) <- "observedAuxiliaryStatistics"
	attr(obsStats,"auxiliaryStatisticName") <-
			deparse(substitute(auxiliaryFunction))
	attr(obsStats,"joint") = join
	attr(obsStats,"time") = ttcObservation
	
	##  Calculate the simulated auxiliary statistics
	if (verbose) cat("Calculating auxiliary statistics for waves",wave,".\n")
	if (join) 
	{
		if (!is.null(cluster)) {
			ttcSimulation <- system.time(tmp <- lapply(wave, function (j) {
								tmp <- parSapply(cluster, 1:iterations, function (i) 
										{ auxiliaryFunction(simsSp[[i]][[j]], missingSp[[j]])})
								dimnames(tmp)[[2]] <-  1:iterations
								t(tmp)
							}))
		} else {
			ttcSimulation <- system.time(tmp <- lapply(wave, function (j) {
								tmp <- sapply(1:iterations, function (i) 
										{ auxiliaryFunction(simsSp[[i]][[j]], missingSp[[j]])})
								dimnames(tmp)[[2]] <-  1:iterations
								t(tmp)
							}))
		}

		simStats <- tmp[[1]]
		if (length(tmp)>1) {
			for (i in 2:length(tmp)) {
				simStats = simStats + tmp[[i]]
			}
		}
		simStats <- list(Joint=simStats)
	} else {
		if (!is.null(cluster)) {
			ttcSimulation <- system.time(simStats <- lapply(wave, function (j) {
								tmp <- parSapply(cluster, 1:iterations, function (i) {
											auxiliaryFunction(simsSp[[i]][[j]], missingSp[[j]])
										})
								dimnames(tmp)[[2]] <-  1:iterations
								t(tmp)
							}))
		} else {
			ttcSimulation <- system.time(simStats <- lapply(wave, function (j) {
								tmp <- sapply(1:iterations, function (i) {
											auxiliaryFunction(simsSp[[i]][[j]], missingSp[[j]])
										})
								dimnames(tmp)[[2]] <-  1:iterations
								t(tmp)
							}))
		}
		names(simStats) <- paste("Wave",wave)
	}
	class(simStats) <- "simulatedAuxiliaryStatistics"
	attr(simStats,"auxiliaryStatisticName") <-
			deparse(substitute(auxiliaryFunction))
	attr(simStats,"joint") = join
	attr(simStats,"time") = ttcSimulation
	
	applyTest <-  function (observed, simulated) 
	{
		if (class(simulated) != "matrix")
		{
			stop("Invalid input.")
		}
		if (class(observed) != "matrix") 
		{
			observed <- matrix(observed,nrow=1)
		}
		if (class(observed) != "matrix") 
		{
			stop("Observation must be a matrix.")
		}
		if (ncol(observed) != ncol(simulated)) 
		{
			stop("Dimensionality of function parameters do not match.")
		}
		observations = nrow(observed)
		simulations=nrow(simulated)
		variates=ncol(simulated)
		if (robust) {
			a <- cov.rob(simulated)$cov
		} else {
			a <- cov(simulated)
		}
		ainv <- ginv(a)
		arank <- rankMatrix(a)
		expectation = apply(simulated, 2, expectationFunction);
		centeredSimulations <- t(sapply(1:simulations, 
				function(i) simulated[i,] - expectation))
		if (variates==1) 
		{
			centeredSimulations <- t(centeredSimulations)
		}
		simTestStat <- sapply(1:simulations, function(i) 
						  centeredSimulations[i,] %*% 
						  ainv %*% centeredSimulations[i,])
		obsTestStat <- sapply(1:observations, function(i) 
					t(observed[i,] - expectation) %*%
							ainv %*% 
							(observed[i,] - expectation))
		if (twoTailed) 
		{
			p <- sapply(1:observations, function (i) 
						1 - abs(1 - 2 * sum(obsTestStat[i] <= 
						simTestStat)/length(simTestStat)) )
		} 
		else
		{
			p <- sapply(1:observations, function (i) 
				sum(obsTestStat[i] <= simTestStat) /length(simTestStat))
		}
		ret <- list( p = p,
				SimulatedTestStat=simTestStat,
				ObservedTestStat=obsTestStat,
				TwoTailed=twoTailed,
				Simulations=simulated,
				Observations=observed,
				Rank=arank)
		class(ret) <- "sienaGofTest"
		attr(ret,"auxiliaryStatisticName") <- 
				attr(obsStats,"auxiliaryStatisticName")
		ret
	}
	applyCumulativeTest <-  function (observed, simulated) 
	{
		if (class(simulated) != "matrix") 
		{
			stop("Invalid input.")
		}
		if (class(observed) != "matrix")
		{
			observed <- matrix(observed,nrow=1)
		}
		if (class(observed) != "matrix")
		{
			stop("Observation must be a matrix.")
		}
		if (ncol(observed) != ncol(simulated))
		{
			stop("Dimensionality of function parameters do not match.")
		}
		observations = nrow(observed)
		simulations=nrow(simulated)
		variates=ncol(simulated)
		nonnormalizedCdf = apply(simulated,2,mean)
		
		simTestStat <- sapply(1:simulations, function(i) 
					max(abs(simulated[i,] - nonnormalizedCdf)))
		obsTestStat <- sapply(1:observations, function(i) 
					max(abs(observed[i,] - nonnormalizedCdf)))
		if (twoTailed) 
		{
			p <- sapply(1:observations, function (i) 
					1 - abs(1 - 2 * sum(obsTestStat[i] <= 
					simTestStat)/length(simTestStat)) )
		} 
		else
		{
			p <- sapply(1:observations, function (i) 
					sum(obsTestStat[i] <= simTestStat) /length(simTestStat))
		}
		ret <- list( p = p,
				SimulatedTestStat=simTestStat,
				ObservedTestStat=obsTestStat,
				TwoTailed=twoTailed,
				Simulations=simulated,
				Observations=observed,
				Rank=NULL)
		class(ret) <- "sienaGofTest"
		attr(ret,"auxiliaryStatisticName") <- 
				attr(obsStats,"auxiliaryStatisticName")
		ret
	}
	if (cumulative)
	{
		testTime <- system.time(res <- lapply(1:length(simStats), function (i) 
					applyCumulativeTest(obsStats[[i]], simStats[[i]]) ))
	}
	else 
	{
		testTime <- system.time(res <- lapply(1:length(simStats), function (i) 
					 applyTest(obsStats[[i]], simStats[[i]]) ))
	}

	names(res) <- names(obsStats)
	class(res) <- "sienaGOF"
	attr(res, "originalSimulations") <- simsSp
	attr(res, "originalObservations") <- observedSp
	attr(res, "originalMissings") <- missingSp
	attr(res,"auxiliaryStatisticName") <-
			attr(obsStats,"auxiliaryStatisticName")
	attr(res, "obsTime") <- attr(obsStats,"time")
	attr(res, "simTime") <- attr(simStats,"time")
	attr(res, "testTime") <- testTime
	attr(res, "twoTailed") <- twoTailed
	attr(res, "joined") <- join
	attr(res, "cumulative") <- cumulative
	res
}

print.sienaGOF <- function (x, ...) {
	## require(Matrix)
	levels <- 1:length(x)
	pVals = sapply(levels, function(i) x[[i]]$p)
	if (attr(x, "cumulative"))
	{
		titleStr= "Monte Carlo Kolmogorov-Smirnov test P-value: "
	} 
	else
	{
		titleStr= "Monte Carlo Mahalanobis distance test P-value: "
	}
	cat("Siena Goodness of Fit (", 
			attr(x,"auxiliaryStatisticName") ,")\n=====\n")
	if (! attr(x,"join"))
	{
		cat(" >",titleStr, "\n")
		for (i in 1:length(pVals))
		{
			cat(" > Wave ", i, ": ", pVals[i], "\n")
		}
		for (i in 1:length(pVals)) 
		{
			if (! attr(x, "cumulative") &&
					x[[i]]$Rank != ncol(x[[i]]$Simulations)) 
			{
				cat(" * Note for wave ",i,
					": Only ", x[[i]]$Rank, " statistics are ",
					"necessary in the auxiliary function.\n")
			}	
		}
	}
	else
	{
		cat(titleStr,pVals[1], "\n")
		if (! attr(x, "cumulative") && 
				x[[1]]$Rank != ncol(x[[1]]$Simulations))
		{
			cat("**Note: Only ", x[[1]]$Rank, " statistics are ",
					"necessary in the auxiliary function.\n")
		}
	}
	
	if ( attr(x, "twoTailed") ) 
	{
		cat("-----\nTwo tailed test used.")	
	} 
	else 
	{
		cat("-----\nOne tailed test used ",
		"(i.e. area under curve for greater distance than observation).")	
	}
	
	cat("\nComputation time for auxiliary statistic calculations on observation: ",
			attr(x, "obsTime")["elapsed"] , "seconds.")
	
	cat("\nComputation time for auxiliary statistic calculations on simulations: ",
			attr(x, "simTime")["elapsed"] , "seconds.")
	
	cat("\nComputation time for distance/test statistic calculations: ",
			attr(x, "testTime")["elapsed"] , "seconds.\n")
}

dyadGOF <- function (x, wave, threshold, ...)
{
	## require(Matrix)
	a <- attr(x,"originalSimulations")
	b <- attr(x,"originalMissings")
	c <- attr(x,"originalObservations")
	truePositive <- a[[1]][[wave]] * 0
	trueNegative <- a[[1]][[wave]] * 0
	falsePositive <- a[[1]][[wave]] * 0
	falseNegative <- a[[1]][[wave]] * 0
	
	for(i in 1:length(a)) 
	{
		truePositive = truePositive + 
			a[[i]][[wave]] * c[[wave]]
		falseNegative = falseNegative + 
			(1-a[[i]][[wave]]) * c[[wave]]
		falsePositive = falsePositive +
			a[[i]][[wave]] * (1-c[[wave]])
		trueNegative = trueNegative +
			(1-a[[i]][[wave]]) * (1-c[[wave]])
	}
	threshold = length(a) * threshold
	
## 	trueNegative[trueNegative < threshold] <- 0
## 	trueNegative[trueNegative > 0] <- 1
## 	truePositive[truePositive < threshold] <- 0
## 	truePositive[truePositive > 0] <- 1
## 	falseNegative[falseNegative < (1-threshold)] <- 0
## 	falseNegative[falseNegative > 0] <- 1
## 	falsePositive[falsePositive < (1-threshold)] <- 0
## 	falsePositive[falsePositive > 0] <- 1

	falseNegative[falseNegative > threshold] <- 1
	falseNegative[falseNegative < 1] <- 0
	falsePositive[falsePositive > threshold] <- 1
	falsePositive[falsePositive < 1] <- 0
	
	plotted = as.matrix(-1*falsePositive + falseNegative)
	if (sum(plotted)==0) 
	{
		plotted <- matrix(0, nrow=nrow(plotted),ncol=ncol(plotted))
	}
	diag(plotted)
	image(z=plotted, zlim=c(-1,1), 
			col=c("red", "white", "blue"), ...)
}

plot.sienaGOF <- function (x, standardize=NA, violin=TRUE, 
			ylim=NULL, xlim=NULL, 
			xlab=NULL, ylab=NULL, key=NULL, 
			perc=.05, wave=1, main=NULL, 
			image=FALSE, imageThreshold=.3, ...) 
	{
	## require(Matrix)
	if (image) 
	{
		dyadGOF(x, wave, imageThreshold)
	}
	else 
	{
		if (is.null(main))
		{
			main=paste("Goodness of Fit of ", 
					attr(x,"auxiliaryStatisticName"))	
		}
		if (!attr(x,"joined"))
		{
			main = paste(main, "Wave", wave)
		}
		cumulative = attr(x,"cumulative")
		if (is.na(standardize)) 
		{
			if (attr(x,"cumulative")) 
			{
				standardize=0
			}
			else 
			{
				standardize=3
			}
		}
		x <- x[[wave]]
		sims <- x$Simulations
		obs <- x$Observations
		itns <- nrow(sims)
		vars <- ncol(sims)

		## Need to check for useless statistics here:
		n.obs <- nrow(obs)
		if (standardize==3) 
		{
			sims.median <- apply(sims, 2, median)
			sims.min <- apply(sims, 2, min)
			sims.max <- apply(sims, 2, max)
			sims <- sapply(1:ncol(sims), function(i)
						(sims[,i] - sims.median[i])/(sims.max[i] - 
									sims.min[i] ) )
			obs <- matrix(sapply(1:ncol(sims), function(i) 
								(obs[,i] - sims.median[i])/(sims.max[i] -
											sims.min[i] ) ), nrow=n.obs )
		} 
		else if (standardize==2) 
		{
			sims.min <- apply(sims, 2, min)
			sims.max <- apply(sims, 2, max)
			sims <- sapply(1:ncol(sims), function(i) (sims[,i] -
									sims.min[i])/(sims.max[i] - sims.min[i] ) )
			obs <- matrix(sapply(1:ncol(sims), function(i) (obs[,i] -
											sims.min[i])/(sims.max[i] - sims.min[i] )
					), nrow=n.obs )
		} 
		else if (standardize==1)
		{
			sims.mean <- apply(sims, 2, mean)
			sims <- sapply(1:ncol(sims), function(i)
						(sims[,i] - sims.mean[i]) )
			obs <- matrix(sapply(1:ncol(sims), function(i)
								(obs[,i] - sims.mean[i]) ), nrow=n.obs )
		}
		
		screen <- sapply(1:ncol(obs),function(i){
			(sum(is.nan(rbind(sims,obs)[,i])) == 0) }) &
			(diag(var(rbind(sims,obs)))!=0)
		sims <- sims[,screen, drop=FALSE]
		obs <- obs[,screen, drop=FALSE]
		obsLabels <- round(x$Observations[,screen, drop=FALSE],3)
		key <- key[screen]
		
		if (is.null(ylim)) 
		{
			ylim = c(min(obs, sims), max(obs, sims))
		}
		if (is.null(xlim)) 
		{
			xlim = c(0, ncol(obs)+1)
		}
		if (is.null(xlab))
		{
			xlab= paste( paste("p:", round(x$p, 3), 
							collapse = " "), collapse = "\n")
		}
		if (is.null(ylab)) 
		{
			ylab = "Statistic Values"
		}
		xAxis <- (1:sum(screen))
		
		plot(obs[1,]~xAxis, col="white", type="p",
				ylim=ylim, xlim=xlim, main=main,
				xlab=xlab, ylab=ylab, axes=FALSE, ...)
		if (!is.null(key)) 
		{
			if (length(key) != ncol(obs)) 
			{
				stop("Key length does not match the number of variates.")
			}
			axis(1, at=xAxis, lab=key)
		} 
		else 
		{
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
		
		if (violin) 
		{
			require(vioplot)
			for (i in 1:ncol(sims))
			{
				vioplot(sims[,i], at=xAxis[i],
						add=TRUE, col="gray", wex=.75, ...)
			}
		} 
		else
		{
			boxplot(as.numeric(sim)~rep(1:vars, each=itns), add=TRUE, ...)
		}
		for(i in 1:nrow(obs))
		{
			lines(obs[i,]~xAxis, col="red", type="l", lwd=1, ...)
			lines(obs[i,]~xAxis, col="red", type="p", lwd=3, pch=19, ...)
			text(xAxis, obs[i,], labels=obsLabels[i,], pos=4)
		}
	}
}