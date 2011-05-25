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

sienaGOF <- function(
		sienaFitObject, 
		auxiliaryFunction, 
		groupName=NULL, 
		varName=NULL, 
		wave=NULL,
		verbose=FALSE, join=TRUE,
		twoTailed=FALSE, 
		cluster=NULL, robust=FALSE, ...) 
	{
	require(MASS)
	#require(Matrix)
	##  Check input
	if (! sienaFitObject$returnDeps) 
	{
		stop("You must instruct siena07 to return the simulated networks")
	}
	iterations <- length(sienaFitObject$sims)
	if (iterations < 1) 
	{
		stop("You need at least one iteration.")
	}
	groups <- length(sienaFitObject$f$groupNames)
	if (verbose) cat("Detected", iterations, "iterations and", groups,
				"groups.\n")
	if (! is.null(groupName) ) {
		groupNumber <- match( groupName, sienaFitObject$f$groupNames)
	} 
	else 
	{
		groupName <- sienaFitObject$f$groupNames[1]
		groupNumber <- 1
	}
	if (is.na(groupNumber))
	{
		stop("Invalid group name.")
	}
	if (verbose) cat("Group", groupName, "corresponds to index",
				groupNumber,".\n")
	if(! is.null(varName)) {
		varNumber <- match( varName, sienaFitObject$f$depNames )
	} 
	else  
	{
		varNumber<-1
		varName <- sienaFitObject$f$depNames[1]
	}
	
	if (varNumber < 1 || varNumber > length(sienaFitObject$f$depNames))
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
		wave <- 1:(dim(sienaFitObject$f[[groupName]]$depvars[[varName]])[3] - 1)
	}
	if (varNumber < 1 || varNumber > 
			length(sienaFitObject$sims[[1]][[groupNumber]]))
	{
		stop("Invalid variable number -- out of bounds.")
	}
	if (min(wave) < 1 || max(wave) > 
			dim(sienaFitObject$f[[groupName]]$depvars[[varName]])[3] - 1
		)
	{
		stop("Invalid wave index -- out of bounds")
	}

	 obsStatsByWave <- lapply(wave, function (j) {
						matrix(
						auxiliaryFunction(NULL,
								sienaFitObject$f, sienaFitObject$sims,
								groupName, varName, wave), nrow=1) })
	if (join) 
	{
		obsStats <- Reduce("+", obsStatsByWave)
		obsStats <- list(Joint=obsStats)
	} 
	else  
	{
		obsStats <- obsStatsByWave
		names(obsStats) <- paste("Wave",wave)
	}
	class(obsStats) <- "observedAuxiliaryStatistics"
	attr(obsStats,"auxiliaryStatisticName") <-
			deparse(substitute(auxiliaryFunction))
	attr(obsStats,"joint") <- join
	
	##  Calculate the simulated auxiliary statistics
	if (verbose) cat("Calculating auxiliary statistics for waves",wave,".\n")
	
	if (!is.null(cluster)) {
		ttcSimulation <- system.time(simStatsByWave <- lapply(wave, 
			function (j) {
				simStatsByWave <- parSapply(cluster, 1:iterations, 
					function (i) 
						{ auxiliaryFunction(i,
									sienaFitObject$f, 
									sienaFitObject$sims,
									groupName, varName, wave)
							if (verbose && (i %% 100 == 0) ) 
								cat("  > Completed ", i,
									" calculations\n")
						})
							simStatsByWave <- matrix(simStatsByWave, 
								ncol=iterations)
							dimnames(simStatsByWave)[[2]] <-  1:iterations
							t(simStatsByWave)
						}))
	} 
	else  
	{
		ttcSimulation <- system.time(simStatsByWave <- lapply(wave,
						function (j) {
							simStatsByWave <- sapply(1:iterations, function (i) 
									{
											if (verbose && (i %% 100 == 0) ) 
											{
												cat("  > Completed ", i, 
													" calculations\n")
											}
											auxiliaryFunction(i,
													sienaFitObject$f, 
													sienaFitObject$sims,
													groupName, varName, wave)
									})
							simStatsByWave <-
									matrix(simStatsByWave, ncol=iterations)
							dimnames(simStatsByWave)[[2]] <-  1:iterations
							t(simStatsByWave)
						}))
	}
	
	## Aggregate by wave if necessary to produce simStats
	if (join) 
	{
		simStats <- Reduce("+", simStatsByWave)
		simStats <- list(Joint=simStats)
	} 
	else  
	{
		simStats <- simStatsByWave
		names(simStats) <- paste("Wave",wave)
	}
	class(simStats) <- "simulatedAuxiliaryStatistics"
	attr(simStats,"auxiliaryStatisticName") <-
			deparse(substitute(auxiliaryFunction))
	attr(simStats,"joint") <- join
	attr(simStats,"time") <- ttcSimulation
	
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
		observations <- nrow(observed)
		simulations<-nrow(simulated)
		variates<-ncol(simulated)
		if (robust) {
			a <- cov.rob(simulated)$cov
		} 
		else  
		{
			a <- cov(simulated)
		}
		ainv <- ginv(a)
		arank <- rankMatrix(a)
		expectation <- apply(simulated, 2, mean);
		centeredSimulations <- scale(simulated, scale=FALSE)
		if (variates==1) 
		{
			centeredSimulations <- t(centeredSimulations)
		}
		mhd <- function(x) 
		{
			x %*% ainv %*% x
		}
		simTestStat <- apply(centeredSimulations, 1, mhd)
		centeredObservations <- observed - expectation
		obsTestStat <- apply(centeredObservations, 1, mhd)
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
				InvCovSimStats=a,
				Rank=arank)
		class(ret) <- "sienaGofTest"
		attr(ret,"auxiliaryStatisticName") <- 
				attr(obsStats,"auxiliaryStatisticName")
		ret
	}

	res <- lapply(1:length(simStats),
					function (i) {
				 applyTest(obsStats[[i]], simStats[[i]]) })
	mhdTemplate <- rep(0, sum(sienaFitObject$test))
	names(mhdTemplate) <- rep(0, sum(sienaFitObject$test))
	if (join) {
		OneStepMHD <- list(mhdTemplate)
		PartialOneStepMHD <- list(mhdTemplate)
		GmmMhdValue <- list(mhdTemplate)
	} 
	else  
	{
		OneStepMHD <- lapply(wave, function(i) (mhdTemplate))
		PartialOneStepMHD <- lapply(wave, function(i) (mhdTemplate))
		GmmMhdValue <- lapply(wave, function(i) (mhdTemplate))
	}
	obsMhd <- NULL
	
	ExpStat <- lapply(wave, function(i) {
				apply(simStatsByWave[[i]], 2, mean)
			})
	OneStepSpecs <- matrix(0, ncol=sum(sienaFitObject$test), 
			nrow=length(sienaFitObject$theta))
	PartialOneStepSpecs <- matrix(0, ncol=sum(sienaFitObject$test), 
			nrow=length(sienaFitObject$theta))
	GmmOneStepSpecs <- matrix(0, ncol=sum(sienaFitObject$test), 
			nrow=length(sienaFitObject$theta))
	if (robust) {
		covInvByWave <- lapply(wave, function(i) ginv(
							cov.rob(simStatsByWave[[i]]) ))
	} 
	else  
	{
		covInvByWave <- lapply(wave, function(i) ginv( 
							cov(simStatsByWave[[i]]) ))
	}
	
	obsMhd <- sapply(wave, function (i) {
				 (obsStatsByWave[[i]] - ExpStat[[i]])  %*%
						covInvByWave[[i]] %*%
						t(obsStatsByWave[[i]] - ExpStat[[i]] )
			})

	if (sum(sienaFitObject$test) > 0) {
		effectsObject <- sienaFitObject$requestedEffects
		nSims <- sienaFitObject$Phase3nits
		if (join) {
			names(OneStepMHD[[1]]) <- 
					effectsObject$effectName[sienaFitObject$test]
			names(PartialOneStepMHD[[1]]) <- 
					effectsObject$effectName[sienaFitObject$test]
			names(GmmMhdValue[[1]]) <- 
					effectsObject$effectName[sienaFitObject$test]
		} 
		else  
		{
			for (i in wave) {
				names(OneStepMHD[[i]]) <- 
						effectsObject$effectName[sienaFitObject$test]
				names(PartialOneStepMHD[[i]]) <-
						effectsObject$effectName[sienaFitObject$test]
				names(GmmMhdValue[[i]]) <-
						effectsObject$effectName[sienaFitObject$test]
			}
		}
		rownames(OneStepSpecs) <- effectsObject$effectName
		colnames(OneStepSpecs) <- effectsObject$effectName[sienaFitObject$test]
		rownames(PartialOneStepSpecs) <- effectsObject$effectName
		colnames(PartialOneStepSpecs) <-
				effectsObject$effectName[sienaFitObject$test]
		rownames(GmmOneStepSpecs) <- effectsObject$effectName
		colnames(GmmOneStepSpecs) <-
				effectsObject$effectName[sienaFitObject$test]
		
		counterTestEffects <- 0
		for(index in which(sienaFitObject$test)) {
			if (verbose) {
				cat("Estimating test statistic for model including ",
						effectsObject$effectName[index], "\n")
			}
			counterTestEffects <- counterTestEffects + 1
			effectsToInclude <- !sienaFitObject$test
			effectsToInclude[index] <- TRUE
			theta0 <- sienaFitObject$theta
			names(theta0) <- effectsObject$effectName
			theta0 <- theta0[effectsToInclude]
			obsSuffStats <- 
					t(sienaFitObject$targets2[effectsToInclude, , drop=FALSE])
			G <- sienaFitObject$sf2[, , effectsToInclude, drop=FALSE] -
					rep(obsSuffStats, each=nSims)
			sigma <- cov(apply(G, c(1, 3), sum))
			SF <- sienaFitObject$ssc[ , , effectsToInclude, drop=FALSE]
			dimnames(SF)[[3]] <- effectsObject$effectName[effectsToInclude]
			dimnames(G) <- dimnames(SF)
			if (!(sienaFitObject$maxlike || sienaFitObject$FinDiff.method))
			{
				D <- RSienaTest:::derivativeFromScoresAndDeviations(SF, G)
			} 
			else  
			{
				DF <- sienaFitObject$
						sdf2[ , , effectsToInclude, effectsToInclude,
						drop=FALSE]
				D <- t(apply(DF, c(3, 4), mean))
			}
			fra <- apply(G, 3, sum) / nSims
			doTests <- rep(FALSE, sum(effectsToInclude))
			names(doTests) <- effectsObject$effectName[effectsToInclude]
			doTests[effectsObject$effectName[index]] <- TRUE
			mmThetaDelta <- as.numeric(
					RSienaTest:::ScoreTest(length(doTests), D,
							sigma, fra, doTests,
							maxlike=sienaFitObject$maxlike)$oneStep )
			mmPartialThetaDelta <- rep(0,length(theta0))
			mmPartialThetaDelta[length(theta0)] <-
					mmThetaDelta[length(theta0)]
			JacobianExpStat <- lapply(wave, function (i) {
				t(SF[,i,]) %*% simStatsByWave[[i]] / nSims
					})
			Gradient <- lapply(wave, function(i) {
						-2  * JacobianExpStat[[i]] %*%
								covInvByWave[[i]] %*% 
								t( obsStatsByWave[[i]] - ExpStat[[i]] ) })
			Hessian <- lapply(wave, function (i) { 
							2 *
							JacobianExpStat[[i]] %*%
							covInvByWave[[i]] %*%
							t(JacobianExpStat[[i]])
					})
			gmmThetaDelta <- -1 * as.numeric( ginv(Reduce("+", Hessian)) %*% 
							Reduce("+", Gradient) )
			OneStepSpecs[effectsToInclude,counterTestEffects] <- theta0 + 
					mmThetaDelta
			PartialOneStepSpecs[effectsToInclude,counterTestEffects] <- 
					theta0 + mmPartialThetaDelta
			GmmOneStepSpecs[effectsToInclude,counterTestEffects] <- theta0 + 
					gmmThetaDelta
			if (join) {
				OneStepMHD[[1]][counterTestEffects] <-
					as.numeric(
					sum(obsMhd) + 
					mmThetaDelta %*% Reduce("+", Gradient) + 0.5 *
					mmThetaDelta %*% Reduce("+", Hessian) %*% 
					mmThetaDelta)
				PartialOneStepMHD[[1]][counterTestEffects] <-
					as.numeric(
					sum(obsMhd) + 
					mmPartialThetaDelta %*% Reduce("+", Gradient) + 0.5 *
					mmPartialThetaDelta %*% Reduce("+", Hessian) %*% 
					mmPartialThetaDelta)
				GmmMhdValue[[1]][counterTestEffects] <-
					as.numeric( obsMhd + 
					gmmThetaDelta %*% 
					Reduce("+", Gradient) + 0.5 *
					gmmThetaDelta %*% 
					Reduce("+", Hessian) %*%
					gmmThetaDelta )
		} 
		else  
		{
				for (i in 1:length(obsMhd)) {
					OneStepMHD[[i]][counterTestEffects] <-  as.numeric(
						obsMhd[i] + 
						mmThetaDelta %*% Gradient[[i]] + 0.5 *
						mmThetaDelta %*% Hessian[[i]] %*% mmThetaDelta)
					GmmMhdValue[[i]][counterTestEffects] <-
							as.numeric( obsMhd[i] + 
							gmmThetaDelta %*% 
							Gradient[[i]] + 0.5 *
							gmmThetaDelta %*% 
							Hessian[[i]] %*%
							gmmThetaDelta )
					PartialOneStepMHD[[1]][counterTestEffects] <-
							as.numeric(
							sum(obsMhd) + 
							mmPartialThetaDelta %*% Reduce("+", Gradient) + 
							0.5 *
							mmPartialThetaDelta %*% Reduce("+", Hessian) %*% 
							mmPartialThetaDelta)
				}
			}
		}
	}

	names(res) <- names(obsStats)
	class(res) <- "sienaGOF"
	attr(res, "scoreTest") <- (sum(sienaFitObject$test) > 0)
	attr(res, "originalMahalanobisDistances") <- obsMhd
	attr(res, "oneStepMahalanobisDistances") <- OneStepMHD
	attr(res, "oneStepSpecs") <- OneStepSpecs
	attr(res, "partialOneStepMahalanobisDistances") <- PartialOneStepMHD
	attr(res, "partialOneStepSpecs") <- PartialOneStepSpecs
	attr(res, "gmmOneStepSpecs") <- GmmOneStepSpecs
	attr(res, "gmmOneStepMahalanobisDistances") <- GmmMhdValue
	attr(res,"auxiliaryStatisticName") <-
			attr(obsStats,"auxiliaryStatisticName")
	attr(res, "simTime") <- attr(simStats,"time")
	attr(res, "twoTailed") <- twoTailed
	attr(res, "joined") <- join
	res
}

print.sienaGOF <- function (x, ...) {
	## require(Matrix)
	levels <- 1:length(x)
	pVals <- sapply(levels, function(i) x[[i]]$p)
	titleStr <- "Monte Carlo Mahalanobis distance test P-value: "
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
			cat(" * Note for wave ",i,
				": Only ", x[[i]]$Rank, " statistics are ",
				"necessary in the auxiliary function.\n")
		}
	}
	else
	{
		cat(titleStr,pVals[1], "\n")
		cat("**Note: Only ", x[[1]]$Rank, " statistics are ",
				"necessary in the auxiliary function.\n")
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
	originalMhd <- attr(x, "originalMahalanobisDistances")
	if (attr(x, "joined")) {
		cat("-----\nCalculated joint MHD = (",
				sum(originalMhd),") for current model.\n")
	} 
	else  
	{
		for (j in 1:length(originalMhd)) {
			cat("-----\nCalculated wave ", j, " MHD = (",
					originalMhd[j],") for current model.")
		}
	}
	if (attr(x, "scoreTest")) {
		oneStepSpecs <- attr(x, "oneStepSpecs")
		oneStepMhd <- attr(x, "oneStepMahalanobisDistances")
		gmmMhd <- attr(x, "gmmOneStepMahalanobisDistances") 
		gmmOneStepSpecs <- attr(x, "gmmOneStepSpecs")
		partialOneStepSpecs <- attr(x, "partialOneStepSpecs")
		partialOneStepMhd <- attr(x, "partialOneStepMahalanobisDistances")
		
		if (attr(x, "joined")) {
			for (i in 1:ncol(oneStepSpecs)) {
				a <- cbind(oneStepSpecs[,i, drop=FALSE],
						partialOneStepSpecs[,i, drop=FALSE],
						gmmOneStepSpecs[,i, drop=FALSE] )
				b <- matrix( c(oneStepMhd[[1]][i], 
								partialOneStepMhd[[1]][i],
								gmmMhd[[1]][i]), ncol=3)
				rownames(b) <- c("MHD")
				a <- rbind(a, b)
				a <- round(a, 3)
				cat("\n**Model", colnames(a)[1], "\n")
				colnames(a) <- c("MM(Full)","MM(Par.)", "GMM")
				print(a)
			}
		} 
		else  
		{
			for (j in 1:length(oneStepMhd)) {
				for (i in 1:ncol(oneStepSpecs)) {
					a <- cbind(oneStepSpecs[,i, drop=FALSE],
							partialOneStepSpecs[,i, drop=FALSE],
							gmmOneStepSpecs[,i, drop=FALSE] )
					b <- matrix( c(oneStepMhd[[j]][i], 
									partialOneStepMhd[[j]][i],
									gmmMhd[[j]][i]), ncol=3)
					rownames(b) <- c("MHD")
					a <- rbind(a, b)
					a <- round(a, 3)
					cat("\n**Model", colnames(a)[1], "\n")
					colnames(a) <- c("MM(Full)","MM(Par.)", "GMM")
					print(a)
				}
			}
		}
		cat("\n-----")
	}
cat("\nComputation time for auxiliary statistic calculations on simulations: ",
			attr(x, "simTime")["elapsed"] , "seconds.")
}

plot.sienaGOF <- function (x, center=FALSE, scale=FALSE, violin=TRUE, 
		key=NULL, perc=.05, wave=1, main=main, ylab=ylab,  ...) 
{
	require(lattice)
	args <- list(...)
	if (is.null(args$main))
	{
		main=paste("Goodness of Fit of", 
				attr(x,"auxiliaryStatisticName"))
		if (!attr(x,"joined"))
		{
			main = paste(main, "Wave", wave)
		}
	} 
	else  
	{
		main=args
	}

	x <- x[[wave]]
	sims <- x$Simulations
	obs <- x$Observations
	itns <- nrow(sims)
	vars <- ncol(sims)
	## Need to check for useless statistics here:
	n.obs <- nrow(obs)
	
	if (center) 
	{
		sims.median <- apply(sims, 2, median)
		sims <- sapply(1:ncol(sims), function(i)
					(sims[,i] - sims.median[i]) )
		obs <- matrix(sapply(1:ncol(sims), function(i) 
							(obs[,i] - sims.median[i])), nrow=n.obs )
	} 
	if (scale) 
	{
		sims.min <- apply(sims, 2, min)
		sims.max <- apply(sims, 2, max)
		sims <- sapply(1:ncol(sims), function(i) sims[,i]/(sims.max[i] -
								sims.min[i] ) )
		obs <- matrix(sapply(1:ncol(sims), function(i) obs[,i] /(sims.max[i] - 
										sims.min[i] )
				), nrow=n.obs )
	} 

	if (is.null(args$ylab)) 
	{
		ylabel = "Statistic"
		if (center && scale) {
			ylabel = "Statistic (centered and scaled)"
		} 
		else if (scale) 
		{
			ylabel = "Statistic (scaled)"
		} 
		else if (center)
		{
			ylabel = "Statistic (center)"
		} 
		else  
		{
			ylabel = "Statistic"
		}
	} 
	else 
	{
		ylabel = args$ylab
	}
	
	screen <- sapply(1:ncol(obs),function(i){
						(sum(is.nan(rbind(sims,obs)[,i])) == 0) }) &
			(diag(var(rbind(sims,obs)))!=0)
	sims <- sims[,screen, drop=FALSE]
	obs <- obs[,screen, drop=FALSE]
	obsLabels <- round(x$Observations[,screen, drop=FALSE],3)
	key <- key[screen]

	if (is.null(args$xlab))
	{
		xlabel = paste( paste("p:", round(x$p, 3), 
						collapse = " "), collapse = "\n")
	} 
	else  
	{
		xlabel = args$xlab
	}

	xAxis <- (1:sum(screen))
	
	if (!is.null(key)) 
	{
		if (length(key) != ncol(obs)) 
		{
			stop("Key length does not match the number of variates.")
		}
	} 
	else 
	{
		key=xAxis
	}
	
	br <- trellis.par.get("box.rectangle")
	br$col <- 1
	trellis.par.set("box.rectangle", br)
	bu <- trellis.par.get("box.umbrella")
	bu$col <- 1
	trellis.par.set("box.umbrella", bu)
	plot.symbol <- trellis.par.get("plot.symbol")
	plot.symbol$col <- "black"
	plot.symbol$pch <- 4
	plot.symbol$cex <- 1
	trellis.par.set("plot.symbol", plot.symbol)
	
	panelFunction <- function(..., x=x, y=y, box.ratio){
		ind.lower = max( round(itns * perc/2), 1)
		ind.upper = round(itns * (1-perc/2))
		yperc.lower = sapply(1:ncol(sims), function(i)
					sort(sims[,i])[ind.lower]  )
		yperc.upper = sapply(1:ncol(sims), function(i)
					sort(sims[,i])[ind.upper]  )
		if (violin) {
			panel.violin(x, y, box.ratio=box.ratio, col = "transparent", ...)
		}
		panel.bwplot(x, y, box.ratio=.1, fill = "gray", ...)
		panel.xyplot(xAxis, yperc.lower, lty=3, col = "gray", lwd=3, type="l", 
				...)
		panel.xyplot(xAxis, yperc.upper, lty=3, col = "gray", lwd=3, type="l", 
				...)
		for(i in 1:nrow(obs))
		{
			panel.xyplot(xAxis, obs[i,],  col="red", type="l", lwd=1, ...)
			panel.xyplot(xAxis, obs[i,],  col="red", type="p", lwd=3, pch=19,
					...)
			panel.text(xAxis, obs[i,], labels=obsLabels[i,], pos=4)
		}
	}
	
	bwplot(as.numeric(sims)~rep(xAxis, each=itns), horizontal=FALSE,
			panel = panelFunction, xlab=xlabel, ylab=ylabel,
			scales=list(x=list(labels=key), y=list(draw=FALSE)),
			main=main, ...)
}

sparseMatrixExtraction <- function (i, data, sims, groupName, varName, wave) {
	#require(Matrix)
	dimsOfDepVar=dim(data[[groupName]]$depvars[[varName]][,,wave])
	missing <- Matrix(is.na(data[[groupName]]$depvars[[varName]][,,wave])*1)
	if (is.null(i)) {
		# sienaGOF wants the observation:
		returnValue <- Matrix(data[[groupName]]$depvars[[varName]][,,wave])
		returnValue[is.na(returnValue)] <- 0
	} 
	else  
	{
		#sienaGOF wants the i-th simulation:
		returnValue <- sparseMatrix(
				sims[[i]][[groupName]][[varName]][[wave]][,1],
				sims[[i]][[groupName]][[varName]][[wave]][,2],
				x=sims[[i]][[groupName]][[varName]][[wave]][,3],
				dims=dimsOfDepVar )
	}
	## Zero missings:
	1*((returnValue - missing) > 0)
}

snaEdgelistExtraction <- function (i, data, sims, groupName, varName, wave) {
	require(sna)
	returnValue <- snaSociomatrixExtraction(i, data, sims, groupName, varName,
			wave)
	as.edgelist.sna(returnValue)
}

snaSociomatrixExtraction <- function (i, data, sims, groupName, varName, wave) {
	require(sna)
	dimsOfDepVar=dim(data[[groupName]]$depvars[[varName]][,,wave])
	missing <- as.sociomatrix.sna(
			as.edgelist.sna(
					is.na(data[[groupName]]$depvars[[varName]][,,wave])))
	if (is.null(i)) {
		# sienaGOF wants the observation:
		returnValue <- data[[groupName]]$depvars[[varName]][,,wave]
		returnValue <- as.sociomatrix.sna(returnValue)
		returnValue[is.na(returnValue)] <- 0
	} 
	else  
	{
		#sienaGOF wants the i-th simulation:
		returnValue <- sims[[i]][[groupName]][[varName]][[wave]]
		attr(returnValue, "n") <- dimsOfDepVar[1]
		returnValue <- as.sociomatrix.sna( returnValue )
	}
	returnValue = 1*((returnValue - missing) > 0)
	returnValue
}

igraphEdgelistExtraction <- function (i, data, sims, groupName, varName, wave) {
	require(igraph)
	returnValue <- snaSociomatrixExtraction(i, data, sims, groupName, varName, 
			wave)
	graph.adjacency(returnValue)
}

OutdegreeDistribution <- function (i, data, sims, groupName, varName, wave,
		levels=0:8, extractor=snaSociomatrixExtraction) {
	x <- extractor(i, data, sims, groupName, varName, wave)
	a <- apply(x, 1, sum)
	sapply(levels, function(i){ sum(a<=i) })
}

IndegreeDistribution <- function (i, data, sims, groupName, varName, wave,
		levels=0:8, extractor=snaSociomatrixExtraction) {
	x <- extractor(i, data, sims, groupName, varName, wave)
	a <- apply(x, 2, sum)
	sapply(levels, function(i){ sum(a<=i) })
}

GeodesicDistribution <- function (i, data, sims, groupName, varName, wave,
		extractor=snaEdgelistExtraction, levels=1:8) {
	require(sna)
	x <- extractor(i, data, sims, groupName, varName, wave)
	a <- geodist(x)$gdist
	sapply(levels, function(i){ sum(a<=i) })
}

TriadCensus <- function (i, data, sims, groupName, varName, wave,
		extractor=snaEdgelistExtraction) {
	require(sna)
	x <- extractor(i, data, sims, groupName, varName, wave)
	triad.census(x)
}

KnnDistribution <- function (i, data, sims, groupName, varName, wave,
		extractor=igraphEdgelistExtraction, levels=0:25) {
	require(igraph)
	x <- extractor(i, data, sims, groupName, varName, wave)
	a <- graph.knn(x)$knn
	a[is.nan(a)] <- 0
	sapply(levels, function(i){ sum(a<=i) })
}
