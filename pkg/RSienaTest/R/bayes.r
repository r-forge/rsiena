##/*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: http://www.stats.ox.ac.uk/~snijders/siena
## *
## * File: bayes.r
## *
## * Description: This file contains the code to run Bayesian simulation.
## * Many functions are defined within others to reduce copying of objects.
## *
## ****************************************************************************/

# see functions trafo, antitrafo, devtrafo; these should be consistent.
# antitrafo is inverse, devtrafo is derivative of trafo.

##@bayes Bayesian fit a Bayesian model, allowing a hierarchical structure
bayes <- function(data, effects, model, nwarm=100, nmain=100, nrunMHBatches=20,
				plotit=FALSE, nbrNodes=1,
				initgain=0.1, storeAll = FALSE,
				priorMu=NULL, priorSigma=NULL, priorDf=NULL, priorKappa=NULL,
				prevAns=NULL, clusterType=c("PSOCK", "FORK"),
				getDocumentation=FALSE)
{
    ##@createStores internal bayes Bayesian set up stores
    createStores <- function()
    {
        npar <- length(z$theta)

		numberRows <- nmain * nrunMHBatches
        z$posteriorTot <<- matrix(0, nrow=z$nGroup, ncol=npar)
        z$posteriorMII <<- array(0, dim=c(z$nGroup, npar, npar))
        z$candidates <<- array(NA, dim=c(numberRows, z$nGroup, npar))
        z$acceptances <<- matrix(NA, ncol=z$nGroup, nrow=numberRows)
		if (storeAll)
		{
			z$MHacceptances <<- array(NA, dim=c(numberRows, z$nGroup,
									 z$nDependentVariables, 9))
			z$MHrejections <<- array(NA, dim=c(numberRows, z$nGroup,
									 z$nDependentVariables, 9))
			z$MHproportions <<- array(NA, dim=c(numberRows, z$nGroup,
									 z$nDependentVariables, 9))
		}
		# johan's new stores with default values

		if (storeAll)
		{
			z$StorePosteriorMu <<- matrix(0, nrow=numberRows, ncol=z$TruNumPars)
			z$StorePosteriorSigma <<-
					array( NA, dim=c(numberRows, z$TruNumPars, z$TruNumPars) )
			# tom's additional stores
			z$StoreLoglik <<- rep(NA, numberRows)
		}
		# tom's additional thinned and warming stores
		z$ThinPosteriorMu <<- matrix(NA, nrow=nmain, ncol=z$TruNumPars)
		z$ThinPosteriorSigma <<-
					array( NA, dim=c(nmain, z$TruNumPars, z$TruNumPars) )
		z$ThinLoglik <<- rep(NA, nmain)
		z$ThinParameters <<- array(NA, dim=c(nmain, z$nGroup, z$TruNumPars))
		z$ThinBayesAcceptances <<- matrix(NA, nrow=nmain, ncol=z$nGroup)
		z$warmingMu <<- matrix(NA, nrow=nwarm, ncol=z$TruNumPars)
		z$warmingBasicRates <<- matrix(NA, nrow=nwarm, ncol=sum(z$basicRate))
	}

    ##@storeData internal bayes; put data in stores
    storeData <- function()
    {
        start <- z$sub + 1
        nrun <- nrow(z$parameters)
        end <- start + nrun - 1
        z$posteriorTot <<- z$posteriorTot + colSums(z$parameters)
		if (storeAll)
		{
			z$acceptances[start:end, ] <<- z$accepts
			z$candidates[start:end,, ] <<- z$parameters
			z$StorePosteriorMu[start:end, ] <<- z$posteriorMu
			z$StorePosteriorSigma[start:end, , ] <<-  z$PosteriorSigma
			z$StoreLoglik[start:end] <<- z$loglikVec
		}
		# Tom: Also store a thinned version of the process,
		# containing only the final values in the main steps:
        z$ThinPosteriorMu[z$iimain, ] <<- z$posteriorMu[nrun,]
		z$ThinPosteriorSigma[z$iimain, , ] <<-  z$PosteriorSigma[nrun,,]
        z$ThinLoglik[z$iimain] <<- z$loglikVec[nrun]
		z$ThinBayesAcceptances[z$iimain, ] <<- z$BayesAcceptances

        for (group in 1:z$nGroup)
        {
			# Drop those elements of the parameters array
			# that are NA because for a given group they refer
			# to rate parameters of other groups.
			z$ThinParameters[z$iimain, group, ] <<-
				z$parameters[nrun, group,!is.na(z$parameters[nrun, group, ])]
            for (i in dim(z$parameters)[1])
            {
                z$posteriorMII[group, , ] <<- z$posteriorMII[group, ,] +
                    outer(z$parameters[i, group, ], z$parameters[i, group, ])
            }
        }
		if (storeAll)
		{
			z$MHacceptances[start:end, , , ] <<- z$MHaccepts
			z$MHrejections[start:end, , , ] <<- z$MHrejects
			z$MHproportions[start:end, , , ] <<- z$MHaccepts /
				(z$MHaccepts + z$MHrejects)
        }
        z$sub <<- z$sub + nrun
    }

    ##@storeWarmingData internal bayes; put data in stores during warmup
    storeWarmingData <- function()
    {
		z$warmingMu[z$iiwarm,] <<- colMeans(z$posteriorMu)
        for (group in 1:z$nGroup)
        {
		# Stores only the first basic rate parameters
			groupStart <- (group-1)*length(z$ratePositions[[1]]) + 1
			groupEnd <- group*length(z$ratePositions[[1]])
			z$warmingBasicRates[z$iiwarm, groupStart:groupEnd] <<-
			    mean(z$parameters[, group, z$ratePositions[[group]]])
        }
    }

	##@improveMH internal bayes; find scale factors
    improveMH <- function(tiny=1.0e-15, desired=40, maxiter=100,
							tolerance=15, getDocumentation=FALSE)
    {
		##@rescaleCGD internal improveMH Bayesian
        rescaleCGD <- function(iter)
        {
			u <- ifelse (actual > desired,
							2 - ((iter - actual) / (iter - desired)),
							1 / (2 - (actual / desired)))
            number <<- ifelse(abs(actual - desired) <= tolerance,
								number + 1, 0 )
            success <<- number >= 2
            u
        }
		if (getDocumentation)
		{
			tt <- getInternals()
			return(tt)
		}
        iter <- 0
        number <- rep(0, z$nGroup)
        success <- rep(FALSE, z$nGroup)
		cat('improveMH\n')
        repeat
        {
            iter <- iter + 1
            MCMCcycle(nrunMH=1, nrunMHBatches=100, change=FALSE)
            actual <- z$BayesAcceptances ## acceptances
            ans <- rescaleCGD(100)
            update <- number < 3
            z$scaleFactors[update] <<- z$scaleFactors[update] * ans[update]
            cat(iter, '.', actual, '\n    ', ans, '\n')
			cat(' ', z$scaleFactors, '\n')
			flush.console()
            if (all(success) || iter == maxiter)
            {
                break
            }
            if (any(z$scaleFactors < tiny))
            {
                cat('scalefactor < tiny\n')
				flush.console()
                browser()
            }
        }
        cat('fine tuning took ', iter, ' iterations.\n')
		cat('Scalefactor:', z$scaleFactors, '\n')
		flush.console()
    }
	##@MCMCcycle internal bayes; some loops of (MH steps and sample parameters)
	MCMCcycle <- function(nrunMH, nrunMHBatches, change=TRUE)
	{
		z$accepts <<- matrix(NA, nrow=z$nGroup, nrunMHBatches)
		z$parameters <<- array(NA, dim=c(nrunMHBatches, z$nGroup, z$pp))
		z$MHaccepts <<- array(NA, dim=c(nrunMHBatches, z$nGroup,
								  z$nDependentVariables, 9))
		z$MHrejects <<- array(NA, dim=c(nrunMHBatches, z$nGroup,
								  z$nDependentVariables, 9))
		z$MHaborts <<- array(NA, dim=c(nrunMHBatches, z$nGroup,
								 z$nDependentVariables, 9))
		z$posteriorMu <<- matrix(0, nrow=nrunMHBatches, ncol=z$TruNumPars)
		z$PosteriorSigma <<- array( NA,
							dim=c(nrunMHBatches, z$TruNumPars, z$TruNumPars) )
		storeNrunMH <- z$nrunMH
		z$nrunMH <<- nrunMH
		z$loglikVec <<- rep(NA, nrunMHBatches)
		for (i in 1:nrunMHBatches)
		{
		#	cc <- proc.time()[1]

			ans <- z$FRAN(z, byGroup=TRUE, returnLoglik=TRUE, onlyLoglik=TRUE)
										#	c2 <- proc.time()[1]
										#	cat ('fran',c2-cc,'\n')
			z$loglik <<- ans$loglik
										#	cc <- proc.time()[1]
			sampleParameters(change)
										#	cc1 <- proc.time()[1]
										#	cat('samp',cc1-cc, '\n')
			# update grand means and covariance
			# do we need to check that there are more than one group?
			if ( z$nGroup > 1)
			{
				sampleMu()
				#print(z$muTemp)
				z$posteriorMu[i, ] <<- z$muTemp
				z$PosteriorSigma[i, ,] <<- z$SigmaTemp
			}
			z$loglikVec[i] <<- z$loglik
			# back to Ruth
			z$accepts[, i] <<- z$accept
			z$parameters[i, , ] <<- z$thetaMat
			z$MHaccepts[i, , , ] <<-
				t(do.call(cbind,
						  tapply(ans$accepts, factor(z$callGrid[, 1]),
								 function(x)Reduce("+", x))))
			z$MHrejects[i, , , ] <<-
				t(do.call(cbind, tapply(ans$rejects, factor(z$callGrid[, 1]),
										function(x)Reduce("+", x))))
			z$MHaborts[i, , , ] <<- t(do.call(cbind,
										tapply(ans$aborts,
												factor(z$callGrid[, 1]),
												function(x)Reduce("+", x))))
		}
		z$BayesAcceptances <<- rowSums(z$accepts)
		z$nrunMH <<- storeNrunMH
	}

	##@sampleParameters internal bayes; propose new parameters
	## and accept them or not
	## Replace accepted new parameters only if change.
	sampleParameters <- function(change=TRUE)
	{
		## get a multivariate normal with covariance matrix proposalCov
		## multiplied by a scale factor which varies between groups
		require(MASS)

	 	# do the fandango
		thetaChanges <- t(sapply(1:z$nGroup, function(i)
			{
				tmp <- z$thetaMat[i, ]
				use <- !is.na(z$thetaMat[i, ])
				tmp[use] <-
					mvrnorm(1, mu=rep(0, sum(use)),
							Sigma=z$scaleFactors[i] *
							z$proposalCov[[i]])
				tmp
			}
			))

		thetaOld <- z$thetaMat
		thetaOld[, z$basicRate] <- trafo(thetaOld[, z$basicRate])
		thetaNew <- thetaOld + thetaChanges
 		priorOld <- sapply(1:z$nGroup, function(i)
				{
				tmp <- thetaOld[i, ]
				   use <- !is.na(tmp)
				   dmvnorm(tmp[use],  mean=z$muTemp,
				   sigma=z$SigmaTemp)
				}
				)
		priorNew <- sapply(1:z$nGroup, function(i)
					{
					tmp <- thetaNew[i, ]
					use <- !is.na(tmp)
					dmvnorm(tmp[use],  mean=z$muTemp,
					sigma=z$SigmaTemp)
					}
					)

		logpOld <- z$loglik

# check that logpOld and getProbabilitiesFromC(z)[[1]]
# give same results; what is z[[1]]
		thetaNew[, z$basicRate] <- antitrafo(thetaNew[, z$basicRate])
		z$thetaMat <<- thetaNew
		logpNew <- getProbabilitiesFromC(z)[[1]]
		proposalProbability <- priorNew - priorOld + logpNew - logpOld
		##cat(proposalProbability, priorNew, priorOld, logpNew, logpOld, '\n')
		z$accept <<- log(runif(length(proposalProbability))) <
			proposalProbability
		thetaOld[, z$basicRate] <- antitrafo(thetaOld[, z$basicRate])
		if (!change)
		{
			z$thetaMat <<- thetaOld
		}
		else
		{
			##print(z$thetaMat)
			z$thetaMat[!z$accept, ] <<- thetaOld[!z$accept, ]
		}
##		print(thetaNew)
	}

	#### to draw the top-level means

	# johans new function for drawing thetas
	# tom changed riwish to rWishart - available from R 2.15.0
	##@sampleMu internal bayes; algorithm to propose new Mu parameter
	sampleMu <- function()
		{
		invWish <- function(v,S){
		# Draw from the inverse Wishart distribution with df v
		# and scale matrix S
		# inefficient if drawn multiply for the same S
			chol2inv(chol(rWishart(1,v,chol2inv(chol(S)))[,,1]))
		}
#		require(MCMCpack)
		Thetas <- matrix( c(0), z$TruNumPars, z$nGroup)
		muhat <- matrix( c(0), z$TruNumPars, 1)
		for ( groupiterator in c(1:z$nGroup) )
		{
			tmp <- z$thetaMat[groupiterator, ]
					use <- !is.na(tmp)
					Thetas[,groupiterator] <- tmp[use]
		}
		muhat <- rowMeans( Thetas )# the average across g
		Q <- matrix( c(0), z$TruNumPars, z$TruNumPars)

		for ( groupiterator in c(1:z$nGroup) )
		{
			Q <- Q + tcrossprod( Thetas[,groupiterator] - muhat)
		}
# riwish replaced by invWish
# prior Lambda = z$priorDf*z$priorSigma
		z$SigmaTemp <<- invWish(z$priorDf + z$nGroup,
				z$priorDf*z$priorSigma + Q +
				(z$priorKappa*z$nGroup/(z$priorKappa + z$nGroup))*
				tcrossprod( muhat - z$priorMu, muhat - z$priorMu ) )
		z$muTemp <<-  t(chol( ( z$priorKappa + z$nGroup )^(-1)*z$SigmaTemp )) %*%
			rnorm( z$TruNumPars , 0 , 1 ) +
			(z$nGroup/( z$priorKappa + z$nGroup ))*muhat +
			(z$priorKappa/( z$priorKappa + z$nGroup ))*z$priorMu
	}

	##@averageTheta internal bayes; algorithm to average past theta values
	averageTheta <- function(lukewarm = TRUE)
	{
		thetaMean <- rep(NA, z$pp)
		if (lukewarm)
		{
			start <- ceiling(dim(z$warmingMu)[1]*0.8)
			if (start < 10){start <- ceiling(dim(z$warmingMu)[1]*0.5)}
			end <- dim(z$warmingMu)[1]
			thetaMean[z$basicRate] <-
						colMeans(z$warmingBasicRates[start:end, , drop=FALSE])
			thetaMean[!z$basicRate] <- colMeans(
				z$warmingMu[start:end, z$generalParametersInGroup, drop=FALSE])
		}
		else
		{
			for (group in 1:z$nGroup)
			{
				thetaMean[z$ratePositions[[group]]] <- 	colMeans(
			z$ThinParameters[, group, !z$generalParametersInGroup, drop=FALSE])
			}
			thetaMean[!z$basicRate] <- colMeans(
				z$ThinPosteriorMu[,z$generalParametersInGroup, drop=FALSE])
		}
		thetaMean
	}

    ## ################################
    ## start of function bayes() proper
    ## ################################
	if (getDocumentation != FALSE)
	{
		if (getDocumentation == TRUE)
		{
			tt <- getInternals()
			return(tt)
		}
		else ## need to run getInternals on the argument value
		{
			targs <- formals(getDocumentation[1])
			targs[1:length(targs)] <- 1
			targs['getDocumentation'] <- TRUE
			if (length(getDocumentation) > 1)
			{
				targs['getDocumentation'] <- getDocumentation[-1]
			}
			return(do.call(getDocumentation[1], targs))
		}
	}
	cat('\n')
	ctime1 <- proc.time()[1]
	ctime <- proc.time()[1]

	z <- initializeBayes(data, effects, model, nbrNodes,
						initgain=initgain,
						priorMu=priorMu, priorSigma=priorSigma,
						priorDf=priorDf, priorKappa=priorKappa,
						prevAns=prevAns, clusterType=clusterType)

    createStores()
    z$sub <- 0

	ctime1 <- proc.time()[1]
	cat(ctime1-ctime,'\n')
	improveMH()
	ctime2<- proc.time()[1]

	cat('improveMH', ctime2-ctime1,'seconds.\n')
	flush.console()

    if (plotit)
    {
        require(lattice)
        dev.new()
        thetaplot = dev.cur()
        dev.new()
        ratesplot = dev.cur()
        dev.new()
        tseriesplot = dev.cur()
        dev.new()
        tseriesratesplot = dev.cur()
	}

    for (ii in 1:nwarm)
    {
		MCMCcycle(nrunMH=z$nrunMH, nrunMHBatches=nrunMHBatches)
		z$iiwarm <- ii
		storeWarmingData()
		cat('Warming step',ii,'(',nwarm,')\n')
		flush.console()
    }
	print('end of warming')
	ctime3<- proc.time()[1]

 	cat('warming took', ctime3-ctime2,'seconds.\n')
	flush.console()

	# Tom added averageTheta, and improveMH after warming up
	z$theta <- averageTheta(lukewarm = TRUE)
	cat('Parameter value after warming up\n')
	print(z$theta)
	cat('\n')
	ctime1 <- proc.time()[1]
	cat('Second ')
	improveMH()
	ctime2<- proc.time()[1]

	cat('Second improveMH', ctime2-ctime1,'seconds.\n')
	flush.console()
	ctime3 <- ctime2

    for (ii in 1:nmain)
    {
		MCMCcycle(nrunMH=z$nrunMH, nrunMHBatches=nrunMHBatches)
		z$iimain <- ii
		storeData()
		ctime4<- proc.time()[1]
		cat('main', ii, '(', nmain, ')', ctime4-ctime3, 'seconds.\n')
		flush.console()
		ctime3 <- ctime4
        if (ii %% 10 == 0 && plotit) ## do some plots
        {
            cat('main after ii', ii, '\n')
            dev.set(thetaplot)
            thetadf <-
                lapply(1:z$nGroup, function(i)
                   {
                       data.frame(Group=rep(i, ii * nrunMHBatches),
                                  z$candidates[1:(ii * nrunMHBatches), i, ])
                   }
                       )
            thetadf <- do.call(rbind, thetadf)
            basicRate <- z$basicRate
            ##thetadf <- data.frame(z$candidates)
            acceptsdf <- data.frame(z$MHproportions,
                                    z$acceptances)
            ratesdf <- thetadf[, -1, drop=FALSE][, z$basicRate, drop=FALSE]
            thetadf <- cbind(Group=thetadf[, 1, drop=FALSE],
							 thetadf[, -1, drop=FALSE][,
										   !z$basicRate, drop=FALSE])
            thetaNames<- paste(z$effects$name[!z$basicRate],
                               z$effects$shortName[!z$basicRate], sep=".")
            rateNames <- paste(z$effects$name[basicRate],
							   z$effects$shortName[basicRate],
							   z$effects$period[basicRate],
							   z$effects$group[basicRate], sep=".")
            names(ratesdf) <- rateNames
            ratesdf <- cbind(Group=thetadf[, 1, drop=FALSE], ratesdf)
            names(thetadf)[-1] <- make.names(thetaNames, unique=TRUE)
            names(acceptsdf) <- c("InsDiag", "CancDiag", "Permute", "InsPerm",
                                  "DelPerm", "InsMissing", "DelMissing",
                                  "BayesAccepts")
            varnames <- paste(names(thetadf)[-1], sep="", collapse= " + ")
			if (z$nGroup > 1)
			{
				varcall <- paste("~ ", varnames,  " | Group", sep="",
								 collapse="")
			}
			else
			{
				varcall <- paste("~ ", varnames,  sep="", collapse="")
			}
            print(histogram(as.formula(varcall), data=thetadf, scales="free",
                            outer=TRUE, breaks=NULL, type="density",
                            panel=function(x, ...)
                        {
                            panel.histogram(x, ...)
                            panel.densityplot(x, darg=list(na.rm=TRUE), ...)
                        }
                            ))
            dev.set(ratesplot)
            varnames <- paste(names(ratesdf)[-1], sep="", collapse= " + ")
            varcall <- paste("~ ", varnames, sep="", collapse="")
            print(histogram(as.formula(varcall), data=ratesdf, scales="free",
                            outer=TRUE, breaks=NULL, type="density",
                            panel=function(x, ...)
								{
								panel.histogram(x, ...)
								panel.densityplot(x, darg=list(na.rm=TRUE), ...)
								}
                            ))
            varnames <- paste(names(thetadf)[-1], sep="", collapse= " + ")
			if (z$nGroup > 1)
			{
				varcall <- paste(varnames,  "~ 1:", ii *
								 nrunMHBatches * z$nGroup,
                             " | Group", sep="", collapse="")
			}
			else
			{
				varcall <- paste(varnames,  "~ 1:", ii *
								 nrunMHBatches * z$nGroup,
								 sep="", collapse="")
			}
            dev.set(tseriesplot)
            print(xyplot(as.formula(varcall), data=thetadf, scales="free",
                         outer=TRUE))
            varnames <- paste(names(ratesdf)[-1], sep="", collapse= " + ")
            varcall <- paste(varnames,  "~ 1:", ii * nrunMHBatches * z$nGroup,
                             sep="", collapse="")
            dev.set(tseriesratesplot)
            print(xyplot(as.formula(varcall), data=ratesdf, scales="free",
                         outer=TRUE))
            ## dev.set(acceptsplot)
            ## varnames <- paste(names(acceptsdf), sep="", collapse= " + ")
            ## varcall <- paste("~ ", varnames,  sep="", collapse="")
            ## print(histogram(as.formula(varcall), data=acceptsdf,
            ##                 scales=list(x="same", y="free"),
            ##                 outer=TRUE, breaks=NULL, type="density",
            ##                 panel=function(x, ...)
            ##             {
            ##                 panel.histogram(x, ...)
            ##                 panel.densityplot(x, darg=list(na.rm=TRUE), ...)
            ##             }))
        }
    }
	cat('Total duration',ctime4-ctime,'seconds.\n')
	z$theta <- averageTheta(lukewarm = FALSE)
    z$FRAN <- NULL
    class(z) <- "sienaBayesFit"
	if (!storeAll)
	{
		z$acceptances <- NULL
		z$candidates <- NULL
	}
	z$OK <- TRUE
    z
}

##@initializeBayes algorithms do set up for Bayesian model
initializeBayes <- function(data, effects, model, nbrNodes, initgain,
						priorMu, priorSigma, priorDf, priorKappa,
						prevAns, clusterType=c("PSOCK", "FORK"))
{
	##@precision get inverse of z$covtheta, avoiding some inversion problems
	# for MoM estimates only
	precision <- function(z){
		require(MASS)
		t(z$dfrac) %*% ginv(z$msfc) %*% z$dfrac
	}

	##@projectEffects function to project model specification
	# modified version of updateTheta
	# model specification of prevAns is transferred to effects
	projectEffects <- function(effects, prevAns){
		prevEffects <- prevAns$requestedEffects
		prevEffects$initialValue <- prevAns$theta
		oldlist <- apply(prevEffects, 1, function(x)
						paste(x[c("name", "shortName",
								"type",
								"interaction1", "interaction2",
								"period")],
							collapse="|"))
		efflist <- apply(effects, 1, function(x)
						paste(x[c("name", "shortName",
								"type",
								"interaction1", "interaction2",
								"period")],
							collapse="|"))
		use <- efflist %in% oldlist
		effects$include[use] <-
			prevEffects$include[match(efflist, oldlist)][use]
		effects
	}
	## start initializeBayes
    ## initialize
    Report(openfiles=TRUE, type="n") #initialize with no file
    z  <-  NULL

    z$Phase <- 1
    z$Deriv <- FALSE
    z$FinDiff.method <- FALSE
    z$maxlike <- TRUE
    model$maxlike <- TRUE
    model$FRANname <- "maxlikec"
    z$print <- FALSE
    z$int <- 1
    z$int2 <- nbrNodes
    model$cconditional <-  FALSE
    if (!is.null(model$randomSeed))
    {
        set.seed(model$randomSeed)
		##seed <- model$randomSeed
    }
	else
	{
		if (exists(".Random.seed"))
		{
			rm(.Random.seed, pos=1)
		}
		newseed <- trunc(runif(1) * 1000000)
		set.seed(newseed)  ## get R to create a random number seed for me.
		##seed <- NULL
	}

	## Determine starting values and covariance matrices
	## for proposal distributions.
	# Starting values are obtained from MoM estimates.
	# First under the assumption that all parameters are the same.
	startupModel <- sienaModelCreate(n3=500, nsub=2, cond=FALSE)
	startupGlobal <- siena07(startupModel, data=data, effects=effects)
	cat("Initial global estimates\n")
	z$initialResults <- print(startupGlobal)
	cat('\n')
	effects <- updateTheta(effects, startupGlobal)
	z$initialGlobalEstimates <- effects[effects$include,]
	z$initialSE <- sqrt(diag(startupGlobal$covtheta))
	# Now per group.
	# Largest acceptable variance for proposal distribution rate parameters
	# on the trafo scale; note that this is before scaling by scaleFactors
	rateVarScale <- 0.25
	# weight for groupwise precision, relative to global precision:
	w <- 0.8
	# number of parameters per group, number of groups:
	npars <- sum(startupGlobal$effects$group==1)
	ngroups <- length(data)
	if ((startupGlobal$pp - npars)%%(ngroups-1) != 0)
	{
		stop("Not all groups have the same number of periods.")
	}
	# number of rate parameters per group:
	nrates <- (startupGlobal$pp - npars)%/%(ngroups-1)
	rateParameters <- matrix(NA, nrates, ngroups)
	# compute covariance matrices for random walk proposals for all groups:
	initialEstimates <- matrix(0, ngroups, startupGlobal$pp)
	z$proposalCov <- list()
	if (initgain <= 0)
	{
		startup1Model <- sienaModelCreate(n3=500, nsub=0, cond=FALSE)
	}
	else
	{
		startup1Model <- sienaModelCreate(n3=500, nsub=1, firstg=initgain,
										cond=FALSE)
	}
    for (i in 1:ngroups)
    {
		# coordinates of the global parameters used for this group:
		use <- rep(FALSE, startupGlobal$pp)
        use[!startupGlobal$effects$basicRate] <- TRUE
		rates <- !use
        use[startupGlobal$effects$group==i] <- TRUE
		# project model specification to the single group:
		rates <- rates[use]
		effects0 <- getEffects(data[[i]])
		effects0 <- projectEffects(effects0,startupGlobal)
		# get the parameter values from the global fit:
		effects0$initialValue[effects0$include] <-
			startupGlobal$theta[use]
		# estimate uncertainties for this group:
		startup1 <- siena07(startup1Model, data=data[[i]],
			effects=effects0)
		initialEstimates[i,use] <- startup1$theta
		rateParameters[,i] <- trafo(startup1$theta[rates])
		covmati <- chol2inv(chol(
			w*precision(startup1) +
			((1-w)/ngroups)*precision(startupGlobal)[use,use]))
		# now covmati is the roughly estimated covariance matrix
		# of the parameter estimate for group i.
		# Transform this to the trafo scale for the rate parameters:
        # double t() in view of recycling rule for matrix arithmetic
		covmati[,rates] <- t(t(covmati[,rates])*devtrafo(startup1$theta[rates]))
		covmati[rates,] <- covmati[rates,]*devtrafo(startup1$theta[rates])
		# truncate in case rate parameters have a too large variance
		largeRates <- rates & (diag(covmati) > rateVarScale)
		if (any(largeRates))
		{
			factor <- sqrt(rateVarScale/(diag(covmati)[largeRates]))
        # double t() in view of recycling rule for matrix arithmetic
			covmati[,largeRates] <- t(t(covmati[,largeRates])*factor)
			covmati[largeRates,] <- covmati[largeRates,]*factor
		}
		z$proposalCov[[i]] <- covmati
	}
	z$effectName <- effects0$effectName[effects0$include]
	z$requestedEffects <- startup1$requestedEffects
	z$initialGroupwiseEstimates <- initialEstimates

   	z$FRAN <- getFromNamespace(model$FRANname, pkgname)
    z <- initializeFRAN(z, model, data=data, effects=effects,
                prevAns=prevAns, initC=FALSE, onlyLoglik=TRUE)
	z$thetaMat <- initialEstimates
	z$basicRate <- z$effects$basicRate
    z$nGroup <- z$f$nGroup
	is.batch(TRUE)

	cat('Initial global parameter value ')
    print(z$theta)

    if (nbrNodes > 1 && z$observations > 1)
    {
        require(parallel)
		clusterType <- match.arg(clusterType)
		if (clusterType == "PSOCK")
		{
        clusterString <- rep("localhost", nbrNodes)
        z$cl <- makeCluster(clusterString, type = "PSOCK",
                            outfile = "cluster.out")
		}
		else
		{
			z$cl <- makeCluster(nbrNodes, type = "FORK",
								outfile = "cluster.out")
		}
        clusterCall(z$cl, library, pkgname, character.only = TRUE)
        clusterCall(z$cl, storeinFRANstore,  FRANstore())
        clusterCall(z$cl, FRANstore)
        clusterCall(z$cl, initializeFRAN, z, model,
                    initC = TRUE, profileData=FALSE, returnDeps=FALSE)
		clusterSetRNGStream(z$cl, iseed = as.integer(runif(1,
								max=.Machine$integer.max)))
    }

    z$scaleFactors <- rep(1, z$nGroup)
    ## z$returnDataFrame <- TRUE # chains come back as data frames not lists
    z$returnChains <- FALSE
  	groupPeriods <- attr(z$f, "groupPeriods")
    netnames <- z$f$depNames
	# Prepare some objects used for bookkeeping:
	z$rateParameterPosition <-
        lapply(1:z$nGroup, function(i, periods, data)
           {
               lapply(1:periods[i], function(j)
                  {
                      rateEffects <-
                          z$effects[z$effects$basicRate &
                                    z$effects$period == j &
                                    z$effects$group == i,]
                      rateEffects <-
                          rateEffects[match(netnames,
                                            rateEffects$name), ]
                      tmp <- as.numeric(row.names(rateEffects))
                      names(tmp) <- netnames
                      tmp
                  }
                      )
           }, periods=groupPeriods - 1, data=z$f[1:z$nGroup]
               )
	z$ratePositions <- lapply(z$rateParameterPosition, unlist)
	# Number of parameters in each group:
	z$TruNumPars <- sum( !z$basicRate ) + length(z$ratePositions[[1]])
	# Indicator of parameters in the parameter vector of length pp,
	# for which the hierarchical model applies,
	# which is the set of all parameters without the basic rate parameters:
	z$generalParameters <- (z$effects$group == 1) & (!z$basicRate)

	# Determination parameters of prior distribution;
	# this specifies the defaults!
	if (is.null(priorMu))
	{
		z$priorMu <- matrix( c(0), z$TruNumPars, 1 )
	}
	else
	{
		if (length(priorMu) != z$TruNumPars)
		{
			stop(paste("priorMu must have length ",z$TruNumPars))
		}
		z$priorMu <- matrix(priorMu, z$TruNumPars, 1 )
	}
	z$priorKappa <- 1
	if (!is.null(priorKappa))
	{
		if (priorKappa > 0)
		{
			z$priorKappa <- priorKappa
		}
	}
	z$priorDf <- z$TruNumPars + 2
	if (!is.null(priorDf))
	{
		if (priorDf > 0)
		{
			z$priorDf <- priorDf
		}
	}
	if (z$priorDf + z$nGroup <= z$TruNumPars)
	{
		z$priorDf <- z$TruNumPars - z$nGroup + 2
	}
	z$priorSigma <- diag(z$TruNumPars)
	if (!is.null(priorSigma))
	{
		if ((class(priorSigma)=="matrix") &
			all(dim(priorSigma)==c(z$TruNumPars,z$TruNumPars)))
		{
		z$priorSigma <- priorSigma
		}
		else
		{
		stop(paste("priorSigma is not a matrix of dimensions",z$TruNumPars))
		}
	}
	# Rate parameters are nuisance parameters, and get a data-induced prior;
	# observed covariance matrix on the trafo scale plus a small ridge.
	z$priorMu[rates] <- rowMeans(rateParameters)
	z$priorSigma[rates,] <- 0
	z$priorSigma[,rates] <- 0
	z$priorSigma[rates,rates] <- cov(t(rateParameters))
	diag(z$priorSigma[rates,rates]) <- diag(z$priorSigma[rates,rates]) + 0.1
	# Construction of indicator of parameters of the hierarchical model
	# within the groupwise parameters:
	vec <- rep(1,z$pp)
	vec[!z$basicRate] <- 2
	vec[z$ratePositions[[1]]] <- 3
	z$generalParametersInGroup <- vec[vec!= 1] == 2
	# Now generalParametersInGroup is a logical vector of length TruNumPars,
	# indicating the non-basicRate parameters.

	if (sum(z$basicRate) != z$nGroup*length(z$ratePositions[[1]]))
		{
		stop("Function bayes does not allow non-basic rate parameters.")
		}
	z$muTemp <- matrix(c(0), z$TruNumPars, 1)
	# Tom wonders about the deeper meaning of the next 10 or so lines.
	z$SigmaTemp <- matrix(c(0), z$TruNumPars, z$TruNumPars)
	z$SigmaTemp <- diag(z$TruNumPars)
    for (i in 1:z$nGroup)
    {
        use <- rep(FALSE, z$pp)
        use[z$ratePositions[[i]]] <- TRUE
        use[!z$basicRate] <- TRUE
        z$thetaMat[i, !use] <- NA
       # z$SigmaTemp <- z$SigmaTemp + z$proposalCov[use, use]
    }
    z
}

##@flattenChains algorithms converts a nested list of chains to a single list
flattenChains <- function(zz)
{
        for (i in 1:length(zz)) ##group
        {
            for (j in 1:length(zz[[i]])) ## period
            {
                attr(zz[[i]][[j]], "group") <- i
                attr(zz[[i]][[j]], "period") <- j
            }
        }
    zz <- do.call(c, zz)
    zz
}

##@dmvnorm algorithms calculates multivariate normal density
##inefficient: should not call mahalanobis and eigen with same sigma repeatedly
dmvnorm <- function(x, mean , sigma)
{
    if (is.vector(x))
    {
        x <- matrix(x, ncol=length(x))
    }
    distval <- mahalanobis(x, center=mean, cov=sigma)
    logdet <- sum(log(eigen(sigma, symmetric=TRUE, only.values=TRUE)$values))
    -(ncol(x) * log(2 * pi) + logdet + distval) / 2
}

##@getProbabilitiesFromC bayes gets loglik from chains in C
getProbabilitiesFromC <- function(z, index=1, getScores=FALSE)
{
	## expects maximum likelihood parallelisations
    f <- FRANstore()

	callGrid <- z$callGrid
    ## z$int2 is the number of processors if iterating by period, so 1 means
    ## we are not. Can only parallelize by period1
    if (nrow(callGrid) == 1)
    {
		theta <- z$thetaMat[1,]
        ans <- .Call("getChainProbabilities", PACKAGE = pkgname, f$pData,
					 f$pModel, as.integer(1), as.integer(1),
					 as.integer(index), f$myeffects, theta, getScores)
        anss <- list(ans)
	}
    else
    {
        if (z$int2 == 1 )
        {
            anss <- apply(callGrid, 1,
						  doGetProbabilitiesFromC, z$thetaMat, index, getScores)
        }
        else
        {
            use <- 1:(min(nrow(callGrid), z$int2))
            anss <- parRapply(z$cl[use], callGrid,
							  doGetProbabilitiesFromC, z$thetaMat, index,
							  getScores)
        }
    }
	ans <- list()
	ans[[1]] <- sum(sapply(anss, "[[", 1))
	if (getScores)
	{
		ans[[2]] <- rowSums(sapply(anss, "[[", 2))
	}
	ans[[3]] <- sapply(anss, "[[", 3)
	ans
}

##@doGetProbabilitiesFromC Maximum likelihood
doGetProbabilitiesFromC <- function(x, thetaMat, index, getScores)
{
    f <- FRANstore()
	theta <- thetaMat[x[1], ]
    .Call("getChainProbabilities", PACKAGE = pkgname, f$pData,
		  f$pModel, as.integer(x[1]), as.integer(x[2]),
		  as.integer(index), f$myeffects, theta, getScores)

}

##@trafo link function rates
# the log link is too strong!
trafo <- function(x){ifelse(x<0.01, 0.1, sqrt(x))}
##@antitrafo inverse link function rates
antitrafo <- function(x){x^2}
##@devtrafo derivative link function rates
devtrafo <- function(x){ifelse(x<0.01, 0.05, 1/(2*sqrt(x)))}
