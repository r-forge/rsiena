##/*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: http://www.stats.ox.ac.uk/~snijders/siena
## *
## * File: sienaBayes.r
## *
## * Description: This file contains the code to run Bayesian simulation.
## * Many functions are defined within others to reduce copying of objects.
## *
## ****************************************************************************/

# see functions trafo, antitrafo, devtrafo; these should be consistent.
# antitrafo is inverse, devtrafo is derivative of trafo.
# These functions also ensure positivity of the basic rate parameters.

##@sienaBayes Bayesian fit a Bayesian model, allowing a hierarchical structure
sienaBayes <- function(data, effects, algo, saveFreq=100,
				initgainGlobal=0.1, initgainGroupwise = 0.02, initfgain=0.2,
				gamma=0.05, initML=FALSE,
				priorMu=NULL, priorSigma=NULL, priorDf=NULL, priorKappa=NULL,
				priorRatesFromData=2,
				frequentist=FALSE, incidentalBasicRates=FALSE,
				reductionFactor=0.5,
				delta=1e-4,
				nwarm=50, nmain=250, nrunMHBatches=20,
				nSampVarying=1, nSampConst=1, nSampRates=0,
				nImproveMH=100, targetMHProb=0.25,
				lengthPhase1=round(nmain/5), lengthPhase3=round(nmain/5),
				storeAll = FALSE, prevAns=NULL, usePrevOnly=TRUE,
				prevBayes = NULL, newProposalFromPrev=TRUE,
				silentstart=TRUE,
				nbrNodes=1, clusterType=c("PSOCK", "FORK"),
				getDocumentation=FALSE)
{
    ##@createStores internal sienaBayes Bayesian set up stores
    createStores <- function()
    {
       # npar <- length(z$theta)
		numberRows <- nmain * nrunMHBatches
       # z$candidates <<- array(NA, dim=c(numberRows, z$nGroup, npar))
       # z$acceptances <<- matrix(NA, nrow=numberRows, ncol=z$nGroup)
		if (storeAll)
		{
			z$acceptances <<- matrix(NA, nrow=numberRows, ncol=z$nGroup)
			z$MHacceptances <<- array(NA, dim=c(numberRows, z$nGroup,
									 z$nDependentVariables, 9))
			z$MHrejections <<- array(NA, dim=c(numberRows, z$nGroup,
									 z$nDependentVariables, 9))
			z$MHproportions <<- array(NA, dim=c(numberRows, z$nGroup,
									 z$nDependentVariables, 9))
			z$StorePosteriorMu <<- matrix(0, nrow=numberRows, ncol=z$p1)
			z$StorePosteriorEta <<- matrix(0, nrow=numberRows, ncol=z$p2)
			z$StorePosteriorSigma <<-
					array( NA, dim=c(numberRows, z$p1, z$p1) )
		}
		# additional thinned stores
		z$ThinPosteriorMu <<- matrix(NA, nrow=nwarm+nmain, ncol=z$p1)
		z$ThinPosteriorEta <<- matrix(NA, nrow=nwarm+nmain, ncol=z$p2)
		z$ThinPosteriorSigma <<-
					array( NA, dim=c(nwarm+nmain, z$p1, z$p1) )
		z$ThinParameters <<- array(NA,
						dim=c(nwarm+nmain, z$nGroup, z$TruNumParsPlus))
		z$ThinBayesAcceptances <<-
						matrix(NA, nrow=nwarm+nmain, ncol=z$nGroup+2)
		z$sumBasicRates <<- matrix(0, nrow=z$nGroup, ncol=sum(z$basicRate))
		if (z$incidentalBasicRates)
		{
			z$ThinScores <<- matrix(NA, nrow=nwarm+nmain, ncol=sum(z$basicRate))
		}
#		z$ThinEtaScores <<- array(NA, dim=c(nwarm+nmain, z$p2, z$nGroup))
#		z$ThinEtaHessian <<- array(NA, dim=c(nwarm+nmain, z$p2, z$p2))
#	to be revived for frequentist:
#		z$ThinP3MuHat <<- matrix(NA, z$lengthPhase3, z$p1)
#		z$ThinP3SigmaHat <<- array(NA,
#							dim=c(z$lengthPhase3, z$p1, z$p1))
#		z$ThinP3EtaScores <<- matrix(NA, z$lengthPhase3, z$p2)
	} # end createStores

    ##@storeData internal sienaBayes; put data in stores
    storeData <- function()
    {
        start <- z$sub + 1
        nrun <- ncol(zm$accepts)
        end <- start + nrun - 1
		if (storeAll)
		{
			z$acceptances[start:end, ] <<- zm$accepts
			z$StorePosteriorMu[start:end, ] <<- zm$posteriorMu
			z$StorePosteriorEta[start:end, ] <<- zm$posteriorEta
			z$StorePosteriorSigma[start:end, , ] <<-  zm$PosteriorSigma
		}
		# Also store a thinned version of the process,
		# containing only the final values in the main steps:
        z$ThinPosteriorMu[z$iimain, ] <<- zm$posteriorMu[nrun,]
        z$ThinPosteriorEta[z$iimain, ] <<- zm$posteriorEta[nrun,]
		z$ThinPosteriorSigma[z$iimain, , ] <<-  zm$PosteriorSigma[nrun,,]
		z$ThinBayesAcceptances[z$iimain, ] <<- zm$BayesAcceptances
		if ((z$incidentalBasicRates) & (!is.null(zm$thescores)))
		{
			z$ThinScores[z$iimain, ] <<- zm$thescores
		}
        for (group in 1:z$nGroup)
        {
			# Drop those elements of the parameters array
			# that are NA because for a given group they refer
			# to rate parameters of other groups.
			z$ThinParameters[z$iimain, group, ] <<-
					z$thetaMat[group, !is.na(z$thetaMat[group, ])]
        }
		if (storeAll)
		{
			z$MHacceptances[start:end, , , ] <<- zm$MHaccepts
			z$MHrejections[start:end, , , ] <<- zm$MHrejects
			z$MHproportions[start:end, , , ] <<- zm$MHaccepts /
				(zm$MHaccepts + zm$MHrejects)
        }
# Dimension of the following wrong for more than 2 waves:
#		z$ThinEtaScores[z$iimain, ,] <<-
#				t(zm$ans.last$fra)[z$set2,,drop=FALSE]
#		z$ThinEtaHessian[z$iimain, ,] <<-
#				as.double(zm$ans.last$dff[z$set2,z$set2,drop=FALSE])
        z$sub <<- z$sub + nrun
    } # end storeData

	##@byM internal sienaBayes; for output of vector to screen in nice lines
	byM <- function(x,M=10)
	{
		m <- length(x)
		m1 <- ceiling(m/M)
		for (i in (1:m1))
		{
			cat(sprintf("%10.3f", x[((i-1)*M+1):min(m, i*M)]), "\n")
		}
		x
	}

	##@improveMH internal sienaBayes; find scale factors
    improveMH <- function(tiny=1.0e-15, totruns=100, target=targetMHProb,
					maxiter=20, tolerance=totruns/20, getDocumentation=FALSE)
    {
	# A rough stochastic approximation algorithm.
	# Steps when at least two iterations of "totruns" runs resulted,
	# for each coordinate of actual, in "desired" acceptances
	# with a tolerance of "tolerance".
		if (getDocumentation)
		{
			tt <- getInternals()
			return(tt)
		}
		desired <- trunc(target*totruns)
        iter <- 0
		nearGoal <- rep(FALSE, z$nGroup+2)
		farFromGoal <- rep(TRUE, z$nGroup+2)
		pastSuccesses <- rep(0, z$nGroup+2)
		groups <- c(rep(TRUE, z$nGroup), FALSE, FALSE)
		cat('improveMH\n')
		cat('Desired acceptances', desired, '.\n')
        repeat
        {
            iter <- iter + 1
            MCMCcycle(nrunMH=totruns, nSampVar=1, nSampCons=1,
							nSampRate=0, change=FALSE, bgain=-1)
            actual <- zm$BayesAcceptances
			## actual is a vector of the expected number of acceptances by group
			## and for the non-varying parameters, in the MH change of
			## parameters out of nrunMHBatches trials
			if (!(any(z$set2) & (!frequentist)))
			{
				actual[(z$nGroup+1) : (z$nGroup+2)] <- desired
			}
			if (any(is.na(actual)))
			{
				cat('is.na(actual) for ',which(is.na(actual)),"\n")
cat("Please report to Tom; and hit a key to continue\n")
# if this never happens this check can be dropped
browser()
				actual[is.na(actual)] <- actual.previous
				# this is OK unless it is the first time,
				# when actual.previous is still undefined
			}
			# if nearGoal (already at the previous step),
			# replace actual by exponentially weighted moving average of actual
			actual <- ifelse(nearGoal, (actual + actual.previous)/2, actual)
			# now update nearGoal
			nearGoal <- ifelse((abs(actual - desired) < tolerance),
								TRUE, nearGoal)
			farFromGoal <- ifelse((abs(actual - desired) < 2*tolerance),
								FALSE, farFromGoal)
			farMult <- ifelse(actual > desired, 2, 0.5)
			farMult <- ifelse(actual > 5*desired, 5, farMult)
			farMult <- ifelse(5*actual < desired, 0.2, farMult)
			actual.previous <- actual
            pastSuccesses <- ifelse(abs(actual - desired) <= tolerance,
									pastSuccesses+1, pastSuccesses)
# The calculation of the update:
			multiplier <- ifelse(farFromGoal, farMult,
				1 + (actual - desired)/
								(sqrt(sqrt(iter)*desired*(totruns - desired))))
			multiplier <- pmin(pmax(multiplier, 0.1), 10)
            z$scaleFactors <<-
					z$scaleFactors * multiplier[groups]
			z$scaleFactor0 <<-
					z$scaleFactor0 * multiplier[z$nGroup+1]
			z$scaleFactorEta <<-
					z$scaleFactorEta * multiplier[z$nGroup+2]
			if (any(is.na(z$scaleFactors)))
			{
				cat('\nany(is.na(z$scaleFactors)) in improveMH\n')
				browser()
			}
            cat('\n',iter, '.         ', sprintf("%5.1f", actual),
				     '\n multipliers  ', sprintf("%5.3f", multiplier),
			         '\n scaleFactors ', sprintf("%5.3f", z$scaleFactors),
					sprintf("%5.3f", z$scaleFactor0),
					sprintf("%5.3f", z$scaleFactorEta),'\n')

			flush.console()
            if ((min(pastSuccesses) >= 2) || (iter == maxiter))
            {
                break
            }
            if (any(z$scaleFactors < tiny))
            {
				cat('tiny: ', which(z$scaleFactors < tiny),'\n')
				cat('This is a problem: scaleFactors too small.\n')
				cat('Try a model with fewer randomly varying effects.\n')
				stop('Tiny scaleFactors.')
            }
        }
        cat('fine tuning took ', iter, ' iterations.\n')
#drop		cat('scaleFactors:', z$scaleFactors, z$scaleFactor0,
#										z$scaleFactorEta, '\n')
		flush.console()
    } # end improveMH


	##@MCMCcycle internal sienaBayes; some loops of
	## (MH steps and sample parameters)
	MCMCcycle <- function(nrunMH, nSampVar, nSampCons, nSampRate,
							change=TRUE, bgain)
	{
	# is for the storage in the object produced
	# and for its by-effects on the MCMC state:
	# thetaMat, muTemp, sigmaTemp
		zm <<- list()
		zm$accepts <<- matrix(NA, nrow=z$nGroup, nrunMH)
		acceptsEta <<- rep(NA, nrunMH)
		zm$MHaccepts <<- array(NA, dim=c(nrunMH, z$nGroup,
							  z$nDependentVariables, 9))
		zm$MHrejects <<- array(NA, dim=c(nrunMH, z$nGroup,
							  z$nDependentVariables, 9))
		zm$MHaborts <<- array(NA, dim=c(nrunMH, z$nGroup,
							 z$nDependentVariables, 9))
		zm$posteriorMu <<- matrix(0, nrow=nrunMH, ncol=z$p1)
		zm$posteriorEta <<- matrix(0, nrow=nrunMH, ncol=z$p2)
		zm$PosteriorSigma <<- array( NA, dim=c(nrunMH, z$p1, z$p1))
		zm$BayesAcceptances <<- rep(NA, z$nGroup+2)
        zsmall <<- getFromNamespace("makeZsmall", pkgname)(z)

		for (i in 1:nrunMH)
		{
		# This is the hot loop.
		# sample MH steps in the network chain:
			if (i%%10 == 0)
			{
				cat(".") # show something is happening
				flush.console()
			}
			if (i < nrunMH)
			{
				ans <- z$FRAN(zsmall, byGroup=TRUE, returnLoglik=FALSE)
			}
			else
			{#abcd
				zsmall$Deriv <<- TRUE
				ans <- z$FRAN(zsmall, byGroup=TRUE, returnLoglik=TRUE, onlyLoglik=FALSE)
				zsmall$Deriv <<- FALSE
			}
			# After the nrunMH the loglik by group is used as logpOld,
			# see below.
			if (z$p1 > 0)
			{
			# sample the group-level parameters:
			# First prepare eigenvalues and inverse of SigmaTemp
			# for use in dmvnorm2 in sampleVaryingParameters
			z$SigmaTemp <<- correctMatrix(z$SigmaTemp, z$delta)
			z$SigmaTemp.evs <<-
				eigen(z$SigmaTemp, symmetric=TRUE, only.values=TRUE)$values
			if (inherits(try(z$SigmaTemp.inv <<- chol2inv(chol(z$SigmaTemp)),
					silent=TRUE), "try-error"))
			{ # practically impossible given correctMatrix
				cat('Covariance matrix z$SigmaTemp is non-invertible;\n')
					if (nmain <= z$p1 + 2)
				{
					cat('this may be because nmain is too small;\n or')
				}
				cat('perhaps the model is not identified; please check.\n')
				stop('Inversion error in Cholesky decomposition.')
			}

			for (i2 in 1:nSampVar)
			{
				mcmcc.accept <- sampleVaryingParameters(change, bgain)
			}
			# extra sampling of basic rate parameters:
			if (nSampRate >= 1)
			{
				for (i2 in 1:nSampRate)
				{
					sampleVaryingParameters(change, bgain, onlyRates=TRUE)
				}
			}
			}
			if (any(z$set2) & (!frequentist))
			{
				for (i2 in 1:nSampCons)
				{
					sampleConstantParameters(change)
				}
			}
			else
			{
				z$acceptEta <<- 2*targetMHProb
			# 0.5 is twice the target probability of acceptance in improveMH,
			# and signifies in improveMH that nothing is to be done
			# with respect to z$scaleFactor0 and z$scaleFactorEta.
			# The "twice" is because
			# the two methods for proposing eta are used in alternation,
			# so only half the number of runs is used.
			}
			# update grand means muTemp and covariance SigmaTemp
			# do we need to check that there are more than one group?
			if ( z$nGroup > 1)
			{
				if (!frequentist)
				{
				# sample the global parameters:
				sampleMuSigma()
				}
				zm$posteriorMu[i, ] <<- z$muTemp
				zm$posteriorEta[i,] <<- z$thetaMat[1,z$set2]
				zm$PosteriorSigma[i, ,] <<- z$SigmaTemp
			}
			zm$accepts[, i] <<- mcmcc.accept
			acceptsEta[i] <- z$acceptEta
			zm$MHaccepts[i, , , ] <<-
				t(do.call(cbind,
						  tapply(ans$accepts, factor(z$callGrid[, 1]),
								 function(x)Reduce("+", x))))
			zm$MHrejects[i, , , ] <<-
				t(do.call(cbind, tapply(ans$rejects, factor(z$callGrid[, 1]),
										function(x)Reduce("+", x))))
			zm$MHaborts[i, , , ] <<- t(do.call(cbind,
										tapply(ans$aborts,
												factor(z$callGrid[, 1]),
												function(x)Reduce("+", x))))
		}
		etaAcceptances0 <- sum(acceptsEta[2*(1:(nrunMH %/% 2))])
		etaAcceptances1 <- sum(acceptsEta[2*(1:(nrunMH %/% 2))-1])
		zm$BayesAcceptances <<- c(rowSums(zm$accepts),
					etaAcceptances0, etaAcceptances1)
		zm$ans.last <<- ans
#if(change)
#{
#cat('\n mcmcc \n', print(z$thetaMat), '\n')
#browser()
#}
		#zm
	} # end MCMCcycle

	##@thetaMati internal sienaBayes; coordinate for group i
	thetaMati <- function(i)
	{
			par.thisgroup  <- union(z$ratePositions[[i]],
									which(z$varyingObjectiveParameters))
			par.thisgroup <- sort(union(par.thisgroup, which(z$set2)))
			z$thetaMat[i,par.thisgroup]
	}

	##@sampleVaryingParameters internal sienaBayes; propose new parameters
	## and accept them or not
	## Replace accepted new parameters only if change.
	## If z$incidentalBasicRates, for the basic rate parameters
	## not a Bayesian MH step but a Robbins Monro step is made.
	sampleVaryingParameters <- function(change=TRUE, bgain, onlyRates=FALSE)
	{
		## require(MASS)
		if (z$incidentalBasicRates)
		{
			if ((change)&&(bgain > 0))
			{
		# Robbins-Monro step for basic rate parameters
		# Basic rate parameters truncated to minimum value of 0.1.
		# Smaller gain parameter is used because the groupwise rate parameters
		# are less stable; and these steps are done much more
		# frequently than those for the global parameters.
				scores <- getProbabilitiesFromC(z, getScores=TRUE)[[2]]
				z$thetaMat[,z$basicRate] <<-
					pmax(z$thetaMat[,z$basicRate] +
					bgain*z$factorsBasicRate*t(scores[z$basicRate,]), 0.1)
		# Drop the structurally zero elements for storage purposes.
				z$thescores <<- sapply(1:z$nGroup,
					function(i){scores[z$ratePositions[[i]], i]})
			}
		# Basic rate parameters must be masked:
# when this will be revived, it should be thoroughly checked;
# the two calls of dmvnorm below used to be
#		dmvnorm(tmp[use],  mean=z$muTemp[restrictProposal],
#					sigma=z$SigmaTemp[restrictProposal,restrictProposal]) #density
			restrictProposal <- z$objectiveInVarying
		}
		else # !z$incidentalBasicRates
		{
			if (onlyRates)
			{
				# sample only rate parameters
				restrictProposal <- !z$objectiveInVarying
			}
			else if (z$priorRatesFromData < 0)
			{
				restrictProposal <- rep(TRUE, sum(z$objectiveInVarying))
			}
			else
			{
				restrictProposal <- rep(TRUE, z$p1) # no restriction
			}
		}
	 	# do the fandango
		## get a multivariate normal with covariance matrix proposalCov
		## multiplied by a scale factor which varies between groups
		#  propose the changes in the varying parameters
			thetaChanges <- t(sapply(1:z$nGroup, function(i)
					{
						tmp <- z$thetaMat[i, z$set1]
						use <- (!is.na(tmp))
						if (onlyRates){use[!z$basicRate[z$set1]] <- FALSE}
		# !is.na to get rid of parameters for other groups
						tmp[use] <- mvrnorm(1, mu=rep(0, sum(use)),
					Sigma=z$scaleFactors[i] *
					z$proposalCov[[i]][restrictProposal,restrictProposal] )
						tmp
					}
				))
		thetaOld <- z$thetaMat
		if (priorRatesFromData >= 0)
		{
		thetaOld[, z$basicRate] <- trafo(thetaOld[, z$basicRate])
		}
		thetaNew <- thetaOld
		thetaNew[, z$set1] <- thetaOld[,z$set1] + thetaChanges

# The following truncation was added by Tom 28/05/14
# This does not produce a symmetric MH rule,
# but will work as long as the posterior distribution
# of the rate parameters stays well above the value 0.01.
# For usual data sets, thetaNew[, z$basicRate] < 0.01 may occur,
# although with a very low probability.
		if (onlyRates)
		{
			restrictProposal <- rep(TRUE, z$p1) # no restriction from here on
		}
		if (priorRatesFromData >= 0)
		{
		thetaNew[, z$basicRate] <- ifelse(thetaNew[, z$basicRate] < 0.01,
							thetaOld[, z$basicRate], thetaNew[, z$basicRate] )
		}
 		priorOld <- sapply(1:z$nGroup, function(i)
				{
					tmp <- thetaOld[i, z$set1]
					use <- !is.na(tmp)
					dmvnorm2(tmp[use],  mean=z$muTemp, evs=z$SigmaTemp.evs,
										sigma.inv=z$SigmaTemp.inv) #density
				}
				)
		priorNew <- sapply(1:z$nGroup, function(i)
				{
					tmp <- thetaNew[i, z$set1]
					use <- !is.na(tmp)
					dmvnorm2(tmp[use],  mean=z$muTemp, evs=z$SigmaTemp.evs,
										sigma.inv=z$SigmaTemp.inv) #density
				}
				)
		# use z still with z$thetaMat = thetaOld:
		logpOld <- getProbabilitiesFromC(z)[[1]]
		if (z$priorRatesFromData >= 0)
		{
			thetaNew[, z$basicRate] <- antitrafo(thetaNew[, z$basicRate])
		}
		z$thetaMat <<- thetaNew
		logpNew <- getProbabilitiesFromC(z)[[1]]

		if (any(is.na(logpNew))) # This was a bug, hopefully does not occur any more
		{
			save(z,file="ErrorPartialBayesResult.RData")
			cat("Current object saved as ErrorPartialBayesResult.RData \n")
			cat(print(z$thetaMat), "\n")
			cat("(is.na(logpNew)) for ", which(is.na(logpNew)), "\n")
			cat("thetaMat for NA coordinates\n")
			for (i in which(is.na(logpNew))){byM(thetaMati(i))}
		}
		proposalProbability <- priorNew - priorOld + logpNew - logpOld
		acceptProposals <- (log(runif(length(proposalProbability))) <
							proposalProbability)
#		acceptProposals[is.na(acceptProposals)] <- FALSE # very rare but it did happen
		if (z$priorRatesFromData >= 0)
		{
			thetaOld[, z$basicRate] <- antitrafo(thetaOld[, z$basicRate])
		}
		if (change)
		{
			z$thetaMat[!acceptProposals, ] <<- thetaOld[!acceptProposals, ]
			proposals.accept <- acceptProposals
		}
		else
		{# used for improveMH, so that it is better to use probabilities
		    	# for proposals.accept than actual acceptances
			z$thetaMat <<- thetaOld
			proposals.accept <- exp(pmax(pmin(proposalProbability,0), -100))
			# truncation at -100 to avoid rare cases of production of NAs
		}
		if (any(is.na(proposals.accept)))
		{
			cat("any(is.na(proposals.accept))\n")
			cat("The group/s with NA proposals.accept is/are ", which(is.na(proposals.accept)),"\n")
			if (initgainGroupwise > 0.0)
			{
				cat("\nPerhaps a case of divergence?\n")
				cat("\nIf so: perhaps use a smaller value of initgainGroupwise.\n\n")
			}
			cat("Hit return to continue....")
			browser()
		}
		proposals.accept
	} # end sampleVaryingParameters


	##@sampleConstantParameters internal sienaBayes; propose new parameters
	## eta (the non-varying i.e. constant parameters), and accept them or not
	## Replace accepted new parameters only if change.
	## (varying and nonvarying)
	sampleConstantParameters <- function(change=TRUE)
	{
	# do the chamarrita
	# get a multivariate normal with covariance matrix
	# proposalC0eta multiplied by a scale factor
		etaChange <- mvrnorm(1, mu=rep(0, sum(z$set2)),
							Sigma=z$scaleFactorEta * z$proposalC0eta)
		thetaOld <- z$thetaMat
		thetaNew <- thetaOld
		thetaNew[, z$set2] <- thetaOld[, z$set2] +
						matrix(etaChange, z$nGroup, z$p2, byrow=TRUE)
		# use z still with z$thetaMat = thetaOld:
		logpOld <- getProbabilitiesFromC(z)[[1]]
		z$thetaMat <<- thetaNew
		logpNew <- getProbabilitiesFromC(z)[[1]]
		proposalProbability <- sum(logpNew - logpOld)
# cat("eta proposal  ", proposalProbability,  ", ", logpNew, ", ", logpOld, "\n")
		constantPar.accept <- (log(runif(1)) < proposalProbability)
		if (is.na(constantPar.accept)){constantPar.accept <- FALSE}
		if (change)
		{
			z$acceptEta <<- constantPar.accept
			if (!constantPar.accept)
			{
			z$thetaMat <<- thetaOld
			}
		}
		else
		{# used for improveMH, so that it is better to use probabilities
		    	# for z$accept than actual acceptances
			z$acceptEta <<- exp(min(proposalProbability, 0))
			z$thetaMat <<- thetaOld
		}
	} # end sampleConstantParameters

	##@sampleMuSigma internal sienaBayes;
	## Gibbs algorithm to sample new Mu and Sigma parameters
	sampleMuSigma <- function()
		{
		invWish <- function(v,S){
		# Draw from the inverse Wishart distribution with df v
		# and scale matrix S
		# inefficient if drawn multiply for the same S
		protectedInverse(rWishart(1,v,protectedInverse(S))[,,1])
		}
		if (priorRatesFromData < 0)
		{
			Thetas <- t(matrix( z$thetaMat[!is.na(z$thetaMat)],
									z$nGroup, z$TruNumParsPlus))[z$varyingGeneralParametersInGroup,]
		}
		else
		{
		Thetas <- t(matrix( z$thetaMat[!is.na(z$thetaMat)],
									z$nGroup, z$TruNumPars))[z$randomParametersInGroup,]
		}
		muhat <- rowMeans( Thetas )# the average across groups
		matQ <- (z$nGroup - 1)*cov(t(Thetas))
# prior Lambda = z$priorDf*z$priorSigma
		z$SigmaTemp <<- invWish(z$priorDf + z$nGroup,
				z$priorDf*z$priorSigma + matQ +
				(z$priorKappa*z$nGroup/(z$priorKappa + z$nGroup))*
				tcrossprod( muhat - z$priorMu, muhat - z$priorMu ) )
		z$muTemp <<- t(chol( ( z$priorKappa + z$nGroup )^(-1)*z$SigmaTemp )) %*%
			rnorm( z$p1 , 0 , 1 ) +
			(z$nGroup/( z$priorKappa + z$nGroup ))*muhat +
			(z$priorKappa/( z$priorKappa + z$nGroup ))*z$priorMu
	}

	##@RobbinsMonro internal sienaBayes;
	## Robbins Monro step for Mu, eta, and Sigma
	RobbinsMonro <- function(bgain, iPhase3, change=TRUE)
	{
		if (z$nGroup <= 1)
		{
			stop("If there is only 1 group, don't use frequentist method.")
		}
		if ((z$incidentalBasicRates) | (priorRatesFromData <0))
		{
		# The RM update for the basic rate parameters is then done in
		# sampleVaryingParameters. Therefore the rate parameters are masked.
			restrictUpdate <- z$varyingGeneralParametersInGroup
			updated <- z$objectiveInVarying
		}
		else
		{
			restrictUpdate <- z$varyingParametersInGroup
			updated <- rep(TRUE, z$p1)
		}

		Thetas <- matrix( z$thetaMat[!is.na(z$thetaMat)],
						z$nGroup, z$TruNumParsPlus)[, restrictUpdate, drop=FALSE]
		muhat <- colMeans(Thetas)
		matQ <- ((z$nGroup - 1)/z$nGroup)*cov(Thetas)

		if (change)
		{
			#update muTemp
			correct <- 1
			maxdif <- max(abs(z$muTemp[updated] - muhat))
			# Truncate if necessary;
			# the threshold is pretty arbitrary of course.
			# It would be better
			# to have truncation dependent on provisional standard errors.
			# Perhaps use standard deviations in phase 1?
			if (bgain*maxdif > 0.3)
			{
				correct <- 0.3/(bgain*maxdif)
				cat("Correction for jump in mu by ",correct,"\n")
			}
			# The update needs to be applied only to the unmasked part,
			# because the masked part plays no role.
			z$muTemp[updated] <<-  z$muTemp[updated] -
		                  bgain * correct * (z$muTemp[updated] - muhat)
# or with a matrix:
#			z$muTemp[updated] - bgain *
#				as.vector(z$dmuinv %*% (z$muTemp[updated] - muhat))
#			standdev <- sqrt(diag(z$SigmaTemp[updated,updated, drop=FALSE]))
#			lsd <- length(standdev)
		}
		else
		{
			z$ThinP3MuHat[iPhase3,] <<- muhat
			z$ThinP3SigmaHat[iPhase3,,] <<- matQ
		}

		if (sum(z$set2) > 0)
		{
		 # update eta
			scores <- getProbabilitiesFromC(z, getScores=TRUE)[[2]]
			etaScores <- rowSums(scores[z$set2,,drop=FALSE])
#browser()
#cat("difference = ", max(abs(scores.ans.last + (scores))))
#cat("  getProbabilitiesFromC\n")
#browser()
#prc <- getProbabilitiesFromC(z, getScores=TRUE)
#str(prc,1)
#prc[[3]]
			if (change)
			{
				correct <- 1
				maxdif <- max(abs(z$factorsEta*etaScores))
				if (bgain*maxdif > 0.3)
				{
					correct <- 0.3/(bgain*maxdif)
					cat("Correction for jump in eta by ",correct,"\n")
				}
				z$thetaMat[,z$set2] <<- z$thetaMat[,z$set2] +
							bgain * correct * z$factorsEta*etaScores
			}
			else
			{
				z$ThinP3EtaScores[iPhase3,] <<- etaScores
			}
#cat("updateEta", bgain * correct * z$factorsEta*etaScores, "\n")
#browser()
		}

		if (change)
		{
			# update SigmaTemp
			maxdif <- max(abs(z$SigmaTemp[updated,updated] -
										matQ[updated,updated]))
			correct <- 1
			if (bgain*maxdif > 0.2)
			{
				correct <- 0.2/(bgain*maxdif)
				cat("Correction for jump in Sigma by ",correct,"\n")
			}
			z$SigmaTemp <<- correctMatrix(z$SigmaTemp -
					bgain * correct * (z$SigmaTemp - matQ), z$delta)
		}
	}# end RobbinsMonro

	correctMatrix <- function(S, delt){
	# Give S eigenvalues at least delt.
		# Adapted from function make.positive.definite in package corpcor
		# which uses a method by Higham (Linear Algebra Appl 1988)
		# but changed to make the matrix positive definite (psd)
		# instead of nonnegative definite.
		# The idea is to left-truncate all eigenvalues to delta0.
		# The construction with tol, not used now,
		# is to ensure positive definiteness given numerical inaccuracy.
		# This is applied to the correlation matrix;
		# as sqrt(diagonal elements) of SigmaTemp used for the standardization,
		# the old (before update) values standdev are used,
		# because those are known to be positive.
#			corrMatrix <- diag(1/standdev, nrow=lsd) %*%
#				S %*% diag(1/standdev, nrow=lsd)
		es <- eigen(S)
#			es <- eigen(corrMatrix)
			esv <- es$values
			# Read the following maths as if delta = .Machine$double.eps = 0.
#			if (min(esv) <= .Machine$double.eps)
		if (min(esv) < delt)
			{
				cat("correctMatrix: Eigenvalues Sigma = ", sort(esv), "\n")
#				tol <- dd*max(abs(esv))*.Machine$double.eps
#				delt <-  2*tol # factor two is just to make sure the resulting
				# matrix passes all numerical tests of positive definiteness
				# TS: I replace this delta by a larger number
				tau <- pmax(0, delt - esv)
				cat("Smallest eigenvalue of Sigma now is ",
								min(esv),"; make posdef.\n")
				dm =  es$vectors %*%
						diag(tau, nrow=length(tau)) %*% t(es$vectors)
#		diag(standdev, nrow=length(tau)) %*% (corrMatrix + dm) %*% diag(standdev, nrow=length(tau))
				S <- S + dm
				z$correctSigma <<- z$correctSigma + 1
			}
		S
		}

	##@averageTheta internal sienaBayes; algorithm to average past theta values
	averageTheta <- function()
	{
		thetaMean <- rep(NA, z$pp)
		for (group in 1:z$nGroup)
		{
			thetaMean[z$ratePositions[[group]]] <- 	colMeans(
				z$ThinParameters[, group, !z$generalParametersInGroup,
						drop=FALSE], na.rm=TRUE)
		}
		if ((priorRatesFromData <0) | incidentalBasicRates)
		{
		# then (z$set1)&(!z$basicRate) == (z$set1); just for clarity we write:
		thetaMean[(z$set1)&(!z$basicRate)] <- colMeans(
					z$ThinPosteriorMu[,, drop=FALSE], na.rm=TRUE)
		}
		else
		{
			thetaMean[(z$set1)&(!z$basicRate)] <- colMeans(
					z$ThinPosteriorMu[,z$objectiveInVarying,
							drop=FALSE], na.rm=TRUE)
		}
		thetaMean[z$set2] <- colMeans(z$ThinPosteriorEta[,, drop=FALSE], na.rm=TRUE)
		thetaMean[z$fix & (!z$basicRate)] <- z$thetaMat[1,z$fix & (!z$basicRate)]
		thetaMean
	}

	##@somePlots make some plots during operation of the function
	# dropped in version 1.1-278.
	# probably requires storeAll
	##@savePartial makes a partial intermediate save as a sienaBayesFit object
	savePartial <- function()
	{# save intermediate results
		saveTheta <- z$theta
		saveSigma <- z$Sigma
		if (frequentist)
		{
			z$theta <<- c(z$muTemp, z$thetaMat[1, z$set2])
			# but this does not duplicate the rate parameters,
			# whereas averageTheta does duplicate. to be changed.
			z$Sigma <<- z$SigmaTemp
			# these have not changed during Phase 3 (check!)
		}
		else
		{
			z$theta <<- averageTheta()
		}
		save(z,file="PartialBayesResult.RData")
		z$theta <<- saveTheta
		z$Sigma <<- saveSigma
	}

##@covtrob internal sienaBayes protected robust covariance matrix, with a small ridge
covtrob <- function(x){
	if (inherits(try(cr <- cov.trob(x)$cov, silent=TRUE), "try-error")){
		cr <- cov(x) + 0.05*diag(diag(cov(x)), nrow = dim(x)[2])
		}
	cr + 0.05*diag(diag(cr), nrow = dim(x)[2])
}
    ## ###################################### ##
    ## start of function sienaBayes() proper  ##
    ## ###################################### ##
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
	ctime1 <- proc.time()
	ctime <- proc.time()
	if (any(effects$test))
	{
		cat("Specification that some effects are tested is canceled.\n")
		effects$test[effects$test] <- FALSE
	}

	if (is.null(prevBayes))
	{
		z <- initializeBayes(data, effects, algo, nbrNodes,
						initgainGlobal=initgainGlobal,
						initgainGroupwise=initgainGroupwise, initML=initML,
						priorMu=priorMu, priorSigma=priorSigma,
						priorDf=priorDf, priorKappa=priorKappa,
						priorRatesFromData=priorRatesFromData,
						frequentist=frequentist,
						incidentalBasicRates=incidentalBasicRates,
						reductionFactor=reductionFactor, gamma=gamma,
						delta=delta, nmain=nmain, nwarm=nwarm,
						nImproveMH=nImproveMH,
						lengthPhase1=lengthPhase1, lengthPhase3=lengthPhase3,
						prevAns=prevAns, usePrevOnly=usePrevOnly,
						silentstart=silentstart, clusterType=clusterType)
		cat("Initial global model estimates\n")
		print(z$initialResults)
		flush.console()
		# Throw out the largest parts of z$f to save memory;
		# some part of z$f is needed for print.summary.
		# An alternative would be to retain just that part.
		z$f$minimalChain <- NULL
		z$f$chain <- NULL
		f <- FRANstore()
#		f[1:z$nGroup] <- NULL
		f$minimalChain <- NULL
		f$chain <- NULL
		FRANstore(f)
		createStores()
		z$sub <- 0
		zm <- list()
		zsmall <- list()
		if (nrunMHBatches >= 2)
			{
				z$nrunMHBatches <- nrunMHBatches
			}
		else
			{
				cat("nrunMHBatches increased to 2.\n")
				z$nrunMHBatches <- 2
				nrunMHBatches <- 2
			}
		z$nSampVarying <- min(nSampVarying, 1)
		if (nSampVarying < 1)
		{
			cat("nSampVarying increased to 1.\n")
		}
		z$nSampRates <- nSampRates

		if ((nSampRates > 0) && ((incidentalBasicRates) | (priorRatesFromData < 0)))
		{
			stop("If incidentalBasicRates or priorRatesFromData < 0, nSampRates should be 0.")
		}
		z$nSampConst <- min(nSampConst, 1)
		if (nSampConst < 1)
		{
			cat("nSampConst increased to 1.\n")
		}
		class(z) <- "sienaBayesFit"
		z$correctSigma <- 0

		# Determine multiplicative constants for proposal distribution
		ctime1 <- proc.time()
		cat(round((ctime1-ctime)[3]),' seconds elapsed\n')
		improveMH(totruns=nImproveMH, target=targetMHProb)
		ctime2<- proc.time()
		cat('improveMH', round((ctime2-ctime1)[3]),' seconds elapsed for improveMH.\n')
		flush.console()
	}
	else #  (!is.null(prevBayes))
	{
		if (inherits(prevBayes, "sienaBayesFit"))
		{
				cat("Continuation from earlier sienaBayes result",
					deparse(substitute(prevBayes)),".\n\n")
		}
		else
		{
			cat(deparse(substitute(prevBayes)),
						"is not a sienaBayesFit object.\n")
			stop("Do not use this object for prevBayes.")
		}
		z <- prevBayes

		if ((nrunMHBatches != z$nrunMHBatches) & (nrunMHBatches >= 2))
		{
			z$nrunMHBatches <- nrunMHBatches
			cat("nrunMHBatches changed to ", nrunMHBatches, ".\n")
		}
		if ((nSampVarying != z$nSampVarying) & (nSampVarying >= 1))
		{
			z$nSampVarying <- nSampVarying
			cat("nSampVarying changed to ", nSampVarying, ".\n")
		}
		if ((nSampRates != z$nSampRates) & (nSampRates >= 0))
		{
			z$nSampRates <- nSampRates
			cat("nSampRates changed to ", nSampRates, ".\n")
		}
		if ((nSampConst != z$nSampConst) & (nSampConst >= 1))
		{
			z$nSampConst <- nSampConst
			cat("nSampConst changed to ", nSampConst, ".\n")
		}
		flush.console()
		correctMatrix

		esv <- eigen(z$priorSigma)$values
		if (min(esv) <= 0)
		{
			cat('The matrix priorSigma is not positive definite.\n')
			dS <- diag(z$priorSigma)
			if (min(dS) <= 0)
			{
				cat('Its diagonal elements are not all positive.\n')
				stop('Incorrect priorSigma for prevBAyes object.')
			}
			cat('Off-diagonal values are set to 0.\n')
			pS <- matrix(0, dim(z$priorSigma)[1], dim(z$priorSigma)[1])
			diag(pS) <- dS
			z$priorSigma <- pS
		}
		z$correctSigma <- 0
	}

# The division into phases 1-3 is relevant only for frequentist estimation.

	if (is.null(prevBayes))
	{
		endPhase1 <- nwarm + z$lengthPhase1
		endPhase2 <- nwarm + nmain - z$lengthPhase3
		z$frequentist <- frequentist

		# Warming phase
		bgain <- initfgain
		for (ii in 1:nwarm)
		{
			MCMCcycle(nrunMH=nrunMHBatches, nSampVar=nSampVarying,
								nSampCons=nSampConst, nSampRate=nSampRates,
								bgain=bgain)
			z$iimain <- ii
			storeData()
			cat('Warming step',ii,'(',nwarm,')\n')
			cat("Accepts ",sum(zm$BayesAcceptances),"/",
				z$nGroup*nrunMHBatches,"\n")
			flush.console()
		}
		print('end of warming')
		ctime3<- proc.time()
		cat('warming took', round((ctime3-ctime2)[3]),'seconds elapsed for warming.\n')
		flush.console()

		# Improve parameter value by averaging
	# Average the groupwise parameters over the last half of Phase 2
#cat('aaa\n')
#browser()
		halfPhaseWarm <- floor(nwarm/2 + 1)
		for (group in 1:z$nGroup)
		{
			if (priorRatesFromData >= 0)
			{
			z$thetaMat[group, z$ratePositions[[group]]] <-
				colMeans(z$ThinParameters[halfPhaseWarm:nwarm, group,
				!z$generalParametersInGroup, drop=FALSE], dims=1)
			}
			z$thetaMat[group, !z$basicRate] <- colMeans(
				z$ThinParameters[halfPhaseWarm:nwarm, group,
				z$generalParametersInGroup, drop=FALSE], dims=1)
		}
		# Average the past groupwise parameters over runs and groups
		# Average the covariance matrices over the runs
		if (priorRatesFromData >= 0)
		{
			z$muTemp <-
				colMeans(z$ThinParameters[halfPhaseWarm:nwarm, ,
						z$randomParametersInGroup, drop=FALSE], dims=2)
			tmp <- lapply(halfPhaseWarm:nwarm, function(i, x)
				{cov(x[i,,z$randomParametersInGroup])},
										x=z$ThinParameters)
		}
		else
		{# or should this also be randomParametersInGroup?
		z$muTemp <-
			colMeans(z$ThinParameters[halfPhaseWarm:nwarm, ,
						z$varyingGeneralParametersInGroup, drop=FALSE], dims=2)
		tmp <- lapply(halfPhaseWarm:nwarm, function(i, x)
				{cov(x[i,,z$varyingGeneralParametersInGroup])},
										x=z$ThinParameters)
		}
		z$SigmaTemp <- Reduce('+', tmp) / (nwarm - halfPhaseWarm + 1)
		cat('Parameter values after warming up\n')
		for (i in (1:z$nGroup))
		{
			cat(i,". ")
			byM(thetaMati(i))
		}
		cat('\n')
		ctime1 <- proc.time()
		cat('Second ')
		improveMH(totruns=nImproveMH, target=targetMHProb)
		ctime4<- proc.time()
		cat('Second improveMH', round((ctime4-ctime1)[3]),'seconds elapsed for second improveMH.\n')
		flush.console()
	}
	else # (!is.null(prevBayes))
	{ # what remains to be done of the initialization
		nwarm <- 0
		z$nwarm <- 0
		z$nmain <- nmain
		createStores()
		z$Phase <- 1
		if (frequentist)
		{
			lengthPhase1 <- max(lengthPhase1, 2)
			lengthPhase3 <- max(lengthPhase3, 2)
		}
		else
		{
			lengthPhase1 <- round(nmain/5)
			lengthPhase3 <- round(nmain/5)
		}
		z$lengthPhase1 <- lengthPhase1
		z$lengthPhase3 <- lengthPhase3
		endPhase1 <- z$lengthPhase1
		endPhase2 <- nmain - z$lengthPhase3

		algo$maxlike <- TRUE
		algo$FRANname <- "maxlikec"

		z$print <- FALSE
		z$int <- 1
		z$int2 <- nbrNodes
		algo$cconditional <-  FALSE
		if (!is.null(algo$randomSeed))
		{
			set.seed(algo$randomSeed)
		}
		else
		{
			if (exists(".Random.seed"))
			{
				rm(.Random.seed, pos=1)
			}
			newseed <- trunc(runif(1) * 1000000)
			set.seed(newseed)  ## get R to create a random number seed for me.
		}
		ctime4<- proc.time()
		prevThetaMat <- z$thetaMat
		z$FRAN <- getFromNamespace(algo$FRANname, pkgname)
		z <- initializeFRAN(z, algo, data=data, effects=effects,
                prevAns=prevAns, initC=FALSE, onlyLoglik=TRUE)
		cat('back in sienaBayes\n')
		z$thetaMat <- prevThetaMat

		if (nbrNodes > 1 && z$observations > 1)
		{
			## require(parallel)
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
			clusterCall(z$cl, initializeFRAN, z, algo,
                    initC = TRUE, profileData=FALSE, returnDeps=FALSE)
			clusterSetRNGStream(z$cl, iseed = as.integer(runif(1,
								max=.Machine$integer.max)))
		}

		if (newProposalFromPrev)
		{
		# effects0 defined as in initializeBayes,
		# used for determining which groupwise parameters are needed.
		# This is the same for all groups, therefore group 1 is used.
			effects0 <- getEffects(data[[1]])
		# projectEffects copies the model specification,
		# except for the basic rate effects for other groups,
		# and it fixes non-random effects;
		# it also copies estimated values, which is superfluous here but does no harm.
			effects0 <- projectEffects(effects0,effects,z$initialResults,1)
			if (is.null(prevBayes$nwarm))
			{
				nfirst <- 1
			}
			else
			{
				nfirst <- prevBayes$nwarm + 1
			}
			ntot <- max(which(!is.na(prevBayes$ThinParameters[,1,1])))
			if (z$p2 > 0)
			{
				z$proposalC0eta <-
					covtrob(prevBayes$ThinPosteriorEta[nfirst:ntot,, drop=FALSE])
			}
			z$proposalCov <- lapply(1:z$nGroup, function(i){
				covtrob(apply(prevBayes$ThinParameters[nfirst:ntot,i,
														!effects0$fix[effects0$include], drop=FALSE],
								c(1,3), function(x)x))})
# the inner apply is used here to get rid of the middle dimension i;
# this dimension is retained because of the drop=FALSE,
# and this is necessary because the third coordinate might have just one value
# (which will hardly occur in practice...)
			scf <- 2.38/sqrt(z$TruNumPars)
# theoretically optimal value according to Roberts & Rosenthal, 2001
			z$scaleFactors <- rep(scf, z$nGroup)
			if (z$p2 > 0)
			{
				z$scaleFactorEta <- 2.38/sqrt(z$p2)
			}
			else
			{
				z$scaleFactorEta <- 1
			}
			z$scaleFactor0 <- z$scaleFactorEta
			ctime01 <- proc.time()
			cat('For prevBayes ')
			improveMH(totruns=nImproveMH, target=targetMHProb)
			ctime04<- proc.time()
			cat('From prevBayes improveMH', round((ctime04-ctime01)[3]),'seconds elapsed.\n')
			flush.console()
		}
	}

	z$OK <- FALSE
	if (saveFreq >= 2)
	{
		savePartial()
	}

	# Main iterations
	# Phase 1
	bgain <- initfgain
    for (ii in (z$nwarm+1):endPhase1)
    {
		MCMCcycle(nrunMH=nrunMHBatches, nSampVar=nSampVarying,
								nSampCons=nSampConst,
								nSampRate=nSampRates, bgain=bgain)
		z$iimain <- ii
		storeData()
		cat('main', ii, '(', nwarm+nmain, ')')
		if (frequentist) cat(', Phase 1')

		cat("\nMu = ", round(z$muTemp,3), "\nEta = ",
					round(z$ThinPosteriorEta[ii,],3), "\nSigma = \n")
		print(round(z$SigmaTemp,4))
		cat("\n")
		cat('main', ii, '(', nwarm+nmain, ')',
			" Accepts ",sum(zm$BayesAcceptances),"/",
				z$nGroup*nrunMHBatches,"\n")
#	cat("thetaMat = \n")
#	print(round(z$thetaMat,2))
#		cat("\n")
		flush.console()
#		scores.ans.last <- t(zm$ans.last$fra) # the scores a z$pp by z$nGroups matrix
#		dfra.ans.last <- zm$ans.last$dff # the Hessian a sparse z$pp by z$pp matrix
		if (frequentist)
		{
			RobbinsMonro(bgain)
		}
		# save intermediate results
		if (saveFreq >= 2)
		{
			if (ii %% saveFreq == 0)
			{
				savePartial()
			}
		}
    }

#	cat("thetaMat = \n")
#	print(round(z$thetaMat,2))
	ctime5<- proc.time()
	if (frequentist) {cat('For phase 1 ', round((ctime5-ctime4)[3]),' seconds elapsed.\n')}
	flush.console()

#	cat('Third ')
#	improveMH()
	flush.console()

	# Phase 2
    for (ii in (endPhase1+1):endPhase2)
    {
		bgain <- initfgain*exp(-gamma*log(ii-endPhase1+1))
		MCMCcycle(nrunMH=nrunMHBatches, nSampVar=nSampVarying,
								nSampCons=nSampConst,
								nSampRate=nSampRates, bgain=bgain)
		z$iimain <- ii
		storeData()
		cat('main', ii, '(', nwarm+nmain, ')')
		if (frequentist) {cat(', Phase 2')}
		cat("\nMu = ", round(z$muTemp,3), "\nEta = ",
				round(z$ThinPosteriorEta[ii,],3), "\nSigma = \n")
		print(round(z$SigmaTemp,3))
		cat("\n")
		cat('main', ii, '(', nwarm+nmain, ')',
			" Accepts ",sum(zm$BayesAcceptances),"/",
				z$nGroup*nrunMHBatches,"\n")
#		print(round(z$SigmaTemp,4))
#		cat("thetaMat = \n")
#		print(round(z$thetaMat,2))
#		cat("\n")
		flush.console()
		if (frequentist)
		{
			RobbinsMonro(bgain)
		}
		# save intermediate results
		if (saveFreq >= 2)
		{
			if (ii %% saveFreq == 0)
			{
				savePartial()
			}
		}
    }

	# Process results Phase 2.
	prop <- 0.2
	partPhase2 <- floor((prop*endPhase1 + (1-prop)*endPhase2) + 1)
	if (frequentist)
	{
	# Average muTemp, eta (in thetaMat) and SigmaTemp
	# over the last half of Phase 2
		z$muTemp <- colMeans(z$ThinPosteriorMu[partPhase2:endPhase2, ])
		z$thetaMat[1:z$nGroup,z$set2] <-
				colMeans(z$ThinPosteriorEta[partPhase2:endPhase2,,drop=FALSE])
		z$SigmaTemp <- colMeans(z$ThinPosteriorSigma[partPhase2:endPhase2,,],
												dims=1)
	}
	if (z$incidentalBasicRates)
	{
	# Average the rate parameters over the last half of Phase 2
		for (group in 1:z$nGroup)
		{
			z$thetaMat[group, z$ratePositions[[group]]] <-
				colMeans(z$ThinParameters[partPhase2:endPhase2, group,
					!z$generalParametersInGroup, drop=FALSE], dims=1)
		}
	}

	ctime6<- proc.time()
	if (frequentist){cat('Phase 2 took', round((ctime6-ctime5)[3]),'seconds elapsed.\n')}
	cat("\nMu = ", round(z$muTemp,3), "\nEta = ",
				round(z$ThinPosteriorEta[ii,],3), "\nSigma = \n")
	print(round(z$SigmaTemp,3))
	cat("\n")
	flush.console()
	savePartial()

# Phase 3
	if (frequentist){cat('Phase 3\n')}
    for (ii in (endPhase2+1):(nwarm+nmain))
    {
		MCMCcycle(nrunMH=nrunMHBatches, nSampVar=nSampVarying,
								nSampCons=nSampConst,
								nSampRate=nSampRates, bgain=0.000001)
		z$iimain <- ii
		storeData()
		cat('main', ii, '(', nwarm+nmain, ')')
		if (frequentist) cat(', Phase 3')
		cat("\nMu = ", round(z$muTemp,3), "\nEta = ",
				round(z$ThinPosteriorEta[ii,],3), "\nSigma = \n")
		print(round(z$SigmaTemp,3))
		cat("\n")
#		flush.console()
		if (frequentist)
		{
			RobbinsMonro(0.0, ii-endPhase2, change=FALSE)
		}
		else
		{
			cat('main', ii, '(', nwarm+nmain, ')',
				" Accepts ",sum(zm$BayesAcceptances),"/",
				z$nGroup*nrunMHBatches,"\n")
		}
	# save intermediate results
		if (saveFreq >= 2)
		{
			if (ii %% saveFreq == 0)
			{
				savePartial()
			}
		}
    }

# Process results of Phase 3
	ctime7<- proc.time()
	cat('Total time elapsed ',round((ctime7-ctime)[3]),'seconds.\n')
	if (frequentist)
	{
		z$theta <- c(z$muTemp, z$thetaMat[1, z$set2])
		# but this does not duplicate the rate parameters,
		# whereas averageTheta does duplicate. to be changed.
		z$Sigma <- z$SigmaTemp
		# these have not changed during Phase 3 (check!)
	}
	else
	{
		z$theta <- averageTheta()
	}
	z$frequentist <- frequentist
    z$FRAN <- NULL
	if (z$correctSigma > 0)
	{
		cat('Warning: corrections of the covariance matrix to keep eigenvalues larger than ', delta,
		 ' were necessary in ', z$correctSigma, 'cases.\n')
	}
    class(z) <- "sienaBayesFit"

	if (nbrNodes > 1 && z$observations > 1)
	{
		stopCluster(z$cl)
	}
# 	if (!storeAll)
# 	{
# 		z$acceptances <- NULL
# 		z$candidates <- NULL
# 	}
	z$OK <- TRUE
    z
} # end sienaBayes

##@initializeBayes algorithms do set up for Bayesian model
initializeBayes <- function(data, effects, algo, nbrNodes,
						initgainGlobal, initgainGroupwise, initML,
						priorMu, priorSigma, priorDf, priorKappa,
						priorRatesFromData,
						frequentist, incidentalBasicRates,
						reductionFactor, gamma, delta,
						nmain, nwarm, nImproveMH,
						lengthPhase1, lengthPhase3,
						prevAns, usePrevOnly,
						silentstart, clusterType=c("PSOCK", "FORK"))
{
	##@precision internal initializeBayes invert z$covtheta
	## avoiding some inversion problems
	## for MoM estimates only
	precision <- function(z)
	{
		## require(MASS)
# Perhaps think of	t(z$dfrac) %*% ginv(z$msfc) %*% z$dfrac
# But the ...c versions are corrected and have identity rows/columns
# for fixed parameters. This correction is not to be employed here.
		t(z$dfra) %*% ginv(z$msf) %*% z$dfra
	}#end precision in initializeBayes

	## start initializeBayes
    ## initialize
	if (!inherits(data,"siena")){stop("The data is not a valid siena object.")}
	allDistances <- lapply(1:length(data),
		function(i){lapply(data[[i]]$depvars, function(y){attr(y,'distance')})})
	matDistances <- matrix(unlist(allDistances), nrow=length(data),
									ncol=length(unlist(allDistances[[1]])), byrow=TRUE)
	colnames(matDistances) <- names(unlist(allDistances[[1]]))

	if (min(matDistances) <= 0)
	{
		cat('In some groups, there is no change for a dependent variable.\n')
		cat('This is not allowed.\n')
		dist.check <- 3
		wh <- which((matDistances <= dist.check), arr.ind=TRUE, useNames=FALSE)
		cat('The following shows the dependent variables and periods\n')
		cat('for which there are less than',dist.check+1, 'changes: \n')
		prmatrix(cbind(wh[,1], colnames(matDistances)[wh[,2]], matDistances[wh]),
						rowlab=rep('',dim(wh)[1]),
						collab=c('group', 'depvar/period', '# changes'), quote=FALSE)
		cat('In any case, there should be no groups with 0 changes.\n')
		stop('No change in some variable: Adapt the data set.')
	}

	imp.change <- lapply(data, checkImpossibleChanges)
	if (sum(unlist(imp.change)) > 0)
	{
		cat('For some groups, there are impossible changes;\n')
		cat('This occurs for groups\n', which(sapply(imp.change,
						function(ic){(sum(ic) > 0)})),'\n')
		cat('This is not permitted for likelihood-based estimation')
		stop('Impossible changes in some variable: Adapt the data set.')
	}
	
	if (!is.null(algo$MaxDegree))
	{
		cat('The MaxDegree option is not implemented for likelihood-based estimation')
		stop('Do not use MaxDegree.')
	}

    Report(openfiles=TRUE, type="n") #initialize with no file
    z  <-  NULL

    z$Phase <- 1
    z$Deriv <- FALSE
    z$FinDiff.method <- FALSE
    z$maxlike <- TRUE
    algo$maxlike <- TRUE
    algo$FRANname <- "maxlikec"
    z$print <- FALSE
    z$int <- 1
    z$int2 <- nbrNodes
	z$mult <- algo$mult
    algo$cconditional <-  FALSE
    if (!is.null(algo$randomSeed))
    {
        set.seed(algo$randomSeed)
		##seed <- algo$randomSeed
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
	# ensure that all basic rates are varying between groups
	effects$randomEffects[effects$basicRate] <- rep(TRUE,sum(effects$basicRate))
	# note: this change is made only within function initializeBayes

	## Determine starting values and covariance matrices
	## for proposal distributions.
	# Starting values are obtained from MoM estimates.
	# First under the assumption that all parameters are the same.
	# Number of groups:
	if (inherits(data,"sienaGroup"))
	{
		ngroups <- length(data)
	}
	else
	{
		ngroups <- 1
	}
	if (ngroups >= 2)
	{
	# Check that the groups do not differ with respect to uponly or downonly
	# for some of the dependent variables, without allowOnly being set to FALSE.
	# This would lead to different numbers of effects between groups.
		for (j in 1:length(data[[1]]$depvars))
		{
			diverse <- ((any(sapply(1:ngroups,
						function(i)attr(data[[i]]$depvars[[j]], "uponly"))))&&
					(!all(sapply(1:ngroups,
						function(i){attr(data[[i]]$depvars[[j]], "uponly")}))))
			diverse <- diverse ||
					((any(sapply(1:ngroups,
						function(i)attr(data[[i]]$depvars[[j]], "downonly"))))&&
					(!all(sapply(1:ngroups,
					 function(i){attr(data[[i]]$depvars[[j]], "downonly")}))))
			if (diverse && any(sapply(1:ngroups,
					function(i){attr(data[[i]]$depvars[[j]], "allowOnly")})))
			{
			cat("\nThere is heterogeneity between the groups with respect to")
			cat("\ngoing up only, or down only, for dependent variable ",
					attr(data[[1]]$depvars[[j]], "name"),".\n")
			cat("Define this variable in sienaDependent() for all groups",
				"with allowOnly=FALSE.\n")
			stop("allowOnly=FALSE should be set in sienaDependent().")
			}
		}
	}

	if (incidentalBasicRates)
	{
		priorRatesFromData <- 0
	}

	# Short specification for initial multigroup model.
	if (initgainGlobal <= 0)
	{
		startupModel <- sienaAlgorithmCreate(n3=500, nsub=0, cond=FALSE,
						lessMem=TRUE,
						diagonalize=0.2, doubleAveraging=0,
						modelType=algo$modelType, behModelType=algo$behModelType,
						MaxDegree=algo$MaxDegree, Offset=algo$universalOffset,
						seed=algo$randomSeed, projname="startup_sienaBayes")
	}
	else
	{
		nsub <- 3
		if (!is.null(prevAns))
		{
			if ((usePrevOnly) &
				(dim(prevAns$requestedEffects)[1] == sum(effects$include)) &
				!prevAns$x$cconditional)
			{
				nsub <- 0
			}
		}
		startupModel <- sienaAlgorithmCreate(n3=500, nsub=nsub, cond=FALSE,
						firstg=initgainGlobal, lessMem=TRUE,
						modelType=algo$modelType, behModelType=algo$behModelType,
						MaxDegree=algo$MaxDegree, Offset=algo$universalOffset,
						seed=algo$randomSeed, projname="startup_sienaBayes")
	}
	cat("Estimate initial global parameters\n")
	if (is.null(prevAns))
	{
		startupGlobal <- siena07(startupModel, data=data, effects=effects,
								batch=TRUE, silent=silentstart,
								useCluster=(nbrNodes >= 2), nbrNodes=nbrNodes,
								clusterType=clusterType)
	}
	else
	{
		startupGlobal <- siena07(startupModel, data=data, effects=effects,
						batch=TRUE, silent=silentstart,
						useCluster=(nbrNodes >= 2), nbrNodes=nbrNodes,
						clusterType=clusterType, prevAns=prevAns)
	}
	cat("Initial global estimates\n")
	z$initialResults <- print(startupGlobal)
	z$initialResults$sf <- NULL
	z$initialResults$scores <- NULL
	cat('\n')
	effects <- updateTheta(effects, startupGlobal)
	z$initialGlobalEstimates <- effects[effects$include,]
	cat("\nmaximum initial global estimate is ",
		max(abs(startupGlobal$theta), na.rm=TRUE), "\n")
	if (max(abs(startupGlobal$theta), na.rm=TRUE) > 400)
	{
		cat("\nThe initial global estimates are\n")
		print(startupGlobal$theta)
		cat("\nThe largest absolute value is ",
				max(abs(startupGlobal$theta), na.rm=TRUE),
				"\nwhich is too large for good operation of this function.\n",
				"\nAdvice: use a smaller value of initgainGlobal;\n",
				"or, if possible, first estimate a multi-group project using siena07\n",
				"and use this as prevAns in sienaBayes with initgainGlobal=0.\n\n")
		save(startupGlobal, file="startupGlobal.RData")
		cat("startupGlobal saved.\n")
		stop("Divergent initial estimate.")
	}
# some checks
	if (with(startupGlobal$requestedEffects, any(fix&type=='rate')))
	{
		if (with(startupGlobal$requestedEffects[startupGlobal$requestedEffects$type=='rate',], !all(fix)))
		{
			print(startupGlobal$requestedEffects)
			stop('Either all, or no rate effects should be fixed.')
		}
		if (incidentalBasicRates)
		{
			stop('Fixed rate effects cannot be combined with incidentalBasicRates.')
		}
		if (priorRatesFromData >= 0)
		{
			warning('Rate parameters are fixed; priorRatesFromData changed to -1.')
			priorRatesFromData <- -1
		}
	}
	else
	{
		if (priorRatesFromData < 0)
		{
			warning('Rate parameters are not fixed, but priorRatesFromData is less than 0???')
		}
	}
	if (with(startupGlobal$requestedEffects, any(fix&randomEffects&(!basicRate))))
	{
		print(startupGlobal$requestedEffects, includeRandoms=TRUE)
		stop('Fixed effects should not be random.')
	}

	z$initialSE <- sqrt(diag(startupGlobal$covtheta))

	# Now per group.
	# Largest acceptable variance for proposal distribution rate parameters
	# on the trafo scale; note that this is before scaling by scaleFactors
	rateVarScale <- 0.25
	# Note that groupwise covtheta might sometimes be singular small groups.
	# This problem is circumvented by taking a weighted average.
	# weight for groupwise precision, relative to global precision:
	w <- 0.8

	# number of parameters per group:
	npars <- sum(startupGlobal$requestedEffects$group==1)
	# requestedEffects must be used, because startupGlobal$effects
	# includes also all main effects corresponding to interactions,
	# which is not necessarily the case for startupGlobal$requestedEffects.
	# number of Dependent variables:
#	nDepvar <- length(unique(startupGlobal$requestedEffects$name))
	if (ngroups >= 2)
	{
	#	if ((startupGlobal$pp - npars)%%(ngroups-1) != 0)
		if (var(sapply(data,function(x){x$observations})) > 1e-6)
		{
			stop("Not all groups have the same number of periods.")
		}
	}
	nBasicRatesPerGroup <-
		sum(startupGlobal$requestedEffects$basicRate) / ngroups
	if (abs(nBasicRatesPerGroup - round(nBasicRatesPerGroup)) > 1e-6)
	{
		stop("Effects object and data object do not correspond.")
	}
	rateParameters <- matrix(NA, nBasicRatesPerGroup, ngroups)
	nActors <- rep(0, ngroups)
	initialEstimates <- matrix(0, ngroups, startupGlobal$pp)
	factorsBasicRate <- matrix(0, ngroups, startupGlobal$pp)
	# Define the partition for the varying non-fixed (set1)
	# and non-varying non-fixed (set2) effects;
	# effects that are not estimated (fix) are
	# excluded from both sets.
	z$set1 <- rep(FALSE, startupGlobal$pp)
	z$set1[startupGlobal$requestedEffects$basicRate] <- TRUE
	z$set1[startupGlobal$requestedEffects$randomEffects &
				(!startupGlobal$requestedEffects$fix)] <- TRUE
	z$set2 <- !z$set1
	z$set2[startupGlobal$requestedEffects$fix] <- FALSE
	if ((priorRatesFromData < 0) || incidentalBasicRates)
	{
		z$set1[startupGlobal$requestedEffects$basicRate] <- FALSE
	}

	# Number of non-varying non-fixed parameters
	# among the z$truNumPars true parameters:
	z$p2 <- sum(z$set2)
	# Number of varying parameters among the z$truNumPars true parameters:
	if (incidentalBasicRates)
	{
		z$p1 <- npars - z$p2 - with(startupGlobal$requestedEffects,
									sum(((type == 'rate') | fix) & (group==1)))
	}
	else
	{
		z$p1 <- npars - z$p2 - with(startupGlobal$requestedEffects, sum(fix & (group==1)))
	}
	# Note that z$p1 = sum(z$set1) + (z$nGroup-1)*(number of rate parameters per group)
	# unless incidentalBasicRates

	# Some checks of input parameters:
	if (!(is.null(priorMu)))
	{
		if (length(priorMu) != z$p1)
		{
			if (incidentalBasicRates)
			{
				cat('For the incidentalBasicRates option, \n',
				    'rate parameters should not be included in the prior.\n')
			}
			stop(paste("priorMu must have length ",z$p1))
		}
	}
	if (!is.null(priorSigma))
	{
	if (!((class(priorSigma)=="matrix") & all(dim(priorSigma)==c(z$p1,z$p1))))
		{
			stop(paste("priorSigma is not a matrix of dimensions",z$p1))
		}
	}

	# compute covariance matrices for random walk proposals for all groups:
	# proposalC0 for all effects;
	# proposalC0eta for eta, i.e., non-varying effects.
	# proposalCov for theta^{(1)}_j, i.e., groupwise effects.
	z$proposalCov <- list()
	prec <- precision(startupGlobal)
	rates <- startupGlobal$requestedEffects$basicRate
	precRate <- sum(diag(prec[rates,rates]))
	prec[rates,rates] <- 0
	diag(prec)[rates] <- precRate
	prec.PA <- prec # PA = perhaps amended
	if (inherits(try(z$proposalC0 <- chol2inv(chol(prec)),
					silent=TRUE), "try-error"))
	{
		cat("precision(startupGlobal) non-invertible.\n")
		cat("This probably means the model is not identifiable.\n")
		cat("A ridge is added to this precision matrix.\n")
#		stop("Non-invertible precision(startupGlobal).")
		prec.amended <- prec + diag(0.5*diag(prec), nrow=dim(prec)[1])
		prec.PA <- prec.amended
		if (inherits(try(z$proposalC0 <- chol2inv(chol(prec.amended)),
					silent=TRUE), "try-error"))
		{
			cat("Even the amended precision(startupGlobal) is non-invertible.\n")
			print(prec.amended)
			cat("\n Please report this to Tom.\n")
			stop("Non-invertible amended precision(startupGlobal).")
		}
	}
	z$proposalC0eta <- z$proposalC0[z$set2, z$set2, drop=FALSE]
	z$forEtaVersion0 <-
		(z$set1|z$set2)&(!(rates|startupGlobal$requestedEffects$fix))
	z$proposalC0 <- z$proposalC0[z$forEtaVersion0, z$forEtaVersion0, drop=FALSE]

	for (i in 1:ngroups)
	{
		cat(paste("Group",i,"\n"))
		if (initgainGroupwise <= 0)
		{
			startup1Model <- sienaAlgorithmCreate(n3=500, nsub=0,  lessMem=TRUE,
					cond=FALSE,
					projname=paste("project",formatC(i,width=2,flag="0"),sep=""),
					seed=algo$randomSeed,
					modelType=algo$modelType, behModelType=algo$behModelType,
					MaxDegree=algo$MaxDegree, Offset=algo$universalOffset)
		}
		else
		{
			startup1Model <- sienaAlgorithmCreate(n3=10, nsub=1, lessMem=TRUE,
					firstg=initgainGroupwise, cond=FALSE, maxlike=initML,
					projname=paste("project",formatC(i,width=2,flag="0"),sep=""),
					seed=algo$randomSeed,
					modelType=algo$modelType, behModelType=algo$behModelType,
					MaxDegree=algo$MaxDegree, Offset=algo$universalOffset)
		}
		if (initML)
		{
			startup1Model3 <- sienaAlgorithmCreate(n3=500, nsub=0, lessMem=TRUE,
					firstg=initgainGroupwise, cond=FALSE, maxlike=TRUE,
					projname=paste("project",formatC(i,width=2,flag="0"),sep=""),
					modelType=algo$modelType, behModelType=algo$behModelType,
					MaxDegree=algo$MaxDegree, Offset=algo$universalOffset,
					localML=algo$localML,
					seed=algo$randomSeed)
		}
		# Give the algorithm object two extra components that can be used
		# internally in siena07 in CalculateDerivative in phase1.r
		# to prevent a prolonged phase 1 in case of
		# non-sensitivity for certain parameters.
		startup1Model$fromBayes <- TRUE
		if (sum(!effects[effects$include,"basicRate"]) >= 2)
		{
			startup1Model$ddfra <-
				diag(startupGlobal$dfra)[!effects[effects$include,"basicRate"]]/ngroups
		}
		else
		{
			indx <- which(!effects[effects$include,"basicRate"]) # has length 1
			startup1Model$ddfra <- startupGlobal$dfra[indx,indx]/ngroups
		}
		flush.console()
	# roughly estimate varying parameters
		nActors[i] <- length(data[[i]]$nodeSets[[1]])
		# "use" will indicate the coordinates of the global parameters
		# used for this group,
		# and "rates" indicates the rate parameters among these.
		use <- rep(TRUE, startupGlobal$pp)
		rates <- startupGlobal$requestedEffects$basicRate
		use[rates] <- FALSE
		use[startupGlobal$requestedEffects$group==i] <- TRUE
		rates <- rates[use]
		# project model specification to the single group:
		if (ngroups <= 1)
		{
			effects0 <- effects
			startup1 <- startupGlobal
		}
		else
		{
			effects0 <- getEffects(data[[i]])
			# The following function copies the model specification,
			# except for the basic rate effects for other groups,
			# copies estimated values in startupGlobal to initial values in effects0,
			# and fixes non-random effects.
			effects0 <- projectEffects(effects0,effects,startupGlobal,i)
			cat("Estimate initial parameters group", i, "\n")
			startup1 <- siena07(startup1Model, data=data[[i]],
					effects=effects0, batch=TRUE, silent=silentstart)
			if (initML)
			{	# hybrid procedure: phase1 ML, phase3 MoM.
				startup1 <- siena07(startup1Model3, data=data[[i]],
					effects=effects0, batch=TRUE, silent=silentstart, prevAns=startup1)
			}
			cat("\nInitial estimate obtained\n",
				noquote(sprintf("%5.3f",startup1$theta)),"\n")
		}
		initialEstimates[i,use] <- startup1$theta
		factorsBasicRate[i,use] <- (1/nActors[i])*startup1$theta
		# non-rate parameters will be dropped from this object, hence the name
		rateParameters[,i] <- trafo(startup1$theta[rates])
		# estimate uncertainties for this group:

		covmati <- chol2inv(chol(
			w*precision(startup1) +
			((1-w)/ngroups)*prec.PA[use,use,drop=FALSE]))
		# now covmati is the roughly estimated covariance matrix
		# of the parameter estimate for group i.
		# Transform this to the trafo scale for the rate parameters:
		# double t() in view of recycling rule for matrix arithmetic
		covmati[,rates] <- t(t(covmati[,rates])*devtrafo(startup1$theta[rates]))
		covmati[rates,] <- covmati[rates,]*devtrafo(startup1$theta[rates])
		covmati[rates,!rates] <- 0
		covmati[!rates,rates] <- 0
		# truncate in case rate parameters have a too large variance
		largeRates <- rates & (diag(covmati) > rateVarScale)
		if (any(largeRates))
		{
			rfactor <- sqrt(rateVarScale/(diag(covmati)[largeRates]))
		# double t() in view of recycling rule for matrix arithmetic
			covmati[,largeRates] <- t(t(covmati[,largeRates])*rfactor)
			covmati[largeRates,] <- covmati[largeRates,]*rfactor
		}
		# For the proposal, we do not need the fixed effects
		# (which now includes the non-varying effects)
		z$proposalCov[[i]] <- covmati[
				!effects0$fix[effects0$include], !effects0$fix[effects0$include]]
	} # end for (i in 1:ngroups)

	z$effectName <- effects0$effectName[effects0$include]
	z$effectDepvar <- effects0$name[effects0$include]
	z$requestedEffects <- startup1$requestedEffects

   	z$FRAN <- getFromNamespace(algo$FRANname, pkgname)
    z <- initializeFRAN(z, algo, data=data, effects=effects,
                prevAns=prevAns, initC=FALSE, onlyLoglik=TRUE)
	z$thetaMat <- initialEstimates
	z$basicRate <- z$requestedEffects$basicRate
	z$fix <- z$requestedEffects$fix
    z$nGroup <- z$f$nGroup

	if (max(abs(initialEstimates[,(!z$basicRate)&(!z$fix)]), na.rm=TRUE) > 40)
	{
		cat("\nThe initial estimates are\n")
		print(initialEstimates)
		cat("\nThe largest absolute value is ",
				max(abs(initialEstimates), na.rm=TRUE),
				"\nwhich is too large for good operation of this function.\n",
				"\nAdvice: use a smaller value of initgainGroupwise (perhaps 0).\n\n")
		save(startupGlobal, initialEstimates, file="startupGlobal.RData")
		cat("startupGlobal and initialEstimates saved.\n")
		stop("Divergent initial estimate.")
	}

	# Further down the initial parameter values for the rate parameters are
	# truncated depending on the prior for the rates.
	# This might also be done for other parameters. Not done now.

  	groupPeriods <- attr(z$f, "groupPeriods")
    netnames <- z$f$depNames
	# Prepare some objects used for bookkeeping:
	# ? Can't the data parameter be omitted in the following:
	z$rateParameterPosition <-
        lapply(1:z$nGroup, function(i, periods, data)
           {
               lapply(1:periods[i], function(j)
                  {
                      rateEffects <-
                          z$requestedEffects[z$requestedEffects$basicRate &
                                    z$requestedEffects$period == j &
                                    z$requestedEffects$group == i,]
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
	# Indicator of parameters in the parameter vector of length pp,
	# for which the hierarchical model applies,
	# which is the set of all parameters without the basic rate parameters:
#	z$generalParameters <- (z$requestedEffects$group == 1) & (!z$basicRate)
	# Dependent variable:
#	z$dependentVariable <- z$requestedEffects$name

	# Number of parameters in each group:
	z$TruNumPars <- sum( !z$basicRate ) + length(z$ratePositions[[1]])
	z$TruNumParsPlus <- z$TruNumPars
	if (priorRatesFromData < 0)
	{
		z$TruNumPars <- sum(!z$basicRate )
	}
	scf <- 2.38/sqrt(z$TruNumPars)
# theoretically optimal value according to Roberts & Rosenthal, 2001
    z$scaleFactors <- rep(scf, z$nGroup)
	if (z$p2 > 0)
	{
		z$scaleFactorEta <- 2.38/sqrt(z$p2)
	}
	else
	{
		z$scaleFactorEta <- 1
	}
	z$scaleFactor0 <- z$scaleFactorEta

	# Construction of indicator of parameters of the hierarchical model
	# within the groupwise parameters:
	# Question: Why does this use effects, and not requestedEffects?

	vec1 <- 4 - effects$randomEffects[effects$include]
	vec1[effects$fix[effects$include]] <- 5
	vec1[z$basicRate] <- 2
	vec1[z$ratePositions[[1]]] <- 1
	# now 1=rates group 1; 2=other rates; 3=randomly varying;
	# 4 = estimated non-varying (eta); 5 = rate & (non-estimated, i.e. fixed).
	# set1 = (1,2,3) set2 = (4)
	# Note that fixed rates still have code 1 or 2.
	z$generalParametersInGroup <- vec1[vec1!= 2] %in% c(3,4,5)
	z$varyingGeneralParametersInGroup <- vec1[vec1!= 2] == 3
	z$varyingParametersInGroup <- vec1[vec1!= 2] %in% c(1,3)
	if (incidentalBasicRates)
	{
		z$randomParametersInGroup <- vec1[vec1!= 2] == 3
	}
	else
	{
		z$randomParametersInGroup <- z$varyingParametersInGroup
	}
	z$varyingNonRateInEstimated <- vec1[vec1 %in% c(1,3,4)] == 3
	z$varyingInEstimated <- vec1[vec1 %in% c(1,3,4)] %in% c(1,3)
	z$objectiveInVarying <- vec1[vec1 %in% c(1,3)] == 3
	z$ratesInVarying <- !z$objectiveInVarying
	z$varyingObjectiveParameters <- vec1 == 3

	if (frequentist)
	{
		if (z$nGroup <= 1)
		{
			cat("\nFrequentist option for sienaBayes is available only for ")
			cat("more than one group.\n")
			stop('Frequentist option unavailable.')
		}
#		else
#		{
#	# Approximate sensitivity of expected values for parameters
#			z$dmuinv <- matrix(0, z$TruNumPars, z$TruNumPars)
#			diag(z$dmuinv)[rates] <- 1
#			z$dmuinv[!rates, !rates] <-
#				chol2inv(chol(prec.PA
#					[z$generalParameters,z$generalParameters]))
#			z$dmuinv <- diag(1, z$TruNumPars, z$TruNumPars)
#		}
	}
	# Determination parameters of prior distribution;
	# this specifies the defaults!
	if (is.null(priorMu))
	{
		z$priorMu <- matrix( c(0), z$p1, 1 )
		if ((priorRatesFromData >= 0) & (!incidentalBasicRates))
		{
			z$priorMu[z$ratesInVarying] <- 2
		}
	}
	else
	{
		if (length(priorMu) != z$p1)
		{
			stop(paste("priorMu must have length ",z$p1))
		}
		z$priorMu <- matrix(priorMu, z$p1, 1 )
	}
	z$priorKappa <- 1
	if (!is.null(priorKappa))
	{
		if (priorKappa > 0)
		{
			z$priorKappa <- priorKappa
		}
	}
	z$priorDf <- z$p1 + 2
	if (!is.null(priorDf))
	{
		if (priorDf > 0)
		{
			z$priorDf <- priorDf
		}
	}
	if (z$priorDf + z$nGroup <= z$p1)
	{
		z$priorDf <- z$p1 - z$nGroup + 2
	}

		z$priorSigma <- diag(z$p1)
	if (!is.null(priorSigma))
	{
		if ((class(priorSigma)=="matrix") &
			all(dim(priorSigma)==c(z$p1,z$p1)))
		{
		z$priorSigma <- priorSigma
	}
	else
	{
		stop(paste("priorSigma is not a matrix of dimensions",z$p1))
		}
	}
	if (priorRatesFromData == 2)
	{
	# Rate parameters are nuisance parameters, and get a data-induced prior;
	# robust estimates for location and scale plus a small ridge.
		## require(MASS)
		if (ncol(rateParameters) <= nrow(rateParameters))
		{
			try.error <- (inherits(try(
						robustRates <- cov.trob(t(rateParameters)),
								silent=TRUE), "try-error"))
		}
		else
		{
			try.error <- (inherits(try(
					robustRates <- cov.mcd(t(rateParameters), nsamp="sample"),
								silent=TRUE), "try-error"))
		}
		if (try.error)
		{
			cat('Condition priorRatesFromData=2 impossible, changed to 1.\n')
			priorRatesFromData <- 1
		}
		else
		{
			z$priorMu[z$ratesInVarying] <- robustRates$center
			z$priorSigma[z$ratesInVarying,] <- 0
			z$priorSigma[,z$ratesInVarying] <- 0
			z$priorSigma[z$ratesInVarying,z$ratesInVarying] <-
				reductionFactor * robustRates$cov
			diag(z$priorSigma)[z$ratesInVarying] <-
				(diag(z$priorSigma)[z$ratesInVarying] + 0.5)
		}
	}
	if (priorRatesFromData == 1)
	{
	# Rate parameters are nuisance parameters, and get a data-induced prior;
	# observed covariance matrix plus a small ridge.
		z$priorMu[z$ratesInVarying] <- rowMeans(rateParameters)
		z$priorSigma[z$ratesInVarying,] <- 0
		z$priorSigma[,z$ratesInVarying] <- 0
		z$priorSigma[z$ratesInVarying,z$ratesInVarying] <-
			reductionFactor * cov(t(rateParameters))
		diag(z$priorSigma)[z$ratesInVarying] <-
			(diag(z$priorSigma)[z$ratesInVarying] + 0.5)
	}
	# Truncate initial parameter values for the rate parameters
	# depending on the prior for the rates:
	if ((priorRatesFromData >= 0) & (!incidentalBasicRates)) # else all rates are fixed or incidental
	{
	z$thetaMat[,z$basicRate] <- pmin(z$thetaMat[,z$basicRate],
			z$priorMu[z$ratesInVarying] + 2*sqrt(diag(z$priorSigma)[z$ratesInVarying]))
	for (group in 1:z$nGroup)
	{
		z$thetaMat[group, z$ratePositions[[group]]] <-
						pmin(z$thetaMat[group, z$ratePositions[[group]]],
		z$priorMu[z$ratesInVarying] + 2*sqrt(diag(z$priorSigma)[z$ratesInVarying]))
	}
	}

	z$incidentalBasicRates <- incidentalBasicRates
	z$priorRatesFromData <- priorRatesFromData
	z$delta <- delta
	z$gamma <- gamma
	z$reductionFactor <- reductionFactor
	if (nwarm < 5){nwarm <<- 5}
	if (nmain < 5 + nwarm){nmain <<- 5+nwarm}

	if (frequentist)
	{
		lengthPhase1 <- max(lengthPhase1, 2)
		lengthPhase3 <- max(lengthPhase3, 2)
	}
	else
	{
		lengthPhase1 <- round(nmain/5)
		lengthPhase3 <- round(nmain/5)
	}
	if ((nwarm + lengthPhase1 + 5 + lengthPhase3) > nmain)
	{
		nwarm <- max(5,
			round(nwarm*nmain/(nwarm + 2*lengthPhase1+lengthPhase3)))
		oldLengthPhase1 <- lengthPhase1
		lengthPhase1 <- max(5,
			round(lengthPhase1*nmain/(nwarm + 2*lengthPhase1+lengthPhase3)))
		lengthPhase3 <- max(5,
			round(lengthPhase3*nmain/(nwarm + 2*oldLengthPhase1+lengthPhase3)))
		nmain <<- nwarm + 2*lengthPhase1 + lengthPhase3
		cat("Iteration numbers adapted:\n")
		cat("nwarm = ", nwarm, "; nmain = ", nmain)
		if (frequentist)
		{
			cat(", lengthPhase1 = ", lengthPhase1,
				"; lengthPhase3 = ", lengthPhase3, ".\n")
		}
		else
		{
			cat(".\n")
		}

	}
	z$nwarm <- nwarm
	z$nmain <- nmain
	z$nImproveMH <- nImproveMH
	z$priorRatesFromData <- priorRatesFromData
	z$lengthPhase1 <- lengthPhase1
	z$lengthPhase3 <- lengthPhase3

	# Now generalParametersInGroup is a logical vector of length TruNumPars,
	# indicating the non-basicRate parameters.

	if (sum(z$basicRate) != z$nGroup*length(z$ratePositions[[1]]))
		{
		stop("Function sienaBayes does not allow non-basic rate parameters.")
		}
	# thetaMat must get some NAs so that in sampleVaryingParameters it will be
	# clear which values are meaningless and will not be operated upon.
	# These have to always stay NA.
    for (i in 1:z$nGroup)
    {
        use <- rep(FALSE, z$pp)
        use[z$ratePositions[[i]]] <- TRUE
        use[!z$basicRate] <- TRUE
        z$thetaMat[i, !use] <- NA
		factorsBasicRate[i, !use] <- NA
    }
	z$factorsBasicRate <- factorsBasicRate[,z$basicRate]
	z$factorsEta <- (z$initialSE[z$set2])^2

#cat("factorsBasicRate\n", z$factorsBasicRate,"\n")
#cat("factorsEta\n", z$factorsEta,"\n")
	meanThetaMat <- colMeans(matrix( z$thetaMat[!is.na(z$thetaMat)],
									z$nGroup, z$TruNumParsPlus))
	if (priorRatesFromData < 0)
	{
		z$muTemp <- matrix(meanThetaMat[z$varyingGeneralParametersInGroup], z$p1, 1)
	}
	else
	{
	meanThetaMat <- colMeans(matrix( z$thetaMat[!is.na(z$thetaMat)],
									z$nGroup, z$TruNumPars))
		z$muTemp <- matrix(meanThetaMat[z$randomParametersInGroup], z$p1, 1)
	}
	z$SigmaTemp <- diag(z$p1)

	is.batch(TRUE)

    if (nbrNodes > 1 && z$observations > 1)
    {
        ## require(parallel)
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
        clusterCall(z$cl, initializeFRAN, z, algo,
                    initC = TRUE, profileData=FALSE, returnDeps=FALSE)
		clusterSetRNGStream(z$cl, iseed = as.integer(runif(1,
								max=.Machine$integer.max)))
    }
    ## z$returnDataFrame <- TRUE # chains come back as data frames not lists
    z$returnChains <- FALSE
    z
} # end initializeBayes

##@projectEffects used in initializeBayes to project model specification
# modified version of updateTheta
# model specification of prevAnss is transferred to effs for group ii
# and !randomEffects are fixed.
# Since interactions are defined by referring for the interacting effects
# to numbers in the effects object,
# and the within-group effects object has a different numbering
# than the multi-group effects object,
# the main effects are treated in a more easy way,
# and the interaction effects are treated specially.
projectEffects <- function(effs, theEffects, prevAnss, ii)
{
	# take out from effs the default effects in the objective function.
	# these will be added later on if they are included in prevAnss$effects
	# (which normally will be the case, but not always).
	effs[effs$shortName=='density','include'] <- FALSE
	effs[effs$shortName=='recip','include'] <- FALSE
	effs[effs$shortName=='linear','include'] <- FALSE
	effs[effs$shortName=='quad','include'] <- FALSE

	prevEffects <- prevAnss$effects
	prevReqEffects <- prevAnss$requestedEffects
	prevefflist <- apply(prevEffects, 1, function(x)
					paste(x[c("name", "shortName",
							"type",
							"interaction1", "interaction2",
							"period")],
						collapse="|"))
	reqefflist <- apply(prevReqEffects, 1, function(x)
					paste(x[c("name", "shortName",
							"type",
							"interaction1", "interaction2",
							"period")],
						collapse="|"))
	includeEff <-  prevefflist %in% reqefflist
# this is the list of requested effects;
# if any interactions are included in prevEffects without some
# corresponding main effects then prevEffects contains, in addition,
# those main effects .
	prevReqEffects$initialValue <- prevAnss$theta
# this is the same as prevReqEffects <- updateTheta(prevEffects, prevAnss)
# but without also including the interactions:
# since interacting effects have identifying numbers (effects1, effects2, effects3)
# that differ between prefReqEffects and effs,
# they must be treated separately.
	oldInteract <- (prevReqEffects$shortName %in% c("unspInt" , "behUnspInt"))
	oldlist <- apply(prevReqEffects[(!oldInteract),], 1, function(x)
					paste(x[c("name", "shortName",
							"type",
							"interaction1", "interaction2",
							"period")],
						collapse="|"))
	efflist <- apply(effs, 1, function(x)
					paste(x[c("name", "shortName",
							"type",
							"interaction1", "interaction2",
							"period")],
						collapse="|"))
	use <- efflist %in% oldlist
	# first define the include column for the main effects
	effs$include[use] <-
		prevReqEffects$include[match(efflist, oldlist)][use]
	# this is all TRUE;
	# these are only for the main effects (!oldInteract);
	# now treat the include column for the interactions
	if (sum(oldInteract) >= 1)
	{
		# get the characteristics of the interactions
		Names <- prevReqEffects$name[oldInteract]
		Types <- prevReqEffects$type[oldInteract]
		Parameters <- prevReqEffects$parm[oldInteract]
		effects1 <- prevReqEffects$effect1[oldInteract]
		effects2 <- prevReqEffects$effect2[oldInteract]
		effects3 <- prevReqEffects$effect3[oldInteract]
		# look up the interacting main effects
		# Note that some of effects3 may be 0,
		# and theEffects[0,] does not exist.
		shortNames1 <- theEffects[effects1,'shortName']
		shortNames2 <- theEffects[effects2,'shortName']
		shortNames3 <- rep('', length(effects3))
		shortNames3[(effects3 > 0)] <-
						theEffects[effects3[(effects3 > 0)],'shortName']
		interactions11 <- theEffects[effects1,'interaction1']
		interactions12 <- theEffects[effects2,'interaction1']
		interactions13 <- rep('', length(effects3))
		interactions13[(effects3 > 0)] <-
						theEffects[effects3[(effects3 > 0)],'interaction1']
		interactions21 <- theEffects[effects1,'interaction2']
		interactions22 <- theEffects[effects2,'interaction2']
		interactions23 <- rep('', length(effects3))
		interactions23[(effects3 > 0)] <-
						theEffects[effects3[(effects3 > 0)],'interaction2']
		# put the interactions in the groupwise effects object
		for (i in seq_along(effects1))
		{
			if (effects3[i] == 0)
			{
				effs <- includeInteraction(effs,
					shortNames1[i], shortNames2[i], name=Names[i],
					interaction1=c(interactions11[i], interactions12[i]),
					interaction2=c(interactions21[i], interactions22[i]),
					type=Types[i], parameter=Parameters[i],
					character=TRUE, verbose=FALSE)
			}
			else
			{
				effs <- includeInteraction(effs,
		shortNames1[i], shortNames2[i], shortNames3[i], name=Names[i],
		interaction1=c(interactions11[i], interactions12[i], interactions13[i]),
		interaction2=c(interactions21[i], interactions22[i], interactions23[i]),
					type=Types[i], parameter=Parameters[i],
					character=TRUE, verbose=FALSE)
			}
		}
	}
	# now continue with columns initialValue and fix effs.
	effs$initialValue[effs$include & (!effs$basicRate)]	<-
			prevReqEffects$initialValue[
		prevReqEffects$include & (!prevReqEffects$basicRate)]
	effs$initialValue[effs$basicRate] <- prevReqEffects$initialValue[
		(prevReqEffects$basicRate) & (prevReqEffects$group == ii)]
	effs$fix[effs$include & (!effs$basicRate)] <-
		(!prevReqEffects$randomEffects | prevReqEffects$fix)[
		prevReqEffects$include & (!prevReqEffects$basicRate)]
	effs$fix[effs$basicRate] <-
		prevReqEffects$fix[prevReqEffects$basicRate & (prevReqEffects$group == ii)]
	effs
} # end projectEffects
##@endBayes algorithms do set up for Bayesian model

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

##@dmvnorm algorithms calculates log multivariate normal density
##inefficient: should not call mahalanobis and eigen with same sigma repeatedly
dmvnorm <- function(x, mean , sigma)
{
    if (is.vector(x))
    {
        x <- matrix(x, ncol=length(x))
    }
	evs <- eigen(sigma, symmetric=TRUE, only.values=TRUE)$values
	if (min(evs)< 1e-8)
	{
		cat('singular covariance matrix in dmvnorm\n')
		dmvn <- -Inf
	}
	else
	{
		distval <- mahalanobis(x, center=mean, cov=sigma)
		logdet <- sum(log(evs))
		dmvn <- -(ncol(x) * log(2 * pi) + logdet + distval) / 2
	}
	dmvn
}

##@dmvnorm2 algorithms calculates log multivariate normal density
##more efficient: uses previously computed ingredients for mahalanobis and eigen
dmvnorm2 <- function(x, mean , evs, sigma.inv)
{
    if (is.vector(x))
    {
        x <- matrix(x, ncol=length(x))
    }
	if (min(evs)< 1e-8)
	{
		cat('singular covariance matrix in dmvnorm2\n')
		dmvn <- -Inf
	}
	else
	{
		distval <- mahalanobis(x, center=mean, cov=sigma.inv, inverted=TRUE)
		logdet <- sum(log(evs))
		dmvn <- -(ncol(x) * log(2 * pi) + logdet + distval) / 2
	}
	dmvn
}

##@getProbabilitiesFromC sienaBayes gets loglik from chains in C
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
        ans <- .Call(C_getChainProbabilities, PACKAGE = pkgname, f$pData,
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
# It was the following - must be wrong
#	ans[[1]] <- sum(sapply(anss, "[[", 1))
	if (nrow(callGrid) != length(anss))
	{
		cat("Error: nrow(callGrid) = ", nrow(callGrid), "; length(anss) = ",
			length(anss),"\n")
		stop("Error in getProbabilitiesFromC")
	}
	# Sum the log probabilities for the periods corresponding to each group
	logprob <- sapply(anss, "[[", 1)
	ans[[1]] <- sapply(1:z$nGroup, function(i){sum(logprob[callGrid[,1]==i])})
	if (getScores)
	{
		ans[[2]] <- sapply(anss, "[[", 2) # it was rowSums of this
	}
	ans[[3]] <- sapply(anss, "[[", 3)
	ans
}

##@doGetProbabilitiesFromC Maximum likelihood
doGetProbabilitiesFromC <- function(x, thetaMat, index, getScores)
{
# thetaMat <- z$thetaMat # [1,]
# getScores <- TRUE
# x <- c(1,1)
    f <- FRANstore()
	theta <- thetaMat[x[1], ]
#	gcp <-
    .Call(C_getChainProbabilities, PACKAGE = pkgname, f$pData,
		  f$pModel, as.integer(x[1]), as.integer(x[2]),
		  as.integer(index), f$myeffects, theta, getScores)
}

##@glueBayes combines two sienaBayesFit into one
glueBayes <- function(z1,z2,nwarm2=0){
	z <- list()
	dif <- FALSE
	nstart2 <- nwarm2+1
	d1 <- sum(!is.na(z1$ThinPosteriorMu[,1]))
	d2 <- sum(!is.na(z2$ThinPosteriorMu[,1]))
	if (nstart2 > d2){stop("Insufficient data in z2.")}
	z$ThinPosteriorMu <-
		rbind(z1$ThinPosteriorMu[1:d1,,drop=FALSE],
				z2$ThinPosteriorMu[nstart2:d2,,drop=FALSE])
	z$ThinPosteriorEta  <-
		rbind(z1$ThinPosteriorEta[1:d1,,drop=FALSE],
				z2$ThinPosteriorEta[nstart2:d2,,drop=FALSE])
	z$ThinPosteriorSigma  <-
		array(NA, c(d1+d2-nwarm2, dim(z1$ThinPosteriorSigma)[c(2,3)]))
	z$ThinPosteriorSigma[1:d1,,] <- z1$ThinPosteriorSigma[1:d1,,,drop=FALSE]
	z$ThinPosteriorSigma[(d1+1):(d1+d2-nwarm2),,] <-
		z2$ThinPosteriorSigma[nstart2:d2,,,drop=FALSE]
	z$ThinParameters  <-
		array(NA, c(d1+d2-nwarm2, dim(z1$ThinParameters)[c(2,3)]))
	z$ThinParameters[1:d1,,] <- z1$ThinParameters[1:d1,,,drop=FALSE]
	z$ThinParameters[(d1+1):(d1+d2-nwarm2),,] <-
		z2$ThinParameters[nstart2:d2,,,drop=FALSE]
	z$ThinBayesAcceptances <-
		rbind(z1$ThinBayesAcceptances[1:d1,,drop=FALSE],
				z2$ThinBayesAcceptances[nstart2:d2,,drop=FALSE])
	z$frequentist <- z1$frequentist
	z$pp <- z1$pp
	z$cconditional <- z1$cconditional
	z$rate <- z1$rate
	z$nImproveMH  <- z1$nImproveMH
	z$nwarm  <- z1$nwarm
	z$nmain  <- z1$nmain + z2$nwarm+z2$nmain - nwarm2
	what <- ''
	if (all(z1$mult == z2$mult))
	{
		z$mult <- z1$mult
	}
	else
	{
		dif <- TRUE
		what <- 'mult;'
	}
	if (z1$nrunMHBatches == z2$nrunMHBatches)
	{
		z$nrunMHBatches <- z1$nrunMHBatches
	}
	else
	{
		dif <- TRUE
		what <- paste(what, 'nrunMHBatches;')
	}
	if (z1$nSampConst == z2$nSampConst)
	{
		z$nSampConst <- z1$nSampConst
	}
	else
	{
		dif <- TRUE
		what <- paste(what, 'nSampConst;')
	}
	if (z1$nSampVarying == z2$nSampVarying)
	{
		z$nSampVarying <- z1$nSampVarying
	}
	else
	{
		dif <- TRUE
		what <- paste(what, 'nSampVarying;')
	}
	if (z1$nSampRates == z2$nSampRates)
	{
		z$nSampRates <- z1$nSampRates
	}
	else
	{
		dif <- TRUE
		what <- paste(what, 'nSampRates;')
	}
	if (all(z1$fix == z2$fix))
	{
		z$fix <- z1$fix
	}
	else
	{
		dif <- TRUE
		what <- paste(what, 'fix')
	}
	if (all(z1$effectName == z2$effectName))
	{
		z$effectName <- z1$effectName
	}
	else
	{
		dif <- TRUE
		what <- paste(what, 'effectName;')
	}
	if (all(z1$effectDepvar == z2$effectDepvar))
	{
		z$effectDepvar <- z1$effectDepvar
	}
	else
	{
		dif <- TRUE
		what <- paste(what, 'effectDepvar;')
	}
	if (all(z1$groupNames == z2$groupNames))
	{
		z$groupNames <- z1$groupNames
	}
	else
	{
		dif <- TRUE
		what <- paste(what, 'groupNames;')
	}
	if (all(z1$ratesInVarying == z2$ratesInVarying))
	{
		z$ratesInVarying <- z1$ratesInVarying
	}
	else
	{
		dif <- TRUE
		what <- paste(what, 'ratesInVarying;')
	}
	if (z1$TruNumPars == z2$TruNumPars)
	{
		z$TruNumPars <- z1$TruNumPars
	}
	else
	{
		dif <- TRUE
		what <- paste(what, 'TruNumPars;')
	}
	if (dif)
	{
		stop(paste("The two objects do not have the same specification. Difference in:", what))
	}
	if (all(z1$priorMu == z2$priorMu))
	{
		z$priorMu <- z1$priorMu
	}
	else
	{
		dif <- TRUE
		what <- paste(what, 'priorMu;')
	}
	if (all(z1$priorSigma == z2$priorSigma))
	{
		z$priorSigma <- z1$priorSigma
	}
	else
	{
		dif <- TRUE
		what <- paste(what, 'priorSigma;')
	}
	if (z1$priorKappa == z2$priorKappa)
	{
		z$priorKappa <- z1$priorKappa
	}
	else
	{
		dif <- TRUE
		what <- paste(what, 'priorKappa;')
	}
	if (z1$priorDf == z2$priorDf)
	{
		z$priorDf <- z1$priorDf
	}
	else
	{
		dif <- TRUE
		what <- paste(what, 'priorDf;')
	}
	if (z1$priorRatesFromData == z2$priorRatesFromData)
	{
		z$priorRatesFromData <- z1$priorRatesFromData
	}
	else
	{
		dif <- TRUE
		what <- paste(what, 'priorRatesFromData;')
	}
	if (z1$incidentalBasicRates == z2$incidentalBasicRates)
	{
		z$incidentalBasicRates <- z1$incidentalBasicRates
	}
	else
	{
		dif <- TRUE
		what <- paste(what, 'incidentalBasicRates;')
	}
	if (dif)
	{
		stop(paste("The two objects do not have the same prior. Difference in:", what))
	}
	z$initialResults <- z1$initialResults
	z$nGroup <- z1$nGroup
	z$set1  <- z1$set1
	z$set2  <- z1$set2
	z$p1  <- z1$p1
	z$p2  <- z1$p2
	z$f <- z1$f
	z$effects <- z1$effects
	z$requestedEffects <- z1$requestedEffects
	z$basicRate  <- z1$basicRate
	z$ratePositions  <- z1$ratePositions
	z$objectiveInVarying  <- z1$objectiveInVarying
	z$generalParametersInGroup <- z1$generalParametersInGroup
	z$varyingParametersInGroup  <- z1$varyingParametersInGroup
	z$randomParametersInGroup  <- z1$randomParametersInGroup
	z$varyingGeneralParametersInGroup  <- z1$varyingGeneralParametersInGroup
	z$varyingObjectiveParameters <- z1$varyingObjectiveParameters
	z$incidentalBasicRates <- z1$incidentalBasicRates
	class(z) <- "sienaBayesFit"
	z
}

##@protectedInverse inverse of p.s.d matrix
protectedInverse <- function(x)
	{
		if (inherits(try(xinv <- chol2inv(chol(x)),
			silent=TRUE), "try-error"))
		{
	# Now make this x positive definite, if it is not.
	# See above for a more extensive treatment of the same.
	# Adapted from function make.positive.definite in package corpcor
	# which uses a method by Higham (Linear Algebra Appl 1988)
	# but changed to make the matrix positive definite (psd)
	# instead of nonnegative definite.
	# The idea is to left-truncate all eigenvalues to delta0.
	# The construction with tol, not used now,
	# is to ensure positive definiteness given numerical inaccuracy.
			es <- eigen(x)
			esv <- es$values
			delta0 <- 1e-6
			cat("protectedInverse: Eigenvalues Sigma = ", sort(esv), "\n")
			if (min(esv) < delta0)
			{
				delta <- delta0
				tau <- pmax(delta, esv)
		#		cat("Smallest eigenvalue of Sigma now is ",
		#					min(esv),"; make posdef.\n")
				xinv <- es$vectors %*%
							diag(1/tau, dim(x)[1]) %*% t(es$vectors)
			}
		}
		xinv
	}


##@trafo link function rates
# the log link is too strong!
# If the trafo is changed, then in print.summary.sienaBayesFit
# the note about the scale of the basic rate parameters
# must be changed correspondingly.
#trafo <- function(x){ifelse(x<0.01, 0.1, sqrt(x))}
trafo <- function(x){ifelse(x<0.01, 0.01, x)}
##@antitrafo inverse link function rates
#antitrafo <- function(x){x^2}
antitrafo <- function(x){x}
##@devtrafo derivative link function rates
#devtrafo <- function(x){ifelse(x<0.01, 0.05, 1/(2*sqrt(x)))}
devtrafo <- function(x){1}

