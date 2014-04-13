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
				priorMu=NULL, priorSigma=NULL, priorDf=NULL, priorKappa=NULL,
				frequentist=FALSE, incidentalBasicRates=FALSE,
				gamma=0.5, delta=1e-4,
				nwarm=50, nmain=250, nrunMHBatches=20,
				lengthPhase1=round(nmain/5), lengthPhase3=round(nmain/5),
				storeAll = FALSE, prevAns=NULL,
				prevBayes = NULL, silentstart=TRUE,
				plotit=FALSE, nbrNodes=1, clusterType=c("PSOCK", "FORK"),
				getDocumentation=FALSE)
{
    ##@createStores internal sienaBayes Bayesian set up stores
    createStores <- function()
    {
        npar <- length(z$theta)
		numberRows <- nmain * nrunMHBatches
        z$candidates <<- array(NA, dim=c(numberRows, z$nGroup, npar))
        z$acceptances <<- matrix(NA, nrow=numberRows, ncol=z$nGroup)
		if (storeAll)
		{
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
			z$StoreLoglik <<- matrix(NA, numberRows, z$nGroup)
		}
		# additional thinned stores
		z$ThinPosteriorMu <<- matrix(NA, nrow=nwarm+nmain, ncol=z$p1)
		z$ThinPosteriorEta <<- matrix(NA, nrow=nwarm+nmain, ncol=z$p2)
		z$ThinPosteriorSigma <<-
					array( NA, dim=c(nwarm+nmain, z$p1, z$p1) )
		z$ThinLoglik <<- matrix(NA, nwarm+nmain, z$nGroup)
		z$ThinParameters <<- array(NA,
						dim=c(nwarm+nmain, z$nGroup, z$TruNumPars))
		z$ThinBayesAcceptances <<-
						matrix(NA, nrow=nwarm+nmain, ncol=z$nGroup+2)
		z$sumBasicRates <<- matrix(0, nrow=z$nGroup, ncol=sum(z$basicRate))
		if (z$incidentalBasicRates)
		{
			z$ThinScores <<- matrix(NA, nrow=nwarm+nmain, ncol=sum(z$basicRate))
		}
		z$ThinEtaScores <<- array(NA, dim=c(nwarm+nmain, z$p2, z$nGroup))
		z$ThinEtaHessian <<- array(NA, dim=c(nwarm+nmain, z$p2, z$p2))
		z$ThinP3MuHat <<- matrix(NA, z$lengthPhase3, z$p1)
		z$ThinP3SigmaHat <<- array(NA,
							dim=c(z$lengthPhase3, z$p1, z$p1))
		z$ThinP3EtaScores <<- matrix(NA, z$lengthPhase3, z$p2)
	} # end createStores

    ##@storeData internal sienaBayes; put data in stores
    storeData <- function(cycleData)
    {
        start <- z$sub + 1
        nrun <- ncol(cycleData$accepts)
        end <- start + nrun - 1
		if (storeAll)
		{
			z$acceptances[start:end, ] <<- cycleData$accepts
			z$StorePosteriorMu[start:end, ] <<- cycleData$posteriorMu
			z$StorePosteriorEta[start:end, ] <<- cycleData$posteriorEta
			z$StorePosteriorSigma[start:end, , ] <<-  cycleData$PosteriorSigma
			z$StoreLoglik[start:end,] <<- cycleData$loglikVec
		}
		# Also store a thinned version of the process,
		# containing only the final values in the main steps:
        z$ThinPosteriorMu[z$iimain, ] <<- cycleData$posteriorMu[nrun,]
        z$ThinPosteriorEta[z$iimain, ] <<- cycleData$posteriorEta[nrun,]
		z$ThinPosteriorSigma[z$iimain, , ] <<-  cycleData$PosteriorSigma[nrun,,]
        z$ThinLoglik[z$iimain,] <<- cycleData$loglikVec[nrun,]
		z$ThinBayesAcceptances[z$iimain, ] <<- cycleData$BayesAcceptances
		if (z$incidentalBasicRates)
		{
			z$ThinScores[z$iimain, ] <<- cycleData$thescores
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
			z$MHacceptances[start:end, , , ] <<- cycleData$MHaccepts
			z$MHrejections[start:end, , , ] <<- cycleData$MHrejects
			z$MHproportions[start:end, , , ] <<- cycleData$MHaccepts /
				(cycleData$MHaccepts + cycleData$MHrejects)
        }
		z$ThinEtaScores[z$iimain, ,] <<-
				t(cycleData$ans.last$fra)[z$set2,,drop=FALSE]
		z$ThinEtaHessian[z$iimain, ,] <<-
				as.double(cycleData$ans.last$dff[z$set2,z$set2,drop=FALSE])
        z$sub <<- z$sub + nrun
    } # end storeData

	##@improveMH internal sienaBayes; find scale factors
    improveMH <- function(tiny=1.0e-15, totruns=100, desired=totruns/4,
							maxiter=20, tolerance=totruns/17, getDocumentation=FALSE)
    {
	# A rough stochastic approximation algorithm.
	# Two consecutive iterations of "totruns" runs must result
	# for each coordinate of actual, in "desired" acceptances
	# with a tolerance of "tolerance".
		##@rescaleCGD internal improveMH Bayesian
        rescaleCGD <- function(iter)
        {
#			u <- ifelse (actual > desired,
#							2 - ((totruns - actual) / (totruns - desired)),
#							1 / (2 - (actual / desired)))
			u <- 1 + (actual - desired)/
								(sqrt(sqrt(iter)*desired*(totruns - desired)))
			u <- pmin(pmax(u, 0.1), 10)
            number <<- ifelse(abs(actual - desired) <= tolerance,
								number + 1, 0 )
            success <<- (number >= 2)
            u
        }
		if (getDocumentation)
		{
			tt <- getInternals()
			return(tt)
		}
        iter <- 0
		nearGoal <- FALSE
        number <- rep(0, z$nGroup+2)
        success <- rep(FALSE, z$nGroup+2)
		groups <- c(rep(TRUE, z$nGroup), FALSE, FALSE)
		cat('improveMH\n')
        repeat
        {
            iter <- iter + 1
            cyc <- MCMCcycle(nrunMH=totruns, change=FALSE, bgain=-1)
            actual <- cyc$BayesAcceptances
			if (nearGoal)
			# replace actual by exponentially weighted moving average of actual
			{
				actual <- (actual + actual.previous)/2
			}
			else
			{
				if (max(abs(actual - desired)) < 2*tolerance){nearGoal <- TRUE}
			}
			actual.previous <- actual
			## this is a vector of the expected number of acceptances by group
			## and for the non-varying parameters, in the MH change of
			## parameters out of nrunMHBatches trials
			if (!(any(z$set2) & (!frequentist)))
			{
				actual[(z$nGroup+1) : (z$nGroup+2)] <- desired
			}
            multiplier <- rescaleCGD(iter)
            mustUpdate <- (number < 3)
			if (any(is.na(mustUpdate)))
			{
				cat("is.na(mustUpdate) in improveMH\n")
				mustUpdate[is.na(mustUpdate)] <- FALSE
# This should not happen any more,
# as it is caught in sampleVaryingParameters()
			}
            z$scaleFactors[mustUpdate[groups]] <<-
					z$scaleFactors[mustUpdate[groups]] *
								multiplier[mustUpdate & groups]
			if (mustUpdate[z$nGroup+1])
			{
				z$scaleFactor0 <<-
					z$scaleFactor0 * multiplier[z$nGroup+1]
			}
			if (mustUpdate[z$nGroup+2])
			{
				z$scaleFactorEta <<-
					z$scaleFactorEta * multiplier[z$nGroup+2]
			}
            cat('\n',iter, '.         ', sprintf("%4.1f", actual),
				     '\n multipliers  ', sprintf("%4.2f", multiplier),
			         '\n scalefactors ', sprintf("%4.2f", z$scaleFactors),
					sprintf("%4.2f", z$scaleFactor0),
					sprintf("%4.2f", z$scaleFactorEta),'\n')

			flush.console()
            if (all(success) || iter == maxiter)
            {
                break
            }
            if (any(z$scaleFactors < tiny))
            {
				cat('tiny: ', which(z$scaleFactors < tiny),'\n')
				cat('This is a problem: scalefactors too small.\n')
				cat('Try a model with fewer randomly varying effects.\n')
				stop('Tiny scalefactors.')
            }
        }
        cat('fine tuning took ', iter, ' iterations.\n')
#drop		cat('scalefactors:', z$scaleFactors, z$scaleFactor0,
#										z$scaleFactorEta, '\n')
		flush.console()
    } # end improveMH


	##@MCMCcycle internal sienaBayes; some loops of
	## (MH steps and sample parameters)
	MCMCcycle <- function(nrunMH, change=TRUE, bgain)
	{
	# is for the storage in the object produced
	# and for its by-effects on the MCMC state:
	# thetaMat, muTemp, sigmaTemp
		zm <- list()
		zm$accepts <- matrix(NA, nrow=z$nGroup, nrunMH)
		acceptsEta <- rep(NA, nrunMH)
		zm$MHaccepts <- array(NA, dim=c(nrunMH, z$nGroup,
							  z$nDependentVariables, 9))
		zm$MHrejects <- array(NA, dim=c(nrunMH, z$nGroup,
							  z$nDependentVariables, 9))
		zm$MHaborts <- array(NA, dim=c(nrunMH, z$nGroup,
							 z$nDependentVariables, 9))
		zm$posteriorMu <- matrix(0, nrow=nrunMH, ncol=z$p1)
		zm$posteriorEta <- matrix(0, nrow=nrunMH, ncol=z$p2)
		zm$PosteriorSigma <- array( NA, dim=c(nrunMH, z$p1, z$p1))
		zm$loglikVec <- matrix(NA, nrunMH, z$nGroup)
		zm$BayesAcceptances <- rep(NA, z$nGroup+2)

		for (i in 1:nrunMH)
		{
		# This is the hot loop.
		# sample MH steps in the network chain:
			if (i%%20 == 0)
			{
				cat(".")
				flush.console()
			}
			if (i < nrunMH)
			{
				ans <- z$FRAN(z, byGroup=TRUE, returnLoglik=TRUE, onlyLoglik=TRUE)
			}
			else
			{#abcd
				z$Deriv <- TRUE
				ans <- z$FRAN(z, byGroup=TRUE, returnLoglik=TRUE, onlyLoglik=FALSE)
				z$Deriv <- FALSE
			}
			loglik <- ans$loglik  #drop?
			# After the nrunMH the loglik by group is used as logpOld,
			# see below.
			# Perhaps efficiency can be gained here
			# by not calculating loglik all the time.
			# sample the group-level parameters:
			mcmcc.accept <- sampleVaryingParameters(change, bgain)
			if (any(z$set2) & (!frequentist))
			{
				sampleConstantParameters(change, i%%2)
			}
			else
			{
				z$acceptEta <- 0.5
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
				zm$posteriorMu[i, ] <- z$muTemp
				zm$posteriorEta[i,] <- z$thetaMat[1,z$set2]
				zm$PosteriorSigma[i, ,] <- z$SigmaTemp
			}
			zm$loglikVec[i,] <- loglik
			zm$accepts[, i] <- mcmcc.accept
			acceptsEta[i] <- z$acceptEta
			zm$MHaccepts[i, , , ] <-
				t(do.call(cbind,
						  tapply(ans$accepts, factor(z$callGrid[, 1]),
								 function(x)Reduce("+", x))))
			zm$MHrejects[i, , , ] <-
				t(do.call(cbind, tapply(ans$rejects, factor(z$callGrid[, 1]),
										function(x)Reduce("+", x))))
			zm$MHaborts[i, , , ] <- t(do.call(cbind,
										tapply(ans$aborts,
												factor(z$callGrid[, 1]),
												function(x)Reduce("+", x))))
		}
		etaAcceptances0 <- sum(acceptsEta[2*(1:(nrunMH %/% 2))])
		etaAcceptances1 <- sum(acceptsEta[2*(1:(nrunMH %/% 2))-1])
		zm$BayesAcceptances <- c(rowSums(zm$accepts),
					etaAcceptances0, etaAcceptances1)
		zm$ans.last <- ans
		zm
	} # end MCMCcycle

	##@sampleVaryingParameters internal sienaBayes; propose new parameters
	## and accept them or not
	## Replace accepted new parameters only if change.
	## If z$incidentalBasicRates, for the basic rate parameters
	## not a Bayesian MH step but a Robbins Monro step is made.
	sampleVaryingParameters <- function(change=TRUE, bgain)
	{
		require(MASS)

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
					0.01*bgain*z$factorsBasicRate*t(scores[z$basicRate,]), 0.1)
		# Drop the structurally zero elements for storage purposes.
				z$thescores <<- sapply(1:z$nGroup,
					function(i){scores[z$ratePositions[[i]], i]})
			}
		# Basic rate parameters must be masked:
			varyingParameters <- z$varyingObjectiveParameters
			restrictProposal <- z$objectiveInVarying
		}
		else
		{
			varyingParameters <- z$set1
			restrictProposal <- rep(TRUE, z$p1) # no restriction
		}

	 	# do the fandango
		## get a multivariate normal with covariance matrix proposalCov
		## multiplied by a scale factor which varies between groups
		#  propose the changes in the varying parameters
			thetaChanges <- t(sapply(1:z$nGroup, function(i)
					{
						tmp <- z$thetaMat[i, z$set1]
						use <- !is.na(tmp)
		# !is.na to get rid of parameters for other groups
						tmp[use] <- mvrnorm(1, mu=rep(0, sum(use)),
					Sigma=z$scaleFactors[i] *
					z$proposalCov[[i]][restrictProposal,restrictProposal] )
						tmp
					}
				))

		thetaOld <- z$thetaMat
		thetaOld[, z$basicRate] <- trafo(thetaOld[, z$basicRate])
		thetaNew <- thetaOld
		thetaNew[, varyingParameters] <- thetaOld[,varyingParameters] + thetaChanges
 		priorOld <- sapply(1:z$nGroup, function(i)
				{
					tmp <- thetaOld[i, z$set1]
					use <- !is.na(tmp)
					dmvnorm(tmp[use],  mean=z$muTemp[restrictProposal],
							sigma=z$SigmaTemp[restrictProposal,restrictProposal]) #density
				}
				)
		priorNew <- sapply(1:z$nGroup, function(i)
				{
					tmp <- thetaNew[i, z$set1]
					use <- !is.na(tmp)
					dmvnorm(tmp[use],  mean=z$muTemp[restrictProposal],
							sigma=z$SigmaTemp[restrictProposal,restrictProposal]) #density
				}
				)
		# use z still with z$thetaMat = thetaOld:
		logpOld <- getProbabilitiesFromC(z)[[1]]
		# note that z$loglik == sum(logpOld) # no I think it is logpOld itself
		if (z$incidentalBasicRates)
		{
			thetaNew[, z$basicRate] <- antitrafo(thetaOld[, z$basicRate])
		}
		else
		{
			thetaNew[, z$basicRate] <- antitrafo(thetaNew[, z$basicRate])
		}
		z$thetaMat <<- thetaNew
		logpNew <- getProbabilitiesFromC(z)[[1]]
		proposalProbability <- priorNew - priorOld + logpNew - logpOld
# cat(proposalProbability, "\n", priorNew, priorOld, "\n",
#									logpNew, "\n", logpOld, "\n")
		acceptProposals <- (log(runif(length(proposalProbability))) <
							proposalProbability)
		acceptProposals[is.na(acceptProposals)] <- FALSE # very rare but it did happen
		thetaOld[, z$basicRate] <- antitrafo(thetaOld[, z$basicRate])
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
			browser()
		}
		proposals.accept
	} # end sampleVaryingParameters


	##@sampleConstantParameters internal sienaBayes; propose new parameters
	## eta (the non-varying i.e. constant parameters), and accept them or not
	## Replace accepted new parameters only if change.
	## MCversion = 1 proposes change only for the eta parameters,
	## MCversion = 0 proposes one vector to change all non-rate parameters
	## (varying and nonvarying)
	sampleConstantParameters <- function(change=TRUE, MCversion)
	{
	# do the chamarrita
	# get a multivariate normal with covariance matrix
	# proposalC0 (for MCversion=0) or proposalC0eta (for MCversion = 1)
	# multiplied by a scale factor
#cat("SampleConstantParameters\n")
#browser()

		if (MCversion <= 0)
		{
			etaChange <- mvrnorm(1, mu=rep(0, sum(z$forEtaVersion0)),
								Sigma=z$scaleFactor0 * z$proposalC0)
		}
		else
		{
			etaChange <- mvrnorm(1, mu=rep(0, sum(z$set2)),
								Sigma=z$scaleFactorEta * z$proposalC0eta)
		}
		thetaOld <- z$thetaMat
		thetaNew <- thetaOld
		if (MCversion <= 0)
		{
			thetaNew[, z$forEtaVersion0] <-
				thetaOld[, z$forEtaVersion0] +
				matrix(etaChange, z$nGroup, sum(z$forEtaVersion0), byrow=TRUE)
		}
		else
		{
			thetaNew[, z$set2] <- thetaOld[, z$set2] +
							matrix(etaChange, z$nGroup, z$p2, byrow=TRUE)
		}
		# use z still with z$thetaMat = thetaOld:
		logpOld <- getProbabilitiesFromC(z)[[1]]
		# note that z$loglik == sum(logpOld) # no I think it is logpOld itself
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
			chol2inv(chol(rWishart(1,v,chol2inv(chol(S)))[,,1]))
		}

		muhat <- matrix( c(0), z$p1, 1)
		Thetas <- t(matrix( z$thetaMat[!is.na(z$thetaMat)],
									z$nGroup, z$TruNumPars))[z$varyingParametersInGroup,]
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
		if (z$incidentalBasicRates)
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
						z$nGroup, z$TruNumPars)[, restrictUpdate, drop=FALSE]
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
			z$SigmaTemp <<- z$SigmaTemp -
						bgain * correct * (z$SigmaTemp - matQ)
		# Now make this matrix positive definite, if it is not.
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
#				z$SigmaTemp[updated,updated] %*% diag(1/standdev, nrow=lsd)
			es <- eigen(z$SigmaTemp[updated,updated])
#			es <- eigen(corrMatrix)
			esv <- es$values
			cat("Eigenvalues Sigma = ", sort(esv), "\n")
			delta0 <- z$delta
			# Read the following maths as if delta = .Machine$double.eps = 0.
#			if (min(esv) <= .Machine$double.eps)
			if (min(esv) < delta0)
			{
				tol <- sum(updated)*max(abs(esv))*.Machine$double.eps
				delta <-  2*tol # factor two is just to make sure the resulting
				# matrix passes all numerical tests of positive definiteness
				# TS: I replace this delta by a larger number
				delta <- delta0
				tau <- pmax(0, delta - esv)
				cat("Smallest eigenvalue of Sigma now is ",
								min(esv),"; make posdef.\n")
				dm =  es$vectors %*%
						diag(tau, sum(updated)) %*% t(es$vectors)
				z$SigmaTemp[updated,updated] <<-
#					diag(standdev, nrow=lsd) %*% (corrMatrix + dm) %*%
#							diag(standdev, nrow=lsd)
				z$SigmaTemp[updated,updated] + dm
			}
		}
	}# end RobbinsMonro

	##@averageTheta internal sienaBayes; algorithm to average past theta values
	averageTheta <- function()
	{
#	cat('averageTheta\n')
#drop browser()
		thetaMean <- rep(NA, z$pp)
		for (group in 1:z$nGroup)
		{
			thetaMean[z$ratePositions[[group]]] <- 	colMeans(
				z$ThinParameters[, group, !z$generalParametersInGroup,
						drop=FALSE])
		}
		thetaMean[(z$set1)&(!z$basicRate)] <- colMeans(
				z$ThinPosteriorMu[,z$objectiveInVarying,
						drop=FALSE])
		if (any(z$set2))
		{
			thetaMean[z$set2] <- colMeans(z$ThinPosteriorEta[,, drop=FALSE])
		}
		thetaMean
	}

	##@somePlots make some plots during operation of the function
	# probably requires storeAll
	somePlots <- function()
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
			 thetadf[, -1, drop=FALSE][, !z$basicRate, drop=FALSE])
        thetaNames<- paste(z$effects$name[!z$basicRate],
                           z$effects$shortName[!z$basicRate], sep=".")
        rateNames <- paste(z$effects$name[basicRate],
			    z$effects$shortName[basicRate], z$effects$period[basicRate],
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
			varcall <- paste(varnames,  "~ 1:", ii * nrunMHBatches * z$nGroup,
                        " | Group", sep="", collapse="")
		}
		else
		{
			varcall <- paste(varnames,  "~ 1:", ii * nrunMHBatches * z$nGroup,
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

	##@savePartial makes a partial intermediate save as a sienaBayesFit object
	savePartial <- function(z)
	{# save intermediate results
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
		save(z,file="PartialBayesResult.RData")
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
	ctime1 <- proc.time()[1]
	ctime <- proc.time()[1]

	if (is.null(prevBayes))
	{
		z <- initializeBayes(data, effects, algo, nbrNodes,
						initgainGlobal=initgainGlobal,
						initgainGroupwise=initgainGroupwise,
						priorMu=priorMu, priorSigma=priorSigma,
						priorDf=priorDf, priorKappa=priorKappa,
						frequentist=frequentist,
						incidentalBasicRates=incidentalBasicRates,
						delta=delta, nmain=nmain, nwarm=nwarm,
						lengthPhase1=lengthPhase1, lengthPhase3=lengthPhase3,
						prevAns=prevAns, silentstart=silentstart,
						clusterType=clusterType)
		flush.console()
		createStores()
		z$sub <- 0
		z$nrunMHBatches <- nrunMHBatches
		class(z) <- "sienaBayesFit"

		# Determine multiplicative constants for proposal distribution
		ctime1 <- proc.time()[1]
		cat(ctime1-ctime,'\n')
		improveMH()
		ctime2<- proc.time()[1]
		cat('improveMH', ctime2-ctime1,'seconds.\n')
		flush.console()
	}

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
	if (!is.null(prevBayes))
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
			stop("Do not use this for prevBayes.")
		}
	}

# The division into phases 1-3 is relevant only for frequentist estimation.

	if (is.null(prevBayes))
	{
		endPhase1 <- nwarm + z$lengthPhase1
		endPhase2 <- nwarm + nmain - z$lengthPhase3
		z$frequentist <- frequentist

		# Warming phase
		bgain <- initfgain*exp(-gamma)
		for (ii in 1:nwarm)
		{
			cyc <- MCMCcycle(nrunMH=nrunMHBatches, bgain=bgain)
			z$iimain <- ii
			storeData(cyc)
			cat('Warming step',ii,'(',nwarm,')\n')
			cat("Accepts ",sum(cyc$BayesAcceptances),"/",
				z$nGroup*nrunMHBatches,"\n")
			flush.console()
		}
		print('end of warming')
		ctime3<- proc.time()[1]
		cat('warming took', ctime3-ctime2,'seconds.\n')
		flush.console()
		# Improve parameter value by averaging
	# Average the groupwise parameters over the last half of Phase 2
		halfPhaseWarm <- floor(nwarm/2 + 1)
		for (group in 1:z$nGroup)
		{
			z$thetaMat[group, z$ratePositions[[group]]] <-
				colMeans(z$ThinParameters[halfPhaseWarm:nwarm, group,
				!z$generalParametersInGroup, drop=FALSE], dims=1)
			z$thetaMat[group, !z$basicRate] <- colMeans(
				z$ThinParameters[halfPhaseWarm:nwarm, group,
				z$generalParametersInGroup, drop=FALSE], dims=1)
		}
		# Average the past groupwise parameters over runs and groups
		z$muTemp <-
			colMeans(z$ThinParameters[halfPhaseWarm:nwarm, ,
						z$varyingParametersInGroup, drop=FALSE], dims=2)
		# Average the covariance matrices over the runs
		tmp <- lapply(halfPhaseWarm:nwarm, function(i, x)
				{cov(x[i,,z$varyingParametersInGroup])},
										x=z$ThinParameters)
		z$SigmaTemp <- Reduce('+', tmp) / (nwarm - halfPhaseWarm + 1)
		cat('Parameter values after warming up\n')
		print(round(z$thetaMat,2))
		cat('\n')
		ctime1 <- proc.time()[1]
		cat('Second ')
		improveMH()
		ctime4<- proc.time()[1]
		cat('Second improveMH', ctime4-ctime3,'seconds.\n')
		flush.console()
	}
	else # (!is.null(prevBayes))
	{ # what remains to be done of the initialization
		z <- prevBayes
		nwarm <- 0
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
			ctime4<- proc.time()[1]
		}

		z$FRAN <- getFromNamespace(algo$FRANname, pkgname)
		z <- initializeFRAN(z, algo, data=data, effects=effects,
                prevAns=prevAns, initC=FALSE, onlyLoglik=TRUE)

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
			clusterCall(z$cl, initializeFRAN, z, algo,
                    initC = TRUE, profileData=FALSE, returnDeps=FALSE)
			clusterSetRNGStream(z$cl, iseed = as.integer(runif(1,
								max=.Machine$integer.max)))
		}
	}

	z$OK <- FALSE

	# Main iterations
	# Phase 1
	bgain <- initfgain*exp(-gamma)
    for (ii in (nwarm+1):endPhase1)
    {
		cyc <- MCMCcycle(nrunMH=nrunMHBatches, bgain=bgain)
		z$iimain <- ii
		storeData(cyc)
		cat('main', ii, '(', nwarm+nmain, ')')
		if (frequentist) cat(', Phase 1')
		cat("\nMu = ", round(z$muTemp,4), "\nEta = ",
					round(z$ThinPosteriorEta[ii,],4), "\nSigma = \n")
		print(round(z$SigmaTemp,4))
		cat("\n")
		cat("Accepts ",sum(cyc$BayesAcceptances),"/",
				z$nGroup*nrunMHBatches,"\n")
	cat("thetaMat = \n")
	print(round(z$thetaMat,2))
		cat("\n")
		flush.console()
#cat("cyc\n")
#browser()
		scores.ans.last <- t(cyc$ans.last$fra) # the scores a z$pp by z$nGroups matrix
		dfra.ans.last <- cyc$ans.last$dff # the Hessian a sparse z$pp by z$pp matrix
		if (frequentist)
		{
			RobbinsMonro(bgain)
		}
		# save intermediate results
		if (saveFreq >= 2)
		{
			if (ii %% saveFreq == 0)
			{
				savePartial(z)
			}
		}
        if (ii %% 10 == 0 && plotit) ## do some plots
        {
		somePlots()
        }
    }

	cat("thetaMat = \n")
	print(round(z$thetaMat,2))
	ctime5<- proc.time()[1]
	cat('Phase 1 took', ctime5-ctime4,'seconds.\n')
	flush.console()

#	cat('Third ')
#	improveMH()
	flush.console()

	# Phase 2
    for (ii in (endPhase1+1):endPhase2)
    {
		bgain <- initfgain*exp(-gamma*log(ii-endPhase1))
		cyc <- MCMCcycle(nrunMH=nrunMHBatches, bgain=bgain)
		z$iimain <- ii
		storeData(cyc)
		cat('main', ii, '(', nwarm+nmain, ')')
		if (frequentist) cat(', Phase 2')
		cat("\nMu = ", round(z$muTemp,4), "\nEta = ",
				round(z$ThinPosteriorEta[ii,],4), "\nSigma = \n")
		print(round(z$SigmaTemp,4))
		cat("\n")
		cat("Accepts ",sum(cyc$BayesAcceptances),"/",
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
				savePartial(z)
			}
		}
        if (ii %% 50 == 0 && plotit) ## do some plots
        {
		somePlots()
        }
    }

	# Process results Phase 2.
	halfPhase2 <- floor((endPhase1 + endPhase2)/2 + 1)
	if (frequentist)
	{
	# Average muTemp, eta (in thetaMat) and SigmaTemp
	# over the last half of Phase 2
		z$muTemp <- colMeans(z$ThinPosteriorMu[halfPhase2:endPhase2, ])
		z$thetaMat[1:z$nGroup,z$set2] <-
				colMeans(z$ThinPosteriorEta[halfPhase2:endPhase2,,drop=FALSE])
		z$SigmaTemp <- colMeans(z$ThinPosteriorSigma[halfPhase2:endPhase2,,],
												dims=1)
	}
	if (z$incidentalBasicRates)
	{
	# Average the rate parameters over the last half of Phase 2
		for (group in 1:z$nGroup)
		{
			z$thetaMat[group, z$ratePositions[[group]]] <-
				colMeans(z$ThinParameters[halfPhase2:endPhase2, group,
					!z$generalParametersInGroup, drop=FALSE], dims=1)
		}
	}

	ctime6<- proc.time()[1]
	cat('Phase 2 took', ctime6-ctime5,'seconds.\n')
	cat("\nMu = ", round(z$muTemp,4), "\nEta = ",
				round(z$ThinPosteriorEta[ii,],4), "\nSigma = \n")
	print(round(z$SigmaTemp,4))
	cat("\n")
	flush.console()

# Phase 3
	cat('Phase 3\n')
    for (ii in (endPhase2+1):(nwarm+nmain))
    {
		cyc <- MCMCcycle(nrunMH=nrunMHBatches, bgain=0.000001)
		z$iimain <- ii
		storeData(cyc)
		cat('main', ii, '(', nwarm+nmain, ')')
		if (frequentist) cat(', Phase 3')
		cat("\nMu = ", round(z$muTemp,4), "\nEta = ",
				round(z$ThinPosteriorEta[ii,],4), "\nSigma = \n")
		print(round(z$SigmaTemp,4))
		cat("\n")
#		flush.console()
		if (frequentist)
		{
			RobbinsMonro(0.0, ii-endPhase2, change=FALSE)
		}
		else
		{
			cat("Accepts ",sum(cyc$BayesAcceptances),"/",
				z$nGroup*nrunMHBatches,"\n")
		}
	# save intermediate results
		if (saveFreq >= 2)
		{
			if (ii %% saveFreq == 0)
			{
				savePartial(z)
			}
		}
        if (ii %% 10 == 0 && plotit) ## do some plots
        {
			somePlots()
        }
    }

# Process results of Phase 3
	ctime7<- proc.time()[1]
	cat('Total duration',ctime7-ctime,'seconds.\n')
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
    class(z) <- "sienaBayesFit"
	if (!storeAll)
	{
		z$acceptances <- NULL
		z$candidates <- NULL
	}
	z$OK <- TRUE
    z
} # end sienaBayes

##@initializeBayes algorithms do set up for Bayesian model
initializeBayes <- function(data, effects, algo, nbrNodes,
						initgainGlobal, initgainGroupwise,
						priorMu, priorSigma, priorDf, priorKappa,
						frequentist, incidentalBasicRates, delta,
						nmain, nwarm,
						lengthPhase1, lengthPhase3,
						prevAns, silentstart, clusterType=c("PSOCK", "FORK"))
{
	##@precision internal initializeBayes invert z$covtheta
	## avoiding some inversion problems
	## for MoM estimates only
	precision <- function(z)
	{
		require(MASS)
#		t(z$dfrac) %*% ginv(z$msfc) %*% z$dfrac
# The ...c versions are corrected and have identity rows/columns
# for fixed parameters. This correction is not to be employed here.
		t(z$dfra) %*% ginv(z$msf) %*% z$dfra
	}#end precision in initializeBayes

	##@runCumSum internal initializeBayes used for projectEffects only
	runCumSum <- function(x)
	{
	# A function that makes cumulative sums of uninterrupted runs
	# for logical vectors
	# e.g., FTTFFFTFTTTTF produces 0120001012340
		cumsum.inrun <- rep(0,length(x))
		if (x[1]) {cumsum.inrun[1] <- 1}
		for (i in 2:length(x))
		{
			if (x[i]){cumsum.inrun[i] <- cumsum.inrun[i-1]+1}
		}
		cumsum.inrun
	}

	##@projectEffects internal initializeBayes project model specification
	# modified version of updateTheta
	# model specification of prevAns is transferred to effs for group ii
	# and !randomEffects are fixed.
	projectEffects <- function(effs, prevAns, ii)
	{
		prevEffects <- prevAns$requestedEffects
		prevEffects$initialValue <- prevAns$theta
		# Some tricks need to be done to get matches on interactions right
		oldInteract <- (prevEffects$shortName %in% c("unspInt" , "behUnspInt"))
		newInteract <- (effs$shortName %in% c("unspInt", "behUnspInt"))
		oldInterNumber <- runCumSum(oldInteract)
		newInterNumber <- (runCumSum(newInteract) + 2) %/% 3
		# the + 2) %/% 3 takes care of the eval/endow/creation multiplicities
		oldlist <- apply(prevEffects, 1, function(x)
						paste(x[c("name", "shortName",
								"type",
								"interaction1", "interaction2",
								"period")],
							collapse="|"))
		oldlist <- paste(oldlist,oldInterNumber)
		efflist <- apply(effs, 1, function(x)
						paste(x[c("name", "shortName",
								"type",
								"interaction1", "interaction2",
								"period")],
							collapse="|"))
		efflist <- paste(efflist,newInterNumber)
		use <- efflist %in% oldlist
		effs$include[use] <-
			prevEffects$include[match(efflist, oldlist)][use]
		effs$initialValue[use] <-
			prevEffects$initialValue[match(efflist, oldlist)][use]
		effs$initialValue[effs$basicRate] <-
			prevEffects$initialValue[
				prevEffects$basicRate & (prevEffects$group == ii)]
		effs$fix[effs$include] <-
			(!prevEffects$randomEffects[na.omit(match(efflist,oldlist))])|
			(prevEffects$fix[na.omit(match(efflist,oldlist))])
# In the full effects object used for prevAns,
# the basic rate parameters are duplicated for all groups,
# whereas in the effs objcet they are there only for
# a single group.
# This means that the effect numbers are different;
# since effect numbers are used to define the interactions
# in effect1, effect2, effect3,
# the different sequence numbers must be taken into account.
# The vector of {differences between sequence numbers
# in the effects object for prevAns and effs} of
# objective function effects is:
# number of basic rate effects per dep. var. =
#    sum(effs$shortName=="Rate")/length(unique(effs$name))
# * (number of groups - 1)
# * sequence number in effs of dependent variable =
#    match(effs$name,unique(effs$name))
		depVarNumber <- sapply(unique(effs$name),
						function(a){match(a, unique(effs$name))})
		extraRates <-  sum(effs$shortName=="Rate") *
				(max(prevAns$effects$group)-1) /
				length(unique(effs$name))
		position.Difference <- as.vector(extraRates * depVarNumber[effs$name])
		effs$effect1[use] <- pmax((prevEffects[prevAns$effects$group==1,])$effect1 -
								position.Difference[use], 0)
		effs$effect2[use] <- pmax((prevEffects[prevAns$effects$group==1,])$effect2 -
								position.Difference[use], 0)
		effs$effect3[use] <- pmax((prevEffects[prevAns$effects$group==1,])$effect3 -
								position.Difference[use], 0)
		effs
	}#end projectEffects in initializeBayes

	## start initializeBayes
    ## initialize
	if (!inherits(data,"siena")){stop("The data is not a valid siena object.")}
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

	# Short specification for initial model.
	# Use diagonalize=1 to increase stability.
	if (initgainGlobal <= 0)
	{
		startupModel <- sienaAlgorithmCreate(n3=500, nsub=0, cond=FALSE,
						"startup_sienaBayes")
	}
	else
	{
		startupModel <- sienaAlgorithmCreate(n3=500, nsub=3, cond=FALSE,
							diagonalize=1.0, firstg=initgainGlobal,
							projname="startup_sienaBayes")
	}
	cat("Estimate initial global parameters\n")
	startupGlobal <- siena07(startupModel, data=data, effects=effects,
								batch=TRUE, silent=silentstart)
	z$startupGlobal <- startupGlobal
	cat("Initial global estimates\n")
	z$initialResults <- print(startupGlobal)
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
				"\nAdvice: use a smaller value of initgainGlobal.\n\n")
		stop("Divergent initial estimate.")
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
	npars <- sum(startupGlobal$effects$group==1)
	# number of Dependent variables:
	nDepvar <- length(unique(startupGlobal$effects$name))
	if (ngroups >= 2)
	{
	#	if ((startupGlobal$pp - npars)%%(ngroups-1) != 0)
		if (var(sapply(data,function(x){x$observations})) > 1e-6)
		{
			stop("Not all groups have the same number of periods.")
		}
	}
	nBasicRatesPerGroup <- sum(startupGlobal$effects$basicRate) / ngroups
	rateParameters <- matrix(NA, nBasicRatesPerGroup, ngroups)
	nActors <- rep(0, ngroups)
	initialEstimates <- matrix(0, ngroups, startupGlobal$pp)
	factorsBasicRate <- matrix(0, ngroups, startupGlobal$pp)
	if (initgainGroupwise <= 0)
	{
		startup1Model <- sienaAlgorithmCreate(n3=500, nsub=0, cond=FALSE)
	}
	else
	{
		startup1Model <- sienaAlgorithmCreate(n3=500, nsub=1,
										firstg=initgainGroupwise, cond=FALSE)
	}
	# Define the partition for the varying non-fixed (set1)
	# and non-varying non-fixed (set2) effects;
	# effects that are not estimated (fix) are
	# excluded from both sets.
	z$set1 <- rep(FALSE, startupGlobal$pp)
	z$set1[startupGlobal$effects$basicRate] <- TRUE
	z$set1[startupGlobal$effects$randomEffects &
				(!startupGlobal$effects$fix)] <- TRUE
	z$set2 <- !z$set1
	z$set2[startupGlobal$effects$fix] <- FALSE

	# Number of non-varying non-fixed parameters
	# among the z$truNumPars true parameters:
	z$p2 <- sum(z$set2)
	# Number of varying parameters among the z$truNumPars true parameters:
	z$p1 <- npars - z$p2 - sum(startupGlobal$effects$fix)

	# compute covariance matrices for random walk proposals for all groups:
	# proposalC0 for all effects;
	# proposalC0eta for eta, i.e., non-varying effects.
	# proposalCov for theta^{(1)}_j, i.e., groupwise effects.
	z$proposalCov <- list()
	prec <- precision(startupGlobal)
	rates <- startupGlobal$effects$basicRate
	precRate <- sum(diag(prec[rates,rates]))
	prec[rates,rates] <- 0
	diag(prec)[rates] <- precRate
	if (inherits(try(z$proposalC0 <- chol2inv(chol(prec)),
					silent=TRUE), "try-error"))
	{
		cat("precision(startupGlobal) non-invertible.\n")
		cat("This probably means the model is not identifiable.\n")
		stop("Non-invertible precision(startupGlobal).")
	}
	z$proposalC0eta <- z$proposalC0[z$set2, z$set2, drop=FALSE]
	z$forEtaVersion0 <- (z$set1|z$set2)&(!(rates|startupGlobal$effects$fix))
	z$proposalC0 <- z$proposalC0[z$forEtaVersion0, z$forEtaVersion0, drop=FALSE]

	for (i in 1:ngroups)
	{
		cat(paste("Group",i,"\n"))
		flush.console()
	# roughly estimate varying parameters
		nActors[i] <- length(data[[i]]$nodeSets[[1]])
		# "use" will indicate the coordinates of the global parameters
		# used for this group,
		# and "rates" indicates the rate parameters among these.
		use <- rep(TRUE, startupGlobal$pp)
		rates <- startupGlobal$effects$basicRate
		use[rates] <- FALSE
		use[startupGlobal$effects$group==i] <- TRUE
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
			effects0 <- projectEffects(effects0,startupGlobal,i)
			# projectEffects also fixes the !randomEffects effects.
			# get the parameter values from the global fit:
			effects0$initialValue[effects0$include] <-
				startupGlobal$theta[use]
			cat("Estimate initial parameters group", i, "\n")
			startup1 <- siena07(startup1Model, data=data[[i]],
					effects=effects0, batch=TRUE, silent=silentstart)
			cat("\nInitial estimate obtained\n", startup1$theta,"\n")
		}
		initialEstimates[i,use] <- startup1$theta
		factorsBasicRate[i,use] <- (1/nActors[i])*startup1$theta
		# non-rate parameters will be dropped from this object, hence the name
		rateParameters[,i] <- trafo(startup1$theta[rates])
		# estimate uncertainties for this group:
		covmati <- chol2inv(chol(
			w*precision(startup1) +
			((1-w)/ngroups)*precision(startupGlobal)[use,use,drop=FALSE]))
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

	if (max(abs(initialEstimates), na.rm=TRUE) > 400)
	{
		cat("\nThe initial estimates are\n")
		print(initialEstimates)
		cat("\nThe largest absolute value is ",
				max(abs(initialEstimates), na.rm=TRUE),
				"\nwhich is too large for good operation of this function.\n",
				"\nAdvice: use a smaller value of initgainGroupwise.\n\n")
		stop("Divergent initial estimate.")
	}
	z$effectName <- effects0$effectName[effects0$include]
	z$requestedEffects <- startup1$requestedEffects

   	z$FRAN <- getFromNamespace(algo$FRANname, pkgname)
    z <- initializeFRAN(z, algo, data=data, effects=effects,
                prevAns=prevAns, initC=FALSE, onlyLoglik=TRUE)
	z$thetaMat <- initialEstimates
	z$basicRate <- z$effects$basicRate
    z$nGroup <- z$f$nGroup
	is.batch(TRUE)

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
        clusterCall(z$cl, initializeFRAN, z, algo,
                    initC = TRUE, profileData=FALSE, returnDeps=FALSE)
		clusterSetRNGStream(z$cl, iseed = as.integer(runif(1,
								max=.Machine$integer.max)))
    }


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
	# Indicator of parameters in the parameter vector of length pp,
	# for which the hierarchical model applies,
	# which is the set of all parameters without the basic rate parameters:
#	z$generalParameters <- (z$effects$group == 1) & (!z$basicRate)
	# Dependent variable:
#	z$dependentVariable <- z$effects$name


	# Number of parameters in each group:
	z$TruNumPars <- sum( !z$basicRate ) + length(z$ratePositions[[1]])
	scf <- 1
	if (z$TruNumPars > 5)
	{
		scf <- 5/z$TruNumPars # this is just a try
	}
    z$scaleFactors <- rep(scf, z$nGroup)
	z$scaleFactor0 <- scf
	z$scaleFactorEta <- scf

	# Construction of indicator of parameters of the hierarchical model
	# within the groupwise parameters:

	vec1 <- 4 - effects$randomEffects[effects$include]
	vec1[z$basicRate] <- 2
	vec1[z$ratePositions[[1]]] <- 1
	vec1[effects$fix[effects$include]] <- 5
	# now 1=rates group 1; 2=other rates; 3=randomly varying;
	# 4 = estimated non-varying (eta); 5 = non-estimated, fixed.
	# set1 = (1,2,3) set2 = (4)
	z$generalParametersInGroup <- vec1[vec1!= 2] %in% c(3,4,5)
	z$varyingGeneralParametersInGroup <- vec1[vec1!= 2] == 3
	z$varyingParametersInGroup <- vec1[vec1!= 2] %in% c(1,3)
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
#				chol2inv(chol(precision(startupGlobal)
#					[z$generalParameters,z$generalParameters]))
#			z$dmuinv <- diag(1, z$TruNumPars, z$TruNumPars)
#		}
	}
	# Determination parameters of prior distribution;
	# this specifies the defaults!
	if (is.null(priorMu))
	{
		z$priorMu <- matrix( c(0), z$p1, 1 )
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
	# Rate parameters are nuisance parameters, and get a data-induced prior;
	# observed covariance matrix on the trafo scale plus a small ridge.

	z$priorMu[z$ratesInVarying] <- rowMeans(rateParameters)
	z$priorSigma[z$ratesInVarying,] <- 0
	z$priorSigma[,z$ratesInVarying] <- 0
	z$priorSigma[z$ratesInVarying,z$ratesInVarying] <- cov(t(rateParameters))
	diag(z$priorSigma)[z$ratesInVarying] <-
				diag(z$priorSigma)[z$ratesInVarying] + 0.1
	z$incidentalBasicRates <- incidentalBasicRates
	z$delta <- delta
	if (nmain < 10){nmain <<- 10}

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
#browser()
	meanThetaMat <- colMeans(matrix( z$thetaMat[!is.na(z$thetaMat)],
									z$nGroup, z$TruNumPars))
	z$muTemp <- matrix(meanThetaMat[z$varyingParametersInGroup], z$p1, 1)
	z$SigmaTemp <- diag(z$p1)
    z
} # end initializeBayes


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

##@dmvnorm algorithms calculates multivariate normal density
##inefficient: should not call mahalanobis and eigen with same sigma repeatedly
dmvnorm <- function(x, mean , sigma)
{
    if (is.vector(x))
    {
        x <- matrix(x, ncol=length(x))
    }
	if (min(eigen(sigma)$values)< 1e-8)
	{
		cat('singular covariance matrix in dmvnorm\n')
		browser()
	}
    distval <- mahalanobis(x, center=mean, cov=sigma)
    logdet <- sum(log(eigen(sigma, symmetric=TRUE, only.values=TRUE)$values))
    -(ncol(x) * log(2 * pi) + logdet + distval) / 2
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
    .Call("getChainProbabilities", PACKAGE = pkgname, f$pData,
		  f$pModel, as.integer(x[1]), as.integer(x[2]),
		  as.integer(index), f$myeffects, theta, getScores)

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

