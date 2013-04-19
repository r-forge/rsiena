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
sienaBayes <- function(data, effects, model,saveFreq=100,
				initgainGlobal=0.1, initgainGroupwise = 0.02, initfgain=0.01,
				priorMu=NULL, priorSigma=NULL, priorDf=NULL, priorKappa=NULL,
#				fixedZeros=NULL,
				frequentist=FALSE, incidentalBasicRates=FALSE,
				gamma=0.5, delta=1e-4,
				nwarm=50, nmain=250, nrunMHBatches=20,
				lengthPhase1=round(nmain/5), lengthPhase3=round(nmain/5),
				storeAll = FALSE, prevAns=NULL,
				plotit=FALSE, nbrNodes=1, clusterType=c("PSOCK", "FORK"),
				getDocumentation=FALSE)
{
    ##@createStores internal sienaBayes Bayesian set up stores
    createStores <- function()
    {
#        npar <- length(z$theta)

		numberRows <- nmain * nrunMHBatches
#        z$posteriorTot <<- matrix(0, nrow=z$nGroup, ncol=npar)
#        z$posteriorMII <<- array(0, dim=c(z$nGroup, npar, npar))
#        z$candidates <<- array(NA, dim=c(numberRows, z$nGroup, npar))
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
		# tom's additional thinned stores
		z$ThinPosteriorMu <<- matrix(NA, nrow=nmain, ncol=z$TruNumPars)
		z$ThinPosteriorSigma <<-
					array( NA, dim=c(nmain, z$TruNumPars, z$TruNumPars) )
		z$ThinLoglik <<- rep(NA, nmain)
		z$ThinParameters <<- array(NA,
						dim=c(max(nwarm, nmain), z$nGroup, z$TruNumPars))
		z$ThinBayesAcceptances <<- matrix(NA, nrow=nmain, ncol=z$nGroup)
		z$sumBasicRates <<- matrix(0, nrow=z$nGroup, ncol=sum(z$basicRate))
		z$ThinScores <<- matrix(NA, nrow=nmain, ncol=sum(z$basicRate))
	}

    ##@storeData internal sienaBayes; put data in stores
    storeData <- function()
    {
        start <- z$sub + 1
#        nrun <- nrow(z$parameters)
        nrun <- ncol(z$accepts)
        end <- start + nrun - 1
#        z$posteriorTot <<- z$posteriorTot + colSums(z$parameters)
		if (storeAll)
		{
			z$acceptances[start:end, ] <<- z$accepts
#			z$candidates[start:end,, ] <<- z$parameters
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
		if (z$incidentalBasicRates)
		{
			z$ThinScores[z$iimain, ] <<- z$thescores
		}
        for (group in 1:z$nGroup)
        {
			# Drop those elements of the parameters array
			# that are NA because for a given group they refer
			# to rate parameters of other groups.
			z$ThinParameters[z$iimain, group, ] <<-
					z$thetaMat[group, !is.na(z$thetaMat[group, ])]
				#z$parameters[nrun, group,!is.na(z$parameters[nrun, group, ])]
#            for (i in dim(z$parameters)[1])
#            {
#                z$posteriorMII[group, , ] <<- z$posteriorMII[group, ,] +
#                    outer(z$parameters[i, group, ], z$parameters[i, group, ])
#            }
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

	##@improveMH internal sienaBayes; find scale factors
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
            success <<- (number >= 2)
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
            MCMCcycle(nrunMH=1, nrunMHBatches=100, change=FALSE, bgain=-1)
            actual <- z$BayesAcceptances
			## z$BayesAcceptances is a vector of the number of
			## acceptances by group in the MH change of parameters
			## out of nrunMHBatches trials
# Julie Kratz got a bug here with non-permitted NAs in
# the update of z$scaleFactors.
# I (TS) do not know how this could arise, and inserted some statements
# to exclude NAs.
			actual[is.na(actual)] <- 0
            ans <- rescaleCGD(100)
            update <- (number < 3)
			update[is.na(update)] <- FALSE # this should not be necessary
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
		z$allScaleFactors[[length(z$allScaleFactors)+1]] <<-
				z$scaleFactors
        cat('fine tuning took ', iter, ' iterations.\n')
		cat('Scalefactor:', z$scaleFactors, '\n')
		flush.console()
    }
	##@MCMCcycle internal sienaBayes; some loops of (MH steps and sample parameters)
	MCMCcycle <- function(nrunMH, nrunMHBatches, change=TRUE, bgain)
	{
		z$accepts <<- matrix(NA, nrow=z$nGroup, nrunMHBatches)
#		z$parameters <<- array(NA, dim=c(nrunMHBatches, z$nGroup, z$pp))
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
		# This is the hot loop.
#			if (z$someFixedZeros)
#			{
#				oldTheta <- z$thetaMat
#				z$thetaMat[z$fixedZeros] <- 0
#			}
		# sample MH steps in the network chain:
			ans <- z$FRAN(z, byGroup=TRUE, returnLoglik=TRUE, onlyLoglik=TRUE)
										#	c2 <- proc.time()[1]
										#	cat ('fran',c2-cc,'\n')
			z$loglik <<- ans$loglik
			# The loglik is not directly used;
			# rather we need after the nrunMHBatches the loglik by group.
			# This is then used as logpOld, see below.
			# Perhaps efficiency can be gained here
			# by not calculating loglik all the time.
			# sample the group-level parameters:
			sampleParameters(change, bgain)
			# update grand means muTemp and covariance SigmaTemp
			# do we need to check that there are more than one group?
			if ( z$nGroup > 1)
			{
				if (!z$frequentist)
				{
				# sample the global parameters:
				sampleMu()
				}
				z$posteriorMu[i, ] <<- z$muTemp
				z$PosteriorSigma[i, ,] <<- z$SigmaTemp
			}
			z$loglikVec[i] <<- z$loglik
#			z$thetaMat[z$fixedZeros] <- oldThetaMat[z$fixedZeros]
			# back to Ruth
			z$accepts[, i] <<- z$accept
#			z$parameters[i, , ] <<- z$thetaMat
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

	##@sampleParameters internal sienaBayes; propose new parameters
	## and accept them or not
	## Replace accepted new parameters only if change.
	## If z$incidentalBasicRates, for the basic rate parameters
	## not a Bayesian MH step but a Robbins Monro step is made.
	sampleParameters <- function(change=TRUE, bgain)
	{
		require(MASS)

		if (z$incidentalBasicRates)
		{
			if ((change)&&(bgain > 0))
			{
		# Robbins-Monro step for basic rate parameters
		# Basic rate parameters truncated to minimum value of 0.1.
		# Smaller gain parameter because the groupwise rate parameters
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
			unmask1 <- !z$basicRate
			unmask2 <- z$generalParametersInGroup
		# These unmasks indicate the vector of the non-(basic rate parameters);
		# unmask1 in the long vector of z$pp effects,
		# unmask2 in the short vector of z$TruNumPars effects.
		}
		else
		{
			unmask1 <- rep(TRUE, z$pp)
			unmask2 <- rep(TRUE, z$TruNumPars)
		}

	 	# do the fandango
		## get a multivariate normal with covariance matrix proposalCov
		## multiplied by a scale factor which varies between groups
		thetaChanges <- t(sapply(1:z$nGroup, function(i)
			{
				tmp <- z$thetaMat[i,unmask1]
				use <- !is.na(tmp)
				tmp[use] <-
					mvrnorm(1, mu=rep(0, sum(use)),
							Sigma=z$scaleFactors[i] *
							z$proposalCov[[i]][unmask2,unmask2])
				tmp
			}
			))
		thetaOld <- z$thetaMat
		thetaOld[, z$basicRate] <- trafo(thetaOld[, z$basicRate])
		thetaNew <- thetaOld
		thetaNew[, unmask1] <- thetaOld[,unmask1] + thetaChanges
 		priorOld <- sapply(1:z$nGroup, function(i)
				{
					tmp <- thetaOld[i, unmask1]
					use <- !is.na(tmp)
					dmvnorm(tmp[use],  mean=z$muTemp[unmask2],
							sigma=z$SigmaTemp[unmask2, unmask2])
				}
				)
		priorNew <- sapply(1:z$nGroup, function(i)
				{
					tmp <- thetaNew[i, unmask1]
					use <- !is.na(tmp)
					dmvnorm(tmp[use],  mean=z$muTemp[unmask2],
							sigma=z$SigmaTemp[unmask2,unmask2])
				}
				)
		logpOld <- getProbabilitiesFromC(z)[[1]]
		# note that z$loglik == sum(logpOld)
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
		accepts <- log(runif(length(proposalProbability))) <
							proposalProbability
		z$accept <<- accepts
		thetaOld[, z$basicRate] <- antitrafo(thetaOld[, z$basicRate])
		if (!change)
		{
			z$thetaMat <<- thetaOld
		}
		else
		{
			z$thetaMat[!accepts, ] <<- thetaOld[!accepts, ]
		}
	}

	# johans new function for drawing thetas (top-level parameters)
	# tom changed riwish to rWishart - available from R 2.15.0
	##@sampleMu internal sienaBayes; Gibbs algorithm to sample new Mu parameter
	sampleMu <- function()
		{
		invWish <- function(v,S){
		# Draw from the inverse Wishart distribution with df v
		# and scale matrix S
		# inefficient if drawn multiply for the same S
			chol2inv(chol(rWishart(1,v,chol2inv(chol(S)))[,,1]))
		}
		Thetas <- matrix( c(0), z$TruNumPars, z$nGroup)
		muhat <- matrix( c(0), z$TruNumPars, 1)
		for ( groupiterator in c(1:z$nGroup) )
		{
			tmp <- z$thetaMat[groupiterator, ]
					use <- !is.na(tmp)
					Thetas[,groupiterator] <- tmp[use]
		}
		muhat <- rowMeans( Thetas )# the average across g
		matQ <- (z$nGroup - 1)*cov(t(Thetas))
# Tom: riwish replaced by invWish
# prior Lambda = z$priorDf*z$priorSigma
		z$SigmaTemp <<- invWish(z$priorDf + z$nGroup,
				z$priorDf*z$priorSigma + matQ +
				(z$priorKappa*z$nGroup/(z$priorKappa + z$nGroup))*
				tcrossprod( muhat - z$priorMu, muhat - z$priorMu ) )
		z$muTemp <<- t(chol( ( z$priorKappa + z$nGroup )^(-1)*z$SigmaTemp )) %*%
			rnorm( z$TruNumPars , 0 , 1 ) +
			(z$nGroup/( z$priorKappa + z$nGroup ))*muhat +
			(z$priorKappa/( z$priorKappa + z$nGroup ))*z$priorMu
	}

	##@RobbinsMonro internal sienaBayes; Robbins Monro step for Mu and Sigma
	RobbinsMonro <- function(bgain)
		{
		if (z$nGroup <= 1)
		{
			stop("Use frequentist method only for more than 1 group.")
		}
		Thetas <- matrix( z$thetaMat[!is.na(z$thetaMat)],
									z$nGroup, z$TruNumPars)
		muhat <- colMeans(Thetas)
		matQ <- ((z$nGroup - 1)/z$nGroup)*cov(Thetas)

		if (z$incidentalBasicRates)
		{
		# Basic rate parameters must be masked:
			unmask <- z$generalParametersInGroup
		}
		else
		{
			unmask <- rep(TRUE, z$TruNumPars)
		}
		correct <- 1
		maxdif <- max(abs(z$muTemp[unmask] - muhat[unmask]))
		# Truncate if necessary; the threshold is pretty arbitrary of course.
		# It would be better
		# to have truncation dependent on provisional standard errors.
		# Perhaps use standard deviations in phase 1?
		if (bgain*maxdif > 0.3)
		{
			correct <- 0.3/(bgain*maxdif)
			cat("Correction for jump in mu by ",correct,"\n")
		}
		# The update does not need to be applied only to the unmasked part,
		# because the masked part plays no role.
		z$muTemp <<-  z$muTemp - bgain * correct * (z$muTemp - muhat)
# or with a matrix:
#	z$muTemp - bgain * as.vector(z$dmuinv %*% (z$muTemp - muhat))
#		standdev <- sqrt(diag(z$SigmaTemp[unmask,unmask, drop=FALSE]))
#		lsd <- length(standdev)
		maxdif <- max(abs(z$SigmaTemp[unmask,unmask] - matQ[unmask,unmask]))
		if (bgain*maxdif > 0.2)
		{
			correct <- 0.2/(bgain*maxdif)
			cat("Correction for jump in Sigma by ",correct,"\n")
		}
		z$SigmaTemp <<- z$SigmaTemp - bgain * correct * (z$SigmaTemp - matQ)
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
#		corrMatrix <- diag(1/standdev, nrow=lsd) %*%
#			z$SigmaTemp[unmask,unmask] %*% diag(1/standdev, nrow=lsd)
		es <- eigen(z$SigmaTemp[unmask,unmask])
#		es <- eigen(corrMatrix)
		esv <- es$values
		cat("Eigenvalues Sigma = ", sort(esv), "\n")
		delta0 <- z$delta
		# Read the following maths as if delta = .Machine$double.eps = 0.
#		if (min(esv) <= .Machine$double.eps)
		if (min(esv) < delta0)
		{
			tol <- z$TruNumPars*max(abs(esv))*.Machine$double.eps
			delta <-  2*tol # factor two is just to make sure the resulting
				# matrix passes all numerical tests of positive definiteness
			# TS: I replace this delta by a larger number
			delta <- delta0
			tau <- pmax(0, delta - esv)
			cat("Smallest eigenvalue of Sigma now is ",
							min(esv),"; make posdef.\n")
			dm =  es$vectors %*%
					diag(tau, sum(unmask)) %*% t(es$vectors)
			z$SigmaTemp[unmask,unmask] <<-
#					diag(standdev, nrow=lsd) %*% (corrMatrix + dm) %*%
#							diag(standdev, nrow=lsd)
			z$SigmaTemp[unmask,unmask] + dm
		}
}

	##@averageTheta internal sienaBayes; algorithm to average past theta values
	averageTheta <- function()
	{
	thetaMean <- rep(NA, z$pp)
	for (group in 1:z$nGroup)
	{
		thetaMean[z$ratePositions[[group]]] <- 	colMeans(
			z$ThinParameters[, group, !z$generalParametersInGroup, drop=FALSE])
	}
	thetaMean[!z$basicRate] <- colMeans(
			z$ThinPosteriorMu[,z$generalParametersInGroup, drop=FALSE])
	thetaMean
	}

	##@somePlots make some plots during operation of the function
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

	z <- initializeBayes(data, effects, model, nbrNodes,
						initgainGlobal=initgainGlobal,
						initgainGroupwise=initgainGroupwise,
						priorMu=priorMu, priorSigma=priorSigma,
						priorDf=priorDf, priorKappa=priorKappa,
#						fixedZeros=fixedZeros,
						frequentist=frequentist,
						incidentalBasicRates=incidentalBasicRates,
						delta=delta, nmain=nmain, nwarm=nwarm,
						lengthPhase1=lengthPhase1, lengthPhase3=lengthPhase3,
						prevAns=prevAns, clusterType=clusterType)
    createStores()
    z$sub <- 0

	# Determine multiplicative constants for proposal distribution
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

	endPhase1 <- z$nwarm + z$lengthPhase1
	endPhase2 <- nmain - z$lengthPhase3

	# Warming phase
	bgain <- initfgain*exp(-gamma)
    for (ii in 1:nwarm)
    {
		MCMCcycle(nrunMH=z$nrunMH, nrunMHBatches=nrunMHBatches, bgain=bgain)
		z$iimain <- ii
		storeData()
		cat('Warming step',ii,'(',nwarm,')\n')
		cat("Accepts ",sum(z$BayesAcceptances),"/",
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
		colMeans(z$ThinParameters[halfPhaseWarm:nwarm, , , drop=FALSE], dims=2)
	# Average the covariance matrices over the runs
	tmp <- lapply(halfPhaseWarm:nwarm, function(i, x)cov(x[i,,]),
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

	# Main iterations
	# Phase 1
	bgain <- initfgain*exp(-gamma)
    for (ii in (nwarm+1):endPhase1)
    {
		MCMCcycle(nrunMH=z$nrunMH, nrunMHBatches=nrunMHBatches, bgain=bgain)
		z$iimain <- ii
		storeData()
		cat('main', ii, '(', nmain, ')')
		if (z$frequentist) cat(', Phase 1')
		cat("\nMu = ", z$muTemp, "\nSigma = \n")
		cat("Accepts ",sum(z$BayesAcceptances),"/",
				z$nGroup*nrunMHBatches,"\n")
#		print(round(z$SigmaTemp,4))
#		cat("thetaMat = \n")
#		print(round(z$thetaMat,4))
		cat("\n")
		flush.console()
		if (((ii - nwarm) %% 200) == 0)
		{
		# Determine step sizes again, the size of good steps may have changed
			improveMH()
		}
		if (z$frequentist)
		{
			RobbinsMonro(bgain)
		}
		# save intermediate results
		if (saveFreq >= 2)
		{
			if (ii %% saveFreq == 0)
			{
			save(z,file="PartialBayesResult.RData")
			}
		}
        if (ii %% 10 == 0 && plotit) ## do some plots
        {
		somePlots()
        }
    }
	ctime5<- proc.time()[1]
	cat('Phase 1 took', ctime5-ctime4,'seconds.\n')
	flush.console()

	cat('Third ')
	improveMH()
	flush.console()

	# Phase 2
    for (ii in (endPhase1+1):endPhase2)
    {
		bgain <- initfgain*exp(-gamma*log(ii-endPhase1))
		MCMCcycle(nrunMH=z$nrunMH, nrunMHBatches=nrunMHBatches, bgain=bgain)
		z$iimain <- ii
		storeData()
		cat('main', ii, '(', nmain, ')')
		if (z$frequentist) cat(', Phase 2')
		cat("\nMu = ", z$muTemp, "\nSigma = \n")
		cat("Accepts ",sum(z$BayesAcceptances),"/",
				z$nGroup*nrunMHBatches,"\n")
#		print(round(z$SigmaTemp,4))
#		cat("thetaMat = \n")
#		print(round(z$thetaMat,4))
#		cat("\n")
		flush.console()
		if ((ii %% 200) == 0)
		{
		# Determine step sizes again, the size of good steps may have changed
#			improveMH()
		}
		if (z$frequentist)
		{
			RobbinsMonro(bgain)
		}
		# save intermediate results
		if (saveFreq >= 2)
		{
			if (ii %% saveFreq == 0)
			{
			save(z,file="PartialBayesResult.RData")
			}
		}
        if (ii %% 10 == 0 && plotit) ## do some plots
        {
		somePlots()
        }
    }

# Process results Phase 2.
	halfPhase2 <- floor((endPhase1 + endPhase2)/2 + 1)
	if (z$frequentist)
	{
	# Average muTemp and SigmaTemp over the last half of Phase 2
		z$muTemp <- colMeans(z$ThinPosteriorMu[halfPhase2:endPhase2, ])
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
	cat("\nMu = ", z$muTemp, "\nSigma = \n")
	print(round(z$SigmaTemp,4))
	cat("\n")
	flush.console()

# Phase 3
	cat('Phase 3\n')
    for (ii in (endPhase2+1):nmain)
    {
		MCMCcycle(nrunMH=z$nrunMH, nrunMHBatches=nrunMHBatches, bgain=0.000001)
		z$iimain <- ii
		storeData()
		cat('main', ii, '(', nmain, ')')
		if (z$frequentist) cat(', Phase 3')
		cat("\nMu = ", z$muTemp, "\nSigma = \n")
		cat("Accepts ",sum(z$BayesAcceptances),"/",
				z$nGroup*nrunMHBatches,"\n")
#		flush.console()
	# save intermediate results
		if (saveFreq >= 2)
		{
			if (ii %% saveFreq == 0)
			{
			save(z,file="PartialBayesResult.RData")
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
	if (z$frequentist)
	{
		z$theta <- z$muTemp
		z$Sigma <- z$SigmaTemp
	}
	else
	{
		z$theta <- averageTheta()
	}
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
initializeBayes <- function(data, effects, model, nbrNodes,
						initgainGlobal, initgainGroupwise,
						priorMu, priorSigma, priorDf, priorKappa,
#						fixedZeros,
						frequentist, incidentalBasicRates, delta,
						nmain, nwarm, lengthPhase1, lengthPhase3,
						prevAns, clusterType=c("PSOCK", "FORK"))
{
	##@precision internal initializeBayes invert z$covtheta
	## avoiding some inversion problems
	## for MoM estimates only
	precision <- function(z){
		require(MASS)
		t(z$dfrac) %*% ginv(z$msfc) %*% z$dfrac
	}

	##@projectEffects internal initializeBayes project model specification
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
	if (!inherits(data,"siena")){stop("The data is not a valid siena object.")}
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
	# Higher than default diagonalize value to increase stability.
	if (ngroups <= 1)
	{
		if (initgainGlobal <= 0)
		{
			startupModel <- sienaAlgorithmCreate(n3=500, nsub=0, cond=FALSE)
		}
		else
		{
			startupModel <- sienaAlgorithmCreate(n3=500, nsub=3, cond=FALSE,
									diagonalize=0.8, firstg=initgainGroupwise)
		}
	}
	else
	{
		if (initgainGlobal <= 0)
		{
			startupModel <- sienaAlgorithmCreate(n3=500, nsub=0, cond=FALSE)
		}
		else
		{
			startupModel <- sienaAlgorithmCreate(n3=500, nsub=3, cond=FALSE,
									diagonalize=0.8, firstg=initgainGlobal)
		}
	}
	cat("Estimate initial global parameters\n")
	startupGlobal <- siena07(startupModel, data=data, effects=effects,
								batch=TRUE, silent=TRUE)
	cat("Initial global estimates\n")
	z$initialResults <- print(startupGlobal)
	cat('\n')
	effects <- updateTheta(effects, startupGlobal)
	z$initialGlobalEstimates <- effects[effects$include,]
	cat("\nmaximum initial global estimate is ",
		max(abs(startupGlobal$theta), na.rm=TRUE), "\n")
	if (max(abs(startupGlobal$theta), na.rm=TRUE) > 40)
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
	# weight for groupwise precision, relative to global precision:
	w <- 0.8
	# number of parameters per group:
	npars <- sum(startupGlobal$effects$group==1)
	if (ngroups >= 2)
	{
		if ((startupGlobal$pp - npars)%%(ngroups-1) != 0)
		{
			stop("Not all groups have the same number of periods.")
		}
		# number of rate parameters per group:
		nrates <- (startupGlobal$pp - npars)%/%(ngroups-1)
	}
	else
	{
		nrates <-  sum(startupGlobal$effects$basicRate)
	}
	rateParameters <- matrix(NA, nrates, ngroups)
	nActors <- rep(0, ngroups)
	# compute covariance matrices for random walk proposals for all groups:
	initialEstimates <- matrix(0, ngroups, startupGlobal$pp)
	factorsBasicRate <- matrix(0, ngroups, startupGlobal$pp)
	z$proposalCov <- list()
	if (initgainGroupwise <= 0)
	{
		startup1Model <- sienaAlgorithmCreate(n3=500, nsub=0, cond=FALSE)
	}
	else
	{
		startup1Model <- sienaAlgorithmCreate(n3=500, nsub=1,
										firstg=initgainGroupwise, cond=FALSE)
	}
	for (i in 1:ngroups)
	{
		nActors[i] <- length(data[[i]]$nodeSets[[1]])
		# coordinates of the global parameters used for this group:
		use <- rep(FALSE, startupGlobal$pp)
		use[!startupGlobal$effects$basicRate] <- TRUE
		rates <- !use
		use[startupGlobal$effects$group==i] <- TRUE
		# project model specification to the single group:
		rates <- rates[use]
		if (ngroups <= 1)
		{
			effects0 <- effects
			startup1 <- startupGlobal
		}
		else
		{
			effects0 <- getEffects(data[[i]])
			effects0 <- projectEffects(effects0,startupGlobal)
			# get the parameter values from the global fit:
			effects0$initialValue[effects0$include] <-
				startupGlobal$theta[use]
			cat("Estimate initial parameters group", i, "\n")
			startup1 <- siena07(startup1Model, data=data[[i]],
					effects=effects0, batch=TRUE, silent=TRUE)
			cat("\nInitial estimate obtained\n", startup1$theta,"\n")
		}
		initialEstimates[i,use] <- startup1$theta
		factorsBasicRate[i,use] <- (1/nActors[i])*startup1$theta
		# non-rate parameters will be dropped from this object, hence the name
		rateParameters[,i] <- trafo(startup1$theta[rates])
		# estimate uncertainties for this group:
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

	if (max(abs(initialEstimates), na.rm=TRUE) > 40)
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

   	z$FRAN <- getFromNamespace(model$FRANname, pkgname)
    z <- initializeFRAN(z, model, data=data, effects=effects,
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
        clusterCall(z$cl, initializeFRAN, z, model,
                    initC = TRUE, profileData=FALSE, returnDeps=FALSE)
		clusterSetRNGStream(z$cl, iseed = as.integer(runif(1,
								max=.Machine$integer.max)))
    }

    z$scaleFactors <- rep(1, z$nGroup)
	z$allScaleFactors <- list()
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
	# Dependent variable:
	z$dependentVariable <- z$effects$name

	if (frequentist)
	{
		if (z$nGroup <= 1)
		{
			cat("\nFrequentist option for sienaBayes is available only for ")
			cat("more than one group.\n")
			frequentist <- FALSE
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
	diag(z$priorSigma)[rates] <- diag(z$priorSigma)[rates] + 0.1
#	if (is.null(fixedZeros))
#	{
#		z$fixedZeros <- matrix(FALSE, z$nGroups, z$TruNumPars)
#		z$someFixedZeros <- FALSE
#	}
#	else
#	{
#		z$fixedZeros <- fixedZeros
#		z$someFixedZeros <- any(fixedZeros)
#	}
	z$frequentist <- frequentist
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
		nwarm <- max(1,
			round(nwarm*nmain/(nwarm + 2*lengthPhase1+lengthPhase3)))
		oldLengthPhase1 <- lengthPhase1
		lengthPhase1 <- max(1,
			round(lengthPhase1*nmain/(nwarm + 2*lengthPhase1+lengthPhase3)))
		lengthPhase3 <- max(1,
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
	z$lengthPhase1 <- lengthPhase1
	z$lengthPhase3 <- lengthPhase3
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
		stop("Function sienaBayes does not allow non-basic rate parameters.")
		}
	# thetaMat must get some NAs so that in sampleparameters it will be clear
	# which values are meaningless and will not be operated upon.
    for (i in 1:z$nGroup)
    {
        use <- rep(FALSE, z$pp)
        use[z$ratePositions[[i]]] <- TRUE
        use[!z$basicRate] <- TRUE
        z$thetaMat[i, !use] <- NA
		factorsBasicRate[i, !use] <- NA
    }
	z$factorsBasicRate <- factorsBasicRate[,z$basicRate]
#cat("factorsBasicRate\n", z$factorsBasicRate,"\n")
#browser()
	meanTruPars <- colMeans(matrix( z$thetaMat[!is.na(z$thetaMat)],
									z$nGroup, z$TruNumPars))
	z$muTemp <- matrix(meanTruPars, z$TruNumPars, 1)
	z$SigmaTemp <- diag(z$TruNumPars)
    z
}


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
    f <- FRANstore()
	theta <- thetaMat[x[1], ]
    .Call("getChainProbabilities", PACKAGE = pkgname, f$pData,
		  f$pModel, as.integer(x[1]), as.integer(x[2]),
		  as.integer(index), f$myeffects, theta, getScores)

}

##@trafo link function rates
# the log link is too strong!
# If the trafo is changed, then in print.summary.sienaBayesFit
# the note about the scale of the basic rate paramaters
# must be changed correspondingly.
#trafo <- function(x){ifelse(x<0.01, 0.1, sqrt(x))}
trafo <- function(x){ifelse(x<0.01, 0.01, x)}
##@antitrafo inverse link function rates
#antitrafo <- function(x){x^2}
antitrafo <- function(x){x}
##@devtrafo derivative link function rates
#devtrafo <- function(x){ifelse(x<0.01, 0.05, 1/(2*sqrt(x)))}
devtrafo <- function(x){1}
#
#createFixedZeroMatrix <- function(x)
#{
#
#}