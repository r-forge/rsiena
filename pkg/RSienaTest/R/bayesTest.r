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


##@multipleBayesTest Tests parameters of sienaBayesFit objects
multipleBayesTest <- function(z, testedPar, nfirst=z$nwarm, tested0=0, ndigits=4){
    theEffects <- z$effects
#	efNames <- format(paste(ifelse(theEffects$type == "creation",
#									"creat", theEffects$type),
#							theEffects$effectName)[!z$basicRate])
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
	# Put mu and eta parameters in their proper places:
	z$ThinObjective <- matrix(NA, nmax,
		dim(z$ThinPosteriorEta)[2] +
			dim(z$ThinPosteriorMu[,z$objectiveInVarying, drop=FALSE])[2])
	vio <- vec1[vec1 %in% c(3,4)] == 3 # could be called z$varyingInObjective
	z$ThinObjective[1:nmax, vio] <-
			z$ThinPosteriorMu[nfirst:ntot, z$objectiveInVarying, drop=FALSE]
	z$ThinObjective[1:nmax, !vio] <- z$ThinPosteriorEta[nfirst:ntot,]
	p <- dim(z$ThinObjective)[2]
	if (inherits(testedPar, "matrix"))
	{
		A <- testedPar
		if (dim(A)[2] != p	)
		{
			stop(paste("If testedPar is a matrix, it should have",p,"columns."))
		}
		effeNames <- efNames
		testedNumber <- dim(A)[1]
		efNames <- rep(NA,testedNumber)
		for (i in (1:testedNumber))
		{
			efNames[i] <-
				paste(paste(A[i,],'*','parameter',1:p)[A[i,]!=0], collapse=" + ")
		}
		posteriorSample <- z$ThinObjective[,, drop=FALSE] %*% t(A)
		usedEffects <- apply(A,2,function(x){any(x!=0)})
		theUsedEffects <- paste(which(usedEffects),
							'. ',effeNames[usedEffects],sep='')
	}
	else
	{
		testedSet <- testedPar
		if ((min(testedSet) < 1) | (max(testedSet) > p))
		{
			stop(paste('testedPar should be a matrix or a set of numbers between 1 and', p))
		}
		testedNumber <- length(testedSet)
		efNames <- efNames[testedSet]
		posteriorSample <- z$ThinObjective[,testedSet, drop=FALSE]
		usedEffects <- (1:p) %in% testedPar
		theUsedEffects <- paste(which(usedEffects),'. ',efNames,sep='')
	}
	if (!(length(tested0) %in% c(1,testedNumber)))
	{
		stop(paste('tested0 must be one number or a vector of',
				testedNumber, 'numbers'))
	}
	if (length(tested0) == 1)
	{
		testedvec <- rep(tested0, testedNumber)
	}
	else
	{
		testedvec <- tested0
	}
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
				nullValue=testedvec, effectNames=efNames, theUsedEffects=theUsedEffects,
				posteriorSample=posteriorSample)
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
	if (p >= 2)
	{
		mydf[,1] <- paste('. ',x$effectNames,"= ")
	}
	else
	{
		mydf[,1] <- paste(x$effectNames,"= ")
	}
	mydf[,2] <- x$nullValue
	names(mydf) <- c(' ', ' ')
	if (!(is.null(x$theUsedEffects)))
	{
		cat('Used effects:\n')
		for (i in 1:length(x$theUsedEffects)){cat(x$theUsedEffects[i],'\n')}
	}
	cat("Tested hypothesis:\n")
	if (p >= 2) {print(mydf)} else {cat(as.matrix(mydf),'\n')}
	cat("\nposterior p-value:", round(x$prob, 4), "\n")
	cat("\ntest statistic:", round(x$chisquared, 2), "d.f. =",p,".\n\n")
	# Construct pattern of all sign combinations of p numbers
	s <- matrix(0,1,0)
	for (i in 1:p){
		pp <- 2^(i-1)
		s <- cbind(rbind(s,s), c(rep(1,pp),rep(-1,pp)))
		}
	propSignPattern <- function(si){
		mean(apply(x$posteriorSample, 1, function(x){identical(sign(x), si)}))
	}
	propSigns <- apply(s,1,propSignPattern)
	cat("Posterior proportion of sign patterns:\n")
	cat("", 1:p, "\n", sep="   ")
	pp <- dim(s)[1]
	# cbind transforms T to 1 and F to 0
	apply(cbind(s, propSigns), 1,
		function(y){cat(' ', ifelse((y[1:p] > 0.5), '>0 ', '<0 '),
				formatC(y[p+1], digits=3),'\n')})
	cat("\n")
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
	xlim[1] <- min(0, xlim[1])
# plot density so that truncation at 0 is clear;
# since the distances are those from the posterior mean,
# it is very unlikely that 0 is not in the support of the distribution.
	d1 <- density(c(post,-post))
	d1$x <- abs(d1$x)
	order.x <- order(d1$x)
	d1$x <- c(0,0,d1$x[order.x])
	d1$y <- c(0,d1$y[1],d1$y[order.x])
	if (is.null(ylim)){ylim <- stretch2(range(d1$y))}
	if (is.null(main)){main <- "posterior distances"}
	if (x$chisquared > xlim[2]){
		cat('note: observed chi-squared =',round(x$chisquared,1),
				'outside plot window.\n')}
	plot(d1, xlim=xlim, ylim=ylim, xlab='distance', ylab="density", col=4,
			cex=4, cex.lab=2, cex.main=2, main=main, lwd=2, ...)
	lines(c(x$chisquared, x$chisquared), c(0, max(d1$y)),lwd=2)
}

##@extract.sienaBayes extracts samples from sienaBayesFit objects
extract.sienaBayes <- function(zlist, nfirst=zlist[[1]]$nwarm+1, extracted,
                   sdLog=TRUE){
	getNames <- function(x){
# effect names without duplicated rate parameters
# and with "unspecified interaction" replaced by
# information about the effects in question
		b <- x$basicRate
		tpar<- rep(NA,length(b))
# True Parameters, i.e., all except rate parameters for groups 2 and up.
		for (i in (2:length(b))){tpar[i] <- !b[i]|(b[i]&!b[i-1])}
		tpar[1] <- TRUE
# Take away the ' (period 1)' in the first rate parameter
		sub(' (period 1)','', x$requestedEffects$effectName[tpar], fixed=TRUE)
	}
	if (!(is.list(zlist)))
	{
		stop('zlist must be a list of sienaBayesFit objects')
	}
	if (!all(sapply(zlist, function(z){inherits(z, "sienaBayesFit")})))
	{
		stop('all elements of zlist must be sienaBayesFit objects')
	}

#browser()
	niter <- sapply(zlist,function(z){dim(z$ThinPosteriorMu)[1]})
	nit <- niter[1] - nfirst + 1

	if (extracted %in% c("all", "rates")){
# rate parameters
		EffName <- zlist[[1]]$effectName
		indices <- sort(which(!zlist[[1]]$generalParametersInGroup))
		nGroups <- zlist[[1]]$nGroup
		nind <- length(indices)
		res1 <- array(dim = c(nit, length(zlist), nGroups*length(indices)))
		fName1 <- rep('',dim(res1)[3])
		for (h in 1:length(zlist)){for (i in 1:nGroups){
			mini <- (i-1)*nind + 1
			maxi <- i*nind
			res1[,h,mini:maxi] <-
				zlist[[h]]$ThinParameters[(niter[h]-nit+1):niter[h],
												i,indices,drop=FALSE]
			fName1[mini:maxi] <- paste(EffName[indices],i)
		}}
	}

	if (extracted %in% c("all", "objective", "varying")){
# mu
		nind <- sum(zlist[[1]]$varyingParametersInGroup)
		if (nind <= 0)
		{
			cat("Note: no varying parameters.\n")
		}
		else
		{
			EffName <- getNames(zlist[[1]])[zlist[[1]]$varyingParametersInGroup]
			res2 <- array(dim = c(nit, length(zlist), (2*nind)))
			for (h in 1:length(zlist)){
			 res2[,h,1:nind] <-
				zlist[[h]]$ThinPosteriorMu[(niter[h]-nit+1):niter[h], ,drop=FALSE]
			 postsig <- zlist[[h]]$ThinPosteriorSigma
			 postsigg <- matrix(NA,dim(postsig)[1],nind)
			 if (sdLog){
				for (k in 1:nind){postsigg[,k] <- 0.5*log(postsig[,k,k])}
			 }
			 else {
				for (k in 1:nind){postsigg[,k] <- sqrt(postsig[,k,k])}
			 }
			 res2[,h,(nind+1):(2*nind)] <-
				postsigg[(niter[h]-nit+1):niter[h],,drop=FALSE]
			}
			fName2 <- rep('',2*nind)
			fName2[1:nind] <- paste(EffName,'mu')
			if (sdLog){
				fName2[(nind+1):(2*nind)] <- paste(EffName,'log(s.d.)')
			}
			else {
				fName2[(nind+1):(2*nind)] <- paste(EffName,'s.d.')
			}
		}
	}

	if (extracted %in% c("all", "objective", "non-varying")){
# eta
		nind <- sum(!zlist[[1]]$varyingParametersInGroup)
		if (nind <= 0)
		{
			cat("Note: no non-varying parameters.\n")
		}
		else
		{
			EffName <- getNames(zlist[[1]])[!zlist[[1]]$varyingParametersInGroup]
			res3 <- array(dim = c(nit, length(zlist), nind))
			for (h in 1:length(zlist)){
			res3[,h,1:nind] <-
				zlist[[h]]$ThinPosteriorEta[(niter[h]-nit+1):niter[h], ,drop=FALSE]
			}
			fName3 <- EffName
		}
	}

	res <- NULL
	switch(extracted,
		'all' = {res <- array(c(res1,res2,res3),
							dim=c(dim(res1)[1], dim(res1)[2],
								dim(res1)[3]+dim(res2)[3]+dim(res3)[3]))
				fName <- c(fName1, fName2, fName3)},
		'rates' = {res <- res1
					fName <- fName1},
		'varying' = {res <- res2
					fName <- fName2},
		'non-varying' = {res <- res3
					fName <- fName3},
		'objective' = {res <- array(c(res2,res3),
							dim=c(dim(res2)[1], dim(res2)[2],
									dim(res2)[3]+dim(res3)[3]))
						fName <- c(fName2, fName3)},
		cat('wrong parameter given: extracted =',extracted,'\n'))
	if (!is.null(res)){
		dimnames(res) <- list(1:dim(res)[1], paste('chain', 1:length(zlist)),
								fName)
		if (any(is.na(res))){
			cat('warning: the extracted array contains NA elements\n')
		}
	}
	res
}