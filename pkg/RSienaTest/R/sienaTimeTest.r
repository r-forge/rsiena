##*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: http://stat.gamma.rug.nl/siena.html
## *
## * File: sienaTimeTest.r
## *
## * Description: This file contains the sienaTimeTest code for testing the
## * significance of additional time dummies interacted with effects, and
## * sienaTimeFix, which is called to set up time dummy interacted effects.
## ****************************************************************************/
##@sienaTimeTest siena07 Does test for time homogeneity of effects
sienaTimeTest <- function (sienaFit)
{
	observations <- sienaFit$f$observations
	if (observations <=2)
	{
		stop("You must have at least three time periods to test
			 for non-heterogeneity across time.")
	}
## Find out which dummies have already been included in the model:
	indRateEffects <- which(sienaFit$effects$shortName=="Rate")
	indDummiedEffects <- grep("Dummy", sienaFit$effects$effectName)
	indBaseEffects <- setdiff(1:nrow(sienaFit$effects), c(indRateEffects,
														  indDummiedEffects))
	toTest <- array(TRUE, dim=c(length(indBaseEffects),
								sienaFit$f$observations - 2))
	dummyByEffect <- array(0, dim=c(length(indBaseEffects),
									sienaFit$f$observations - 2))
	rownames(toTest) <- sienaFit$effects$effectNumber[indBaseEffects]
	colnames(toTest) <- 2:(sienaFit$f$observations - 1)
	dimnames(dummyByEffect) <- dimnames(toTest)
## Must be able to screen out DummyX ego effects that are fixed:
	dscreen <- which(sienaFit$effects$shortName=='egoX' & sienaFit$effects$fix
					 & length(grep("Dummy", sienaFit$effects$effectName)) > 0)
	if (length(dscreen)==0)
	{
		dscreen <- 99999
	}
	for (i in sienaFit$effects$effectNumber[sienaFit$effects$timeDummy != ',']){
		tmp <- toString(sienaFit$effects$timeDummy[
					   sienaFit$effects$effectNumber == i])
		tmp <- strsplit(tmp, split=",", fixed=TRUE)[[1]]
		if (length(which(!tmp == '')) > 0)
		{
			if (tmp[1]=='isDummy' & !(i %in% sienaFit$effects$
									  effectNumber[dscreen]))
			{
				toTest[rownames(toTest)==as.numeric(tmp[3]),
				colnames(toTest)==as.numeric(tmp[2])] <- FALSE
				dummyByEffect[rownames(toTest)==as.numeric(tmp[3]),
				colnames(toTest)==as.numeric(tmp[2])]  <-
				which(sienaFit$effects$effectNumber[-dscreen]==i)
			}

		}
		else
		{
			next
		}
	}
	nEffects <- length(indBaseEffects) + sum(!toTest)
	nSims <- sienaFit$n3
	nameslist <- list(
					Iteration=paste("it", 1:nSims, sep=""),
					Wave=paste("Wave", 1:(observations - 1), sep=""),
					Effect=sienaFit$effects$effectName
					)
	nDummies <- sum(toTest)
	nTotalEffects <- nDummies + nEffects
	obsStats <- t(sienaFit$targets2[-dscreen, ])
	moment <- sienaFit$sf2[, , -dscreen] - rep(obsStats, each=nSims)
	dummyNames <- rep("", nDummies)
	G <- array(0, dim=c(nSims, observations - 1, nEffects + nDummies))
	G[, , 1:nEffects] <- moment
	inc <- nEffects
	for (i in 1:nrow(toTest))
	{
		for (j in 1:ncol(toTest))
		{
			if (toTest[i, j])
			{
				inc <- inc + 1
				G[, j + 1, inc] <- moment[, j + 1, i]
				dummyNames[inc-nEffects] <- paste("(*)Dummy", j + 1, ":",
												  nameslist$Effect[i], sep="")
			}
		}
	}
	dimnames(G) <- list(nameslist$Iteration, nameslist$Wave,
						c(nameslist$Effect[-dscreen], dummyNames))
	sigma <- cov(apply(G, c(1, 3), sum))
	SF <- array(0, dim=c(nSims, observations - 1, nEffects + nDummies))
	SF[, , 1:nEffects] <- sienaFit$ssc[ , , -dscreen]
	inc <- nEffects
## Add any dummies that have been previously estimated:
	dummyProps <- list()
	for (i in 1:nrow(toTest))
	{
		for (j in 1:ncol(toTest))
		{
			if (toTest[i, j])
			{
				inc <- inc + 1
				SF[, j + 1, inc] <- sienaFit$ssc[, j + 1, i]
				dummyNames[inc-nEffects] <- paste("(*)Dummy", j + 1, ":",
												  nameslist$Effect[i], sep="")
				dummyByEffect[i, j]=inc
				dummyProps$shortName[inc] <- sienaFit$effects$shortName[i]
				dummyProps$type[inc] <- sienaFit$effects$type[i]
				dummyProps$period[inc] <- j + 1
			}
		}
	}
	dimnames(SF) <- list(nameslist$Iteration, nameslist$Wave,
						 c(nameslist$Effect[-dscreen], dummyNames))
	D <- derivativeFromScoresAndDeviations(SF, G)
	fra <- apply(G, 3, sum) / nSims
	doTests <- c(rep(FALSE, nEffects), rep(TRUE, nDummies))
	jointTest <- ScoreTest(nTotalEffects, D, sigma, fra, doTests, maxlike=FALSE)
	jointTestP <- 1 - pchisq(jointTest$testresOverall, nDummies)
	individualTest <- jointTest$testresulto[1:nDummies]
	individualTestP <- 2 * (1-pnorm(abs(individualTest))[1:nDummies])
	rownames(jointTestP) <- c("Joint Significant Test")
	colnames(jointTestP) <- c("p-Val")
	thetaOneStep <- c(sienaFit$theta[-dscreen], rep(0, nDummies)) +
	jointTest$oneStep
	effectTest <- sapply(1:length(indBaseEffects), function (i)
					{
						 doTests <- rep(FALSE, nEffects + nDummies)
						 tmp <- which(dummyProps$shortName ==
									  sienaFit$effects$shortName[i])
						 if (length(tmp) > 0)
						 {
							doTests[tmp] <- TRUE
							test <- ScoreTest(nTotalEffects, D, sigma, fra,
										   doTests, FALSE)
							test$testresOverall
						 }
						 else
						 {
							NA
						 }
					})
	dim(effectTest) <- c(length(indBaseEffects), 1)
	effectTestP <- round(1 - pchisq(effectTest, apply(toTest, 1, sum)), 5)
	rownames(effectTestP) <- nameslist$Effect[indBaseEffects]
	colnames(effectTestP) <- c("p-Val")
	thetaStar <- cbind(c(sienaFit$theta[-dscreen], rep(0, nDummies)),
					   thetaOneStep,
					   round(c(2-2 * pnorm(abs(sienaFit$theta[-dscreen]/
											 sqrt(diag(sienaFit$covtheta)[-dscreen]))),
							   individualTestP), 5))
	colnames(thetaStar) <- c("Initial Est.", "One Step Est.", "p-Value")
	rownames(thetaStar) <- dimnames(SF)[[3]]
	returnObj <- list(
					  JointTest=jointTestP,
					  EffectTest=effectTestP,
					  IndividualTest=thetaStar,
					  JointTestStatistics=jointTest,
					  EffectTestStatistics=effectTest,
					  IndividualTestStatistics=individualTest,
					  CovDummyEst=jointTest$covMatrix,
					  Moments=G,
					  NonRateIndices=indBaseEffects,
					  Waves=dim(G)[2],
					  Sims=dim(G)[1],
					  Effects=dim(G)[3],
					  DummyIndexByEffect=dummyByEffect,
					  DummyStdErr=sqrt(diag(jointTest$covMatrix)),
					  OriginalEffects=nEffects,
					  OriginalThetaStderr=sqrt(diag(sienaFit$covtheta))[-dscreen],
					  SienaFit=sienaFit,
					  DummyProps=dummyProps,
					  ToTest=toTest
					  )
	class(returnObj) <- "sienaTimeTest"
	returnObj
}
summary.sienaTimeTest <- function(object, ...)
{
	if (!inherits(object, "sienaTimeTest"))
	{
		stop("not a legitimate Siena time test object")
	}
	class(object) <- c("summary.sienaTimeTest", class(object))
	object
}
print.summary.sienaTimeTest <- function(x, ...)
{
	if (!inherits(x, "summary.sienaTimeTest"))
	{
		stop("not a legitimate Siena time test summary object")
	}
	print.sienaTimeTest(x)
## Additional output to the print will go in here:
	cat("\nIndividual significance tests and one-step estimators:\n")
	print(x$IndividualTest)
	cat("\nParameter-wise joint significance tests (i.e. each
		parameter across all dummies):\n")
	print(x$EffectTest)
	if (x$Waves <=2)
	{
		cat("\n\nNote that these parameter-wise tests have a different
			form than the individual tests, thus testing with 3 observations
			may yield different individual and parameter-wise values.\n\n")
	}
	tmp <- paste(" (", 1:length(rownames(x$IndividualTest)), ") ",
				 rownames(x$IndividualTest), "\n", sep="")
	cat("\nUse the following indices for plotting the pairwise moment
		correlations:\n", tmp)
	tmp <- paste(" (", 1:length(x$NonRateIndices), ") ",
				 rownames(x$IndividualTest)[x$NonRateIndices], "\n", sep="")
	cat("\nUse the following indices for plotting the effect-wise
		fitted parameters:\n", tmp)
	effectNames <- rownames(x$IndividualTest)
	dummies <- grepl("Dummy", effectNames)
	dummyIndex <- paste(" (", 1:sum(dummies), ") ", effectNames[dummies],
						"\n", sep="")
	cat("\nIf you would like to fit time dummies to your model, use the
		following indices:\n", dummyIndex)
	cat("\nType \"?sienaTimeTest\" for more information on this output.\n")
	invisible(x)
}
print.sienaTimeTest <- function(x, ...)
{
	if (!inherits(x, "sienaTimeTest"))
	{
		stop("not a legitimate Siena time test object")
	}
	effectNames <- rownames(x$IndividualTest)
	dummies <- grepl("Dummy", effectNames)
	dummyIndex <- paste(" (", 1:sum(dummies), ") ", effectNames[dummies],
						"\n", sep="")
	cat("Joint significance test of the dummy parameters:\np-Val = ",
		x$JointTest,
		", \nWhere H0: The following parameters are zero:\n",
		dummyIndex
		)
	invisible(x)
}
plot.sienaTimeTest <- function(x, pairwise=FALSE, effects=1:2,
	dims=c(2, 1), scale=.2, plevels=c(.1, .05, .025),
	multiplot=FALSE, ...)
{
	require(lattice)
	timetest <- x
	if (pairwise)
	{
## On a "pairwise" call, print a pairwise plot of moments
		if (length(intersect(effects, 1:timetest$OriginalEffects))!=
			length(effects))
		{
			cat("Detected an error with the effects included. For a
				parameter-plot, use the following indices:")
			tmp <- paste(" (", 1:length(rownames(timetest$IndividualTest)),
						 ") ", rownames(timetest$IndividualTest), "\n", sep="")
			cat("\nUse the following indices for plotting the pairwise
				moment correlations:\n", tmp)
			stop(" ")
		}
		if (length(effects)==0)
		{
			x <- timetest$Moments

		}
		else
		{
			if (class(effects)!="integer")
			{
				stop("Effects is not a vector of integers.")
			}
			x <- timetest$Moments[, , effects]
		}
		panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
		{
			usr <- par("usr"); on.exit(par(usr))
			par(usr = c(0, 1, 0, 1))
			r <- abs(cor(x, y))
			txt <- format(c(r, 0.123456789), digits=digits)[1]
			txt <- paste(prefix, txt, sep="")
			if(missing(cex.cor)) cex.cor <- 0.8 / strwidth(txt)
			text(0.5, 0.5, txt, cex = cex.cor * r)
		}
		pairs(apply(x, c(1, 3), sum),
			  lower.panel=panel.smooth,
			  upper.panel=panel.cor,
			  pch=5, ...)

	}
	else
	{
## Otherwise, make the parameter plots:
		if (length(intersect(effects, 1:1:nrow(x$ToTest)))!=length(effects))
		{
			cat("Detected an error with the effects included. For a parameter
				-plot, use the following indices:")
			tmp <- paste(" (", 1:length(timetest$NonRateIndices), ") ",
						 rownames(timetest$IndividualTest)
						 [timetest$NonRateIndices], "\n", sep="")
			cat("\nUse the following indices for plotting the effect-wise
				fitted parameters:\n", tmp)
			stop(" ")
		}
		if (multiplot)
		{
			dims=c(1,1)
			nplots=length(effects)
		}
		else
		{
		nplots=(dims[1] * dims[2])
		}
		if (length(effects) > nplots)
		{
			stop("You have included space for ", nplots, " plots, but have
				 requested ",
				 length(effects), " effects to be plotted. Use dims=c(x, y).")
		}
		xaxis <- 1:timetest$Waves
		if (length(effects)==1)
		{
			yaxis <- timetest$IndividualTest[as.vector(c(effects,
														 timetest$DummyIndexByEffect[effects, ])), 2]
			dim(yaxis) <- c(1, timetest$Waves)

		}
		else
		{
			yaxis <- timetest$IndividualTest[as.vector(t(cbind(effects,
															   timetest$DummyIndexByEffect[effects, ]))), 2]
			yaxis <- matrix(yaxis, nrow=length(effects), ncol=timetest$Waves,
							byrow=TRUE)
		}
		rownames(yaxis) <- rownames(timetest$IndividualTest)[effects]
		colnames(yaxis) <- 1:timetest$Waves
		vals <- yaxis
		basevals <- array(yaxis[, 1], dim=dim(yaxis))
		basevals[, 1] <- 0
		yaxis <- yaxis + basevals
		pvals <- timetest$IndividualTest[c(effects, as.vector(
															  t(timetest$DummyIndexByEffect[effects, ]))), 3]
		dummysd <- abs(c(vals / qnorm(1 - pvals / 2)))
		dummysd[effects] <- timetest$OriginalThetaStderr[effects]
		dim(dummysd) <- c(length(effects), timetest$Waves)
		dim(pvals) <- dim(dummysd)
		dim(vals) <- dim(dummysd)
		rownames(dummysd) <- rownames(timetest$IndividualTest)[effects]
		colnames(dummysd) <- 1:timetest$Waves
##Function to print a panel:
		makeplot <- function (i)
		{
			ymin=min(yaxis[i, ] - scale * abs(yaxis[i, ]))
			ymax=max(yaxis[i, ] + scale * abs(yaxis[i, ]))
			xyplot(yaxis[i, ] ~ xaxis,
				   type = "p", main = rownames(timetest$IndividualTest)
				   [timetest$NonRateIndices[effects[i]]],
				   sub=paste("p=", timetest$EffectTest[effects[i]]), bty="n",
				   xlab="Wave", ylab="Parameter Value", auto.key=TRUE,
				   ylim=c(ymin, ymax), xlim=c(0, length(xaxis) + 1),
				   panel=function(x, y){
                       for (j in 1:length(x))
                       {
                           if ( c(FALSE, timetest$ToTest[effects[i], ])[j] )
                           {
                               tmp="red"

                           }
                           else
                           {
                               tmp="gray"
                           }
                           l <- yaxis[i, j] - abs(qnorm(plevels[1]  /  2,
                                                        sd=dummysd[i, j]))
                           u <- yaxis[i, j] + abs(qnorm(plevels[1] / 2,
                                                        sd=dummysd[i, j]))
                           panel.xyplot(c(x[j], x[j]), c(l, u), reference=TRUE,
                                        col=tmp, alpha=.50,
                                        type="l", lend=1, lwd=10)
                           panel.xyplot(c(x[j], x[j]), c(l, u), reference=TRUE,
                                        col=tmp, alpha=.75,
                                        type="p", pch=45, cex=3)
                           l <- yaxis[i, j] - abs(qnorm(plevels[2] / 2,
                                                        sd=dummysd[i, j]))
                           u <- yaxis[i, j] + abs(qnorm(plevels[2] / 2,
                                                        sd=dummysd[i, j]))
                           panel.xyplot(c(x[j], x[j]), c(l, u), reference=TRUE,
                                        col=tmp, alpha=.50,
                                        type="l", lend=1, lwd=10)
                           panel.xyplot(c(x[j], x[j]), c(l, u), reference=TRUE,
                                        col=tmp, alpha=.75,
                                        type="p", pch=45, cex=3)
                           l <- yaxis[i, j] - abs(qnorm(plevels[3] / 2,
                                                        sd=dummysd[i, j]))
                           u <- yaxis[i, j] + abs(qnorm(plevels[3] / 2,
                                                        sd=dummysd[i, j]))
                           panel.xyplot(c(x[j], x[j]), c(l, u), reference=TRUE,
                                        col=tmp, alpha=.25,
                                        type="l", lend=1, lwd=10)
                       }
                       panel.xyplot(x, y, type="s", reference=TRUE,
                                    col="black", alpha=.75, pch=2)
                       panel.xyplot(x, y, type="p", reference=TRUE, pch=20,
                                    col=1)
                       panel.abline(a=timetest$IndividualTest[effects[i], 1],
                                    reference=TRUE,
                                    col="black", lwd=2, alpha=.75)

				   }, ...)
		}

		if (length(effects) > 1 & !multiplot)
		{
            print(makeplot(1), newpage=TRUE, more=TRUE, split=c(1, 1, dims[1],
                                                        dims[2]))
			if(dims[1] > dims[2])
			{
				col=1
				row=2

			}
			else
			{
				row=1
				col=2
			}
			if (length(effects) > 2)
			{
				for (i in 2:(length(effects)-1))
				{
					print(makeplot(i), more=TRUE, split=c(row, col, dims[1],
														  dims[2]))
					col <- col + 1
					if(col > dims[2])
					{
						col=1
						row <- row + 1
					}
				}
			}
			print(makeplot(length(effects)), split=c(row, col, dims[1],
                                             dims[2]))
		}

		else if (length(effects) > 1 & multiplot)
		{
			for (i in 1:(length(effects)))
			{
				dev.new()
				print(makeplot(i), newpage=TRUE, more=FALSE,
                      split=c(1, 1, 1, 1))
			}
		}
		else
		{
			print(makeplot(1), newpage=TRUE, more=FALSE, split=c(1, 1, 1, 1))
		}
	}
}
##@sienaTimeFix siena07 Adds time dummy terms to the effects object
sienaTimeFix <- function(effects, data)
{
    if (inherits(data, "sienaGroup"))
    {
        warning("Time dummy not implemented for multi-group projects")
        effects$timeDummy <- ","
    }
    else
    {
        observations <- data$observations - 1
        if (observations < 2 && any(effects$timeDummy != ","))
        {
            warning("Time dummies not relevant with only 2 periods")
            effects$timeDummy <- ","
        }
    }
    use <- effects$name == effects$name[1]
    if (any(effects$timeDummy[!use] != ","))
    {
        warning("Time dummy only implemented for first dependent variable")
        effects$timeDummy[!use] <- ","
    }
    if (length(unique(effects$groupName)) > 1)
    {
        warning("Time dummy not implemented for multi-group projects")
        effects$timeDummy <- ","
    }
    covar <- effects$interaction1 != ""
    if (any(effects$timeDummy[covar] != ","))
    {
        warning("Time dummy not implemented for covariate effects")
        effects$timeDummy[covar] <- ","
    }
    eval <- effects$type =="eval"
	if (any(effects$timeDummy[!eval] !=','))
	{
		warning("Time dummy effects are only implemented",
                " for one mode network effects of type eval.")
        effects$timeDummy[!eval] <- ","
	}
	if (all(effects$timeDummy==',') )
	{
##		No time dummy interactions to add, so kick the inputs back.
		return(list (effects=effects, data=data))
	}
	else
	{
		alreadyDummied <- grep("isDummy", effects$timeDummy)
		effects$timeDummy[effects$timeDummy=="all"]  <-
            paste(2:(data$observations-1), collapse = ",")
		if (length(alreadyDummied)  >  0)
		{
## Just remove those effects that have already been dummied so as to
## not messy things up. The assumption is that the user will retain
## all of the previous dummied effects within the column.
			effects <- effects[-alreadyDummied, ]
		}
		dummiedEffects <- effects$effectNumber[effects$timeDummy != ',']
		covToAdd <- NULL
		dummyCombos <- list()
		ctr=1
## This might need to be changed for sienaGroup:
		nact=dim(data$depvars[[1]])[1]
		nper=dim(data$depvars[[1]])[3]
		for (i in dummiedEffects)
		{
## Get the time periods that we want dummied for effect i:
			tmp <- toString(effects$timeDummy[effects$effectNumber == i])
			tmp <- strsplit(tmp, split=",", fixed=TRUE)[[1]]
			if (length(which(!tmp == '')) > 0)
			{
				tmp=as.numeric(tmp)
				tmp=tmp[tmp > 1 & tmp  <  nper]

			}
			else
			{
				next
			}
			if (length(which(!is.numeric(tmp))) > 0)
			{
				stop("Invalid input for time dummy column of effects object:", tmp)
			}
			if (length(tmp) > 0)
			{
				dummyCombos[[ctr]]=list(effectNumber=i, periods=tmp)
				ctr=ctr + 1
				covToAdd <- unique(c(covToAdd, tmp))
			}
		}
## Add the required covariate effects to the effect objects
		ctr <- length(data$vCovars) + 1
		for (i in covToAdd)
		{
			dname <- paste("Dummy", i, sep='')
			tmp <- array(0, c(nact, nper-1))
			tmp[, i]=1
			tmp <- varCovar(tmp)
			tmp <- addAttributes.varCovar(tmp, name=dname)
			data$vCovars[[ctr]] <- tmp
			names(data$vCovars)[ctr] <- dname
			ctr <- ctr + 1
			tmprow <- allEffects[allEffects$functionName==
			'Sum of outdegrees x xxxxxx' & allEffects$type=='eval'
			& allEffects$effectGroup=='covarNonSymmetricObjective', ]
			tmprow$name <- effects$name[effects$shortName=='density' &
			effects$type=='eval'][1]
			tmprow$effectFn <- 'NULL'
			tmprow$statisticFn <- 'NULL'
			tmprow$netType <- 'oneMode'
			tmprow$groupName <- 'Group1'
			tmprow$group <- 1
			tmprow$fix <- TRUE
			tmprow$include <- TRUE
			tmprow$effectNumber <- max(effects$effectNumber) + 1
			tmprow <- tmprow[, colnames(effects)]
			tmprow$effectName <- gsub('xxxxxx', dname, tmprow$effectName)
			tmprow$functionName <- gsub('xxxxxx', dname, tmprow$functionName)
			tmprow$interaction1 <- dname
			tmprow$timeDummy <- paste('isDummy', i,
									  effects$effectNumber[effects$shortName=='density' &
									  effects$type=='eval'], sep=',')
			rownames(tmprow) <- dname
			effects <- rbind(effects, tmprow)
		}
		for (i in 1:length(dummyCombos))
		{
			baseNum=dummyCombos[[i]]$effectNumber
			for (j in 1:length(dummyCombos[[i]]$periods))
			{
				dname <- paste("Dummy", dummyCombos[[i]]$periods[j], sep="")
				dummyNum <- effects$effectNumber[rownames(effects)==dname]
				if (effects$shortName[baseNum] != 'density')
				{
## Make a user specified interaction
## for the time dummy interacted effect
					tmprow <- allEffects[allEffects$shortName=='unspInt'
					& allEffects$type=='eval'
					& allEffects$effectGroup=='unspecifiedNetInteraction', ]
					tmprow$name <- effects$name[effects$effectNumber==baseNum]
					tmprow$effectFn <- 'NULL'
					tmprow$statisticFn <- 'NULL'
					tmprow$netType <- 'oneMode'
					tmprow$groupName <- 'Group1'
					tmprow$group <- 1
					tmprow$fix <- FALSE
					tmprow$include <- TRUE
					tmprow$effectNumber <- max(effects$effectNumber) + 1
					tmprow <- tmprow[, colnames(effects)]
					tmprow$effectName <- 'unspecified interaction effect'
					tmprow$functionName <- 'unspecified interaction statistic'
					rownames(tmprow) <- paste(dname, baseNum, sep='.')
					tmprow$effect1 <- baseNum
					tmprow$effect2 <- dummyNum
					tmprow$timeDummy <- paste('isDummy',
											  dummyCombos[[i]]$periods[j], baseNum, sep=',')
					effects <- rbind(effects, tmprow)

				}
				else
				{
					effects$fix[effects$effectNumber==dummyNum] <- FALSE
				}
			}
		}
		list(effects=effects, data=data)
	}
}

