#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snijders/siena/
# *
# * File: sienaRI.r
# *
# * Description: Used to determine, print, and plots relative importances of effects
# * for potential decisions of actors at observation moments.
# *****************************************************************************/

##@sienaRI
sienaRI <- function(data, ans=NULL, theta=NULL, algorithm=NULL, effects=NULL)
{
	if (!inherits(data, "siena"))
	{
		stop("no a legitimate Siena data specification")
	}
	if(!is.null(ans))
	{
		if (!inherits(ans, "sienaFit"))
		{
			stop(paste("ans is not a legitimate Siena fit object", sep=""))
		}
		if(!is.null(algorithm)||!is.null(theta)||!is.null(effects))
		{
			warning(paste("some information are multiply defined \n",
		"results will be based on 'theta', 'algorithm', and 'effects'\n",
		"stored in 'ans' (as 'ans$theta', 'ans$x', 'ans$effects')\n", sep=""))
		}
		if (sum(ans$effects$include==TRUE &
			(ans$effects$type =="endow"|ans$effects$type =="creation")) > 0)
			{
stop("sienaRI does not yet work for models that contain endowment or creation effects")
			}
		contributions <- getChangeContributions(algorithm = ans$x, data = data,
										effects = ans$effects)
		RI <- expectedRelativeImportance(conts = contributions,
								effects = ans$effects, theta =ans$theta)
	}else{
		if (!inherits(algorithm, "sienaAlgorithm"))
		{
	stop(paste("algorithm is not a legitimate Siena algorithm specification", sep=""))
		}
		algo <- algorithm
		if (!inherits(effects, "sienaEffects"))
		{
			stop(paste("effects is not a legitimate Siena effects object", sep=""))
		}
		if(sum(effects$include==TRUE &
					(effects$type =="endow"|effects$type =="creation")) > 0)
		{
	stop("sienaRI does not yet work for models containinf endowment or creation effects")
		}
		effs <- effects
		if (!is.numeric(theta))
		{
			stop("theta is not a legitimate parameter vector")
		}
		if(length(theta) != sum(effs$include==TRUE & effs$type!="rate"))
		{
			if(length(theta) != sum(effs$include==TRUE))
			{
	stop("theta is not a legitimate parameter vector \n number of parameters has to match number of effects")
			}
	warning(paste("length of theta does not match the number of objective function effects\n",
				"theta is treated as if containing rate parameters"))
			paras <- theta
			## all necessary information available
			## call getChangeContributions
			contributions <- getChangeContributions(algorithm = algo,
									data = data, effects = effs)
			RI <- expectedRelativeImportance(conts = contributions,
									effects = effs, theta = paras)
		}else{
			paras <- theta
			## all necessary information available
			## call getChangeContributions
			contributions <- getChangeContributions(algorithm = algo,
									data = data, effects = effs)
			RI <- expectedRelativeImportance(conts = contributions,
									effects = effs, theta = paras)
		}
	}
	RI
}

##@getChangeContributions. Use as RSiena:::getChangeContributions
getChangeContributions <- function(algorithm, data, effects)
{
	## The following initializations data, effects, and model
	## for calling "getTargets" in "siena07.setup.h"
	## is more or less copied from "getTargets" in "getTargets.r".
	## However, some modifications have been necessary to get it to work.
	f <- unpackData(data,algorithm)

	effects <- effects[effects$include,]
	if (!is.null(algorithm$settings))
	{
		stop('not implemented: RI together with settings')
		# effects <- addSettingsEffects(effects, algorithm)
	}
	else
	{
		effects$setting <- rep("", nrow(effects))
	}
	pData <- .Call('setupData', PACKAGE=pkgname,
			list(as.integer(f$observations)),
			list(f$nodeSets))
	## register a finalizer
	ans <- reg.finalizer(pData, clearData, onexit = FALSE)
	ans<- .Call('OneMode', PACKAGE=pkgname,
			pData, list(f$nets))
	ans<- .Call('Behavior', PACKAGE=pkgname, pData,
			list(f$behavs))
	ans<-.Call('ConstantCovariates', PACKAGE=pkgname,
			pData, list(f$cCovars))
	ans<-.Call('ChangingCovariates',PACKAGE=pkgname,
			pData,list(f$vCovars))
	ans<-.Call('DyadicCovariates',PACKAGE=pkgname,
			pData,list(f$dycCovars))
	ans<-.Call('ChangingDyadicCovariates',PACKAGE=pkgname,
			pData, list(f$dyvCovars))

	storage.mode(effects$parm) <- 'integer'
	storage.mode(effects$group) <- 'integer'
	storage.mode(effects$period) <- 'integer'

	effects$effectPtr <- rep(NA, nrow(effects))
	depvarnames <- names(data$depvars)
	tmpeffects <- split(effects, effects$name)
	myeffectsOrder <- match(depvarnames, names(tmpeffects))
	ans <- .Call("effects", PACKAGE=pkgname, pData, tmpeffects)
	pModel <- ans[[1]][[1]]
	for (i in 1:length(ans[[2]]))
	{
		effectPtr <- ans[[2]][[i]]
		tmpeffects[[i]]$effectPtr <- effectPtr
	}
	myeffects <- tmpeffects
	for(i in 1:length(myeffectsOrder)){
		myeffects[[i]]<-tmpeffects[[myeffectsOrder[i]]]
	}
	ans <- .Call("getTargets", PACKAGE=pkgname, pData, pModel, myeffects,
			parallelrun=TRUE, returnActorStatistics=FALSE,
			returnStaticChangeContributions=TRUE)
	ans
}

expectedRelativeImportance <- function(conts, effects, theta,
												effectNames = NULL)
{
	waves <- length(conts[[1]])
	effects <- effects[effects$include == TRUE,]
	noRate <- effects$type != "rate"
	effects <- effects[noRate,]
	if(sum(noRate)!=length(theta))
	{
		theta <- theta[noRate]
	}
	effectNa <- attr(conts,"effectNames")
	effectTypes <- attr(conts,"effectTypes")
	networkNames <- attr(conts,"networkNames")
	networkTypes <- attr(conts,"networkTypes")
	networkInteraction <- effects$interaction1
	effectIds <- paste(effectNa,effectTypes,networkInteraction, sep = ".")

	currentDepName <- ""
	depNumber <- 0
	for(eff in 1:length(effectIds))
	{
		if(networkNames[eff] != currentDepName)
		{
			currentDepName <- networkNames[eff]
			actors <- length(conts[[1]][[1]][[1]])
			if(networkTypes[eff] == "oneMode")
			{
				choices <- actors
			}else if(networkTypes[eff] == "behavior"){
				choices <- 3
			}else{
	stop("so far, sienaRI works only for dependent variables of type 'oneMode' or 'behavior'")
			}
			depNumber <- depNumber + 1
			currentDepEffs <- effects$name == currentDepName
			effNumber <- sum(currentDepEffs)

#			RIs <- data.frame(row.names = effectIds[currentDepEffs])
#			RIs <- cbind(RIs, matrix(0, nrow=effNumber, ncol = actors))
			entropies <- vector(mode="numeric", length = actors)

#			currentDepObjEffsNames <- paste(effects$shortName[currentDepEffs],
#				effects$type[currentDepEffs],effects$interaction1[currentDepEffs],sep=".")
#			otherObjEffsNames <- paste(effects$shortName[!currentDepEffs],
#				effects$type[!currentDepEffs],effects$interaction1[!currentDepEffs],sep=".")

			expectedRI <- list()
			expectedI <- list()
			RIActors <- list()
			IActors <- list()
			absoluteSumActors <- list()
			entropyActors <-list()
			for(w in 1:waves)
			{
				currentDepEffectContributions <- conts[[1]][[w]][currentDepEffs]
				currentDepEffectContributions <-
					sapply(lapply(currentDepEffectContributions, unlist),
						matrix, nrow=actors, ncol=choices, byrow=TRUE,
														simplify="array")

				distributions <-
					apply(apply(currentDepEffectContributions, c(2,1), as.matrix),
						3, calculateDistributions, theta[which(currentDepEffs)])
				distributions <-
					lapply(apply(distributions, 2, list),
						function(x){matrix(x[[1]], nrow=effNumber+1,
											ncol=choices, byrow=F)})

				entropy_vector <- unlist(lapply(distributions,
										function(x){entropy(x[1,])}))
				## If one wishes another measure than the
				## L^1-difference between distributions, here is
				## the right place to call some new function instead of "L1D".
				RIs_list <- lapply(distributions,function(x){L1D(x[1,],
														x[2:dim(x)[1],])})
				RIs_matrix <-(matrix(unlist(RIs_list),nrow=effNumber,
														ncol=actors, byrow=F))

#				RIs <- RIs_matrix
				entropies <- entropy_vector
				# divide by column sums:
				RIActors[[w]] <- apply(RIs_matrix, 2, function(x){x/sum(x)})
				absoluteSumActors[[w]] <- colSums(RIs_matrix)
				entropyActors[[w]] <- entropies
				expectedRI[[w]] <- rowSums(RIActors[[w]] )/dim(RIActors[[w]])[2]
				IActors[[w]] <- RIs_matrix
				expectedI[[w]] <- rowMeans(RIs_matrix)
			}
			RItmp <- NULL
			RItmp$dependentVariable <- currentDepName
			RItmp$expectedRI <- expectedRI
			RItmp$RIActors <- RIActors
			RItmp$expectedI <- expectedI
			RItmp$IActors <- IActors
			RItmp$absoluteSumActors <- absoluteSumActors
			RItmp$entropyActors <- entropyActors
			if(!is.null(effectNames))
			{
				RItmp$effectNames <- effectNames[currentDepEffs]
			}else{
				RItmp$effectNames <-
					paste(effectTypes[currentDepEffs], " ",
						effects$effectName[currentDepEffs], sep="")
			}
			class(RItmp) <- "sienaRI"
			if(depNumber == 1){
				RI <- RItmp
			}else if(depNumber == 2){
				RItmp1 <- RI
				RI <- list()
				RI[[1]]<-RItmp1
				RI[[2]]<-RItmp
			}else{
				RI[[depNumber]]<-RItmp
			}
		}
	}
	if(depNumber>1)
	{
		warning(paste("more than one dependent variable\n",
				"return value is therefore not of class 'sienaRI'\n",
				"but a list of objects of class 'sienaRI'.\n"))
	}
	RI
}

calculateDistributions <- function(effectContributions = NULL, theta = NULL)
{
	effects <- dim(effectContributions)[1]
	choices <- dim(effectContributions)[2]
	effectContributions[effectContributions=="NaN"]<-0
	distributions <-  array(dim = c(effects+1,choices))
	distributions[1,] <-
		exp(colSums(theta*effectContributions))/
							sum(exp(colSums(theta*effectContributions)))
	for(eff in 1:effects)
	{
		t <- theta
		t[eff] <- 0
		distributions[eff+1,] <-
			exp(colSums(t*effectContributions))/
								sum(exp(colSums(t*effectContributions)))
	}
	distributions
}

entropy <- function(distribution = NULL)
{
	entropy <- -1*(distribution %*% log(distribution)/log(length(distribution)))
	certainty <- 1-entropy
	certainty
}

KLD <- function(referenz = NULL, distributions = NULL)
{
	if(is.vector(distributions))
	{
		kld <- (referenz %*%
					(log(referenz)-log(distributions)))/log(length(referenz))
	}
	else
	{
		kld <- colSums(referenz *
					(log(referenz)-t(log(distributions))))/log(length(referenz))
	}
	kld
}

## calculates the L^1-differenz between distribution "reference"
## (which is a vector of length n)
## and each row of distributions (which is a matrix with n columns)
L1D <- function(referenz = NULL, distributions = NULL)
{
	if(is.vector(distributions))
	{
		l1d <- sum(abs(referenz-distributions))
	}
	else
	{
		l1d <- colSums(abs(referenz-t(distributions)))
	}
	l1d
}

##@print.sienaRI Methods
print.sienaRI <- function(x, ...){
	if (!inherits(x, "sienaRI"))
	{
		if (inherits(x[[1]], "sienaRI"))
		{
			cat("The components of this object ")
			cat("are Siena relative importance of effects objects.\n")
			cat("Apply the print function to the separate components.\n")
		}
		stop("not a legitimate Siena relative importance of effects object")
	}
	cat(paste("\n  Expected relative importance of effects for dependent variable '",
					x$dependentVariable,"' at observation moments:\n\n\n", sep=""))
	waves <- length(x$expectedRI)
	effs <- length(x$effectNames)
	colNames = paste("wave ", 1:waves, sep="")
	line1 <- format("", width =63)
	line2 <- paste(format(1:effs,width=3), '. ',
						format(x$effectNames, width = 56),sep="")
	line3 <- line2
	line4 <- format("  Entropy", width = 61)
	for(w in 1:length(colNames))
	{
		line1 <- paste(line1, format(colNames[w], width=8),"  ", sep = "")
		line2 <- paste(line2, format(round(x$expectedRI[[w]], 4),
									width=8, nsmall=4),"  ",sep="")
		line3 <- paste(line3, format(round(x$expectedI[[w]], 4),
									width=8, nsmall=4),"  ",sep="")
		line4 <- paste(line4, format(round(mean(x$entropyActors[[w]]), 4),
									width=8, nsmall=4),"  ",sep="")
	}
	line2 <- paste(line2, rep('\n',effs), sep="")
	line3 <- paste(line3, rep('\n',effs), sep="")
	cat(as.matrix(line1),'\n \n', sep='')
	cat(as.matrix(line2),'\n', sep='')
	cat("\n  Expected importance of effects for this dependent variable:\n\n")
	cat(as.matrix(line3),'\n\n', sep='')
	cat(as.matrix(line4),'\n', sep='')
	invisible(x)
}

##@summary.sienaRI Methods
summary.sienaRI <- function(object, ...)
{
	if (!inherits(object, "sienaRI"))
	{
		stop("not a legitimate Siena relative importance of effects object")
	}
	class(object) <- c("summary.sienaRI", class(object))
	object
}
##@print.summary.sienaRI Methods
print.summary.sienaRI <- function(x, ...)
{
	if (!inherits(x, "summary.sienaRI"))
	{
		stop("not a legitimate summary of a Siena relative importance of effects object")
	}
	print.sienaRI(x)
	invisible(x)
}


##@plot.sienaRI Methods
plot.sienaRI <- function(x, actors = NULL, col = NULL, addPieChart = FALSE,
	radius = 1, width = NULL, height = NULL, legend = TRUE,
	legendColumns = NULL, legendHeight = NULL, cex.legend = NULL,
	cex.names = NULL, ...)
{
	if (!inherits(x, "sienaRI"))
	{
		stop("not a legitimate Siena relative importance of effects object")
	}
	waves <- length(x$expectedRI)
	if (is.null(actors))
	{
		nactors <- dim(x$RIActors[[1]])[2]
		actors <- (1:nactors)
	}
	else
	{
		if ((!inherits(actors,"integer")) ||
			(min(actors) < 1) || (max(actors) > dim(x$RIActors[[1]])[2]))
		{
			stop(paste("parameter <actors> must be a set of integers from 1 to",
					dim(x$RIActors[[1]])[2]))
		}
		nactors <- length(actors)
	}
	if(legend)
	{
		if(!is.null(legendColumns))
		{
			if(is.numeric(legendColumns))
			{
				legendColumns <- as.integer(legendColumns)
			}else{
				legendColumns <- NULL
	warning("legendColumns has to be of type 'numeric' \n used default settings")
			}
		}
		if(is.null(legendColumns))
		{
			legendColumns <-floor((nactors+2)/11)
		}
		if(!is.null(legendHeight))
		{
			if(is.numeric(legendHeight))
			{
				legendHeight <- legendHeight
			}else{
				legendHeight <- NULL
	warning("legendHeight has to be of type 'numeric' \n used default settings")
			}
		}
		if(is.null(legendHeight))
		{
			legendHeight <-
				max(0.8,ceiling(length(x$effectNames)/legendColumns)*0.2)
		}
	}
	if(!is.null(height))
	{
		if(is.numeric(height))
		{
			height <- height
		}else{
			height <- NULL
			warning("height has to be of type 'numeric' \n used default settings")
		}
	}
	if(is.null(height))
	{
		height <- 1
	}

	if(!is.null(width))
	{
		if(is.numeric(width))
		{
			width <- width
		}else{
			width <- NULL
			warning("width has to be of type 'numeric' \n used default settings")
		}
	}
	if(is.null(width))
	{
		if(addPieChart)
		{
			width = (nactors/3+4)
		}else{
			width = (nactors/3+3)
		}
	}
	if(!is.null(cex.legend))
	{
		if(is.numeric(cex.legend))
		{
			cex.legend <- cex.legend
		}else{
			cex.legend <- NULL
			warning("cex.legend has to be of type 'numeric' \n used default settings")
		}
	}
	if(is.null(cex.legend))
	{
		cex.legend <- 1.3
	}

	if(!is.null(cex.names))
	{
		if(is.numeric(cex.names))
		{
			cex.names <- cex.names
		}else{
			cex.names <- NULL
			warning("cex.names has to be of type 'numeric' \n used default settings")
		}
	}
	if(is.null(cex.names))
	{
		cex.names <- 1
	}

	if(!is.null(radius))
	{
		if(is.numeric(radius))
		{
			rad <- radius
		}else{
			rad <- NULL
			warning("radius has to be of type 'numeric' \n used default settings")
		}
	}
	if(is.null(radius))
	{
		rad <- 1
	}

	if(!is.null(col))
	{
		cl <- col
	}else{
		alph <- 175
		green <- rgb(127, 201, 127,alph, maxColorValue = 255)
		lila <-rgb(190, 174, 212,alph, maxColorValue = 255)
		orange <- rgb(253, 192, 134,alph, maxColorValue = 255)
		yellow <- rgb(255, 255, 153,alph, maxColorValue = 255)
		blue <- rgb(56, 108, 176,alph, maxColorValue = 255)
		lightgray <- rgb(184,184,184,alph, maxColorValue = 255)
		darkgray <- rgb(56,56,56,alph, maxColorValue = 255)
		gray <- rgb(120,120,120,alph, maxColorValue = 255)
		pink <- rgb(240,2,127,alph, maxColorValue = 255)
		brown <- rgb(191,91,23,alph, maxColorValue = 255)
		cl <- c(green,lila,orange,yellow,blue,lightgray,darkgray,gray,pink,brown)
		while(length(cl)<length(x$effectNames)){
			alph <- (alph+75)%%255
			green <- rgb(127, 201, 127,alph, maxColorValue = 255)
			lila <-rgb(190, 174, 212,alph, maxColorValue = 255)
			orange <- rgb(253, 192, 134,alph, maxColorValue = 255)
			yellow <- rgb(255, 255, 153,alph, maxColorValue = 255)
			blue <- rgb(56, 108, 176,alph, maxColorValue = 255)
			lightgray <- rgb(184,184,184,alph, maxColorValue = 255)
			darkgray <- rgb(56,56,56,alph, maxColorValue = 255)
			gray <- rgb(120,120,120,alph, maxColorValue = 255)
			pink <- rgb(240,2,127,alph, maxColorValue = 255)
			brown <- rgb(191,91,23,alph, maxColorValue = 255)
			cl <- c(cl,green,lila,orange,yellow,blue,lightgray,darkgray,gray,pink,brown)
		}
	}
	bordergrey <-"gray25"

	if(addPieChart)
	{
		if(legend)
		{
			layoutMatrix <- matrix(c(1:(2*waves+1),(2*waves+1)), byrow= TRUE,
										ncol=2, nrow=(waves+1))
			layout(layoutMatrix,widths= c((nactors/6)+10,3.5+2.5*(rad^2)),
									heights=c(rep(height,waves),legendHeight))
		}else{
			layoutMatrix <- matrix(c(1:(2*waves)), byrow= TRUE,
												ncol=2, nrow=waves)
			layout(layoutMatrix,widths = c((nactors/6)+10,7+2.5*(rad^2)),
											heights=rep(height,waves))
		}
	}else{
		if(legend)
		{
			layoutMatrix <- matrix(c(1:(waves+1)), byrow= TRUE,
										ncol=1, nrow=(waves+1))
			layout(layoutMatrix)
		}else{
			layoutMatrix <- matrix(c(1:waves), byrow= TRUE, ncol=1, nrow=waves)
			layout(layoutMatrix, heights=2*rep(height,waves))
	# no widths, because these are only relative numbers, 
	# so requiring constant widths is redundant
		}
		par( oma = c( 1, 1, 2, 1 ), xpd=T , cex = 0.75, no.readonly = TRUE )
	}
	par(mar = c(3,3,1,1))
	for(w in 1:waves)
	{
		barplot(cbind(x$RIActors[[w]][,actors], x$expectedRI[[w]]),
			space=c(rep(0.1,nactors),1.5),width=c(rep(1,nactors),1),
			beside =FALSE, yaxt = "n", xlab="Actor", cex.names = cex.names,
			ylab=paste("wave ", w, sep=""),border=bordergrey,
			col = cl, names.arg=c(actors,"exp. rel. imp."))
		axis(2, at=c(0,0.25,0.5,0.75,1),labels=c("0","","0.5","","1"))
		axis(4, at=c(0,0.25,0.5,0.75,1),labels=c("0","","0.5","","1"))
		if(addPieChart)
		{
			pie(x$expectedRI[[w]], col = cl, labels=NA, border = bordergrey,
												radius = rad)
			mtext("exp. rel. imp.",side = 1, line = 1, cex=cex.names*0.75)
		}
	}
	if(legend)
	{
		plot(c(0,1), c(0,1), col=rgb(0,0,0,0),axes=FALSE, ylab = "", xlab = "")
		legend(0, 1, x$effectNames, fill=cl, ncol = legendColumns,
													bty = "n", cex=cex.legend)
	}
	invisible(cl)
}

