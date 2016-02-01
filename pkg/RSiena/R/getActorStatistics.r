##@getActorStatistics. Use as RSiena:::getActorStatistics
getActorStatistics <- function(algorithm, data, effects)
{
	## The following initializations data, effects, and model
	## for calling "getTargets" in "siena07.setup.h"
	## is more or less copied from "getTargets" in "getTargets.r".
	## However, some modifications have been necessary to get it to work.
	f <- unpackData(data,algorithm)
	
	effects <- effects[effects$include,]
	effects$setting <- rep("", nrow(effects))
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
	ans <- .Call("getTargets", PACKAGE=pkgname, pData, pModel, myeffects, parallelrun=TRUE, returnActorStatistics=TRUE, returnStaticChangeContributions=FALSE)
	ans	
}

