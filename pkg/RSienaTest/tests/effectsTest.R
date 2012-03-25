#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snidjers/siena
# *
# * File: tests/effectsTest.R
# *
# * Description: This module does a set of tests of the package.
# * Not run automatically, as excluded from the tar ball.
# *****************************************************************************/
## R CMD BATCH --no-save
## "--args startNo=1 endNo=65 doRun=1 useMult=1 useRSiena=1" effectsTest.R
## (without the carriage return)
## in the tests directory of your source will run tests 1 to 65 (all at present)
## Use third argument doRun=0 to source the code within R instead.
## Then can run test i
## by
## dotest(i, testdata)

## Other arguments: useMult = 1 to use multi-processors as directed by the tests:
## turn this off  to use valgrind, but remember the results will differ.
## Use useRSiena=1 for RSiena, 0 for RSienaTest

## Not sure how to do RCMD BATCH with valgrind: use
## R -d "valgrind --tool=memcheck --track-origins=yes --leak-check-full"
## --nosave < effectsTest.R 2> test.out
## with no carriage return. Some shells won't do the split redirected output:
## bash will. You will probably have to turn off multiple processes and will
## therefore get different answers!

## tests 1 to 7: normal, no missings, uncond, cond1, cond2, fd uncond,
## fd cond1 fd cond2, ml. Each is followed by some user defined
## interactions and sienaTimeTests.
## tests 8 to 14 same with bipartite nets and missings
## tests 15 to 19, symmetric networks types 2 to 6
## tests 20 to 24, symmetric networks, fd, models 2 to 6
## tests 25 to 29, multi groups no missings
## tests 30 to 36, mult groups bipartite, missings (only one dep var so fewer)
## tests 37 to 43, multiple network effects, constraints.
## tests 44, 45, bayes (normal, multigroup)
## tests 46 to 65 composition change with different options,
## uncond, cond, fd uncond, fd cond, ml
## siena08 gets tested at start,

## pick up the input arguments
xx <- commandArgs(trailing=TRUE)
eval(parse(text=xx))
if (!exists("startNo"))
{
	startNo <- 1
}
if (!exists("doRun"))
{
	doRun <- 0
}
if (!exists("useMult"))
{
	useMult <- 1
}
if (!exists("useRSiena"))
{
	useRSiena <- 1
}

## initialise R session
if (useRSiena)
{
	library(RSiena)
} else
{
	library(RSienaTest)
}
print(sessionInfo())

Sys.setenv(RSIENATESTING=TRUE)
set.seed(100)
options(digits=5)

## create a clean directory to work in
unlink("testdir", recursive=TRUE, force=TRUE)
dir.create("testdir")
setwd("testdir")

## copy over files from the examples directory of RSiena, which must be installed
datafiles <- system.file("examples", package="RSiena")
files1 <- list.files(datafiles, pattern="\\.dat$", full.names=TRUE)
file.copy(files1, ".", overwrite=TRUE)
files1 <- list.files(datafiles, pattern="\\.DAT$", full.names=TRUE)
file.copy(files1, ".", overwrite=TRUE)
files1 <- list.files(datafiles, pattern="\\.csv$", full.names=TRUE)
file.copy(files1, ".", overwrite=TRUE)
files1 <- list.files(datafiles, pattern="\\.sim$", full.names=TRUE)
file.copy(files1, ".", overwrite=TRUE)
files1 <- list.files(datafiles, pattern="\\.txt$", full.names=TRUE)
file.copy(files1, ".", overwrite=TRUE)

## create basic data object
use <- 1:20
mynet <- sienaNet(array(c(s501[use, use], s502[use, use], s503[use, use]),
						dim=c(length(use), length(use), 3)))
mybeh <- sienaNet(s50a[use, ], type='behavior')
myccov <- coCovar(rnorm(50)[use])
myvcov <- varCovar(matrix(rnorm(100)[1:(2 * length(use))], nrow=length(use)))
mydyccov <- coDyadCovar(matrix(rnorm(2500)[1:(length(use) * length(use))],
							   nrow=length(use)))
mydyvcov <- varDyadCovar(array(rnorm(5000)[1:(length(use) * length(use) * 2)],
							   dim=c(length(use), length(use), 2)))
mydata <- sienaDataCreate(mynet, alcohol=mybeh, myccov, myvcov, mydyccov,
						  mydyvcov)
myeff <- getEffects(mydata)
print01Report(mydata, myeff, "Siena01")

## crate data set to use for ml (could be different size)
use <- 1:20
mynetml <- sienaNet(array(c(s501[use, use], s502[use, use], s503[use, use]),
						  dim=c(20, 20, 3)))
mybehml <- sienaNet(s50a[use, ])
myccovml <- coCovar(rnorm(20))
myvcovml <- varCovar(matrix(rnorm(40), nrow=20))
mydyccovml <- coDyadCovar(matrix(rnorm(400), nrow=20))
mydyvcovml <- varDyadCovar(array(rnorm(800),dim=c(20, 20, 2)))
mydataml <- sienaDataCreate(mynetml,mybehml, myccovml, myvcovml, mydyccovml,
							mydyvcovml)
myeffml <- getEffects(mydataml)

## create a variety of models
mymodel <- sienaModelCreate(nsub=1, n3=50, seed=1)
mymodelcond1 <- sienaModelCreate(nsub=1, n3=50, cond=TRUE, condvar=1,
								 seed=2)
mymodelcond2 <- sienaModelCreate(nsub=1, n3=50, cond=TRUE, condvar=2,
								 seed=3)
mymodelfd <- sienaModelCreate(nsub=1, n3=50, findiff=TRUE, seed=4)
mymodelfdcond1 <- sienaModelCreate(nsub=1, n3=50, findiff=TRUE, cond=TRUE,
								   condvar=1, seed=5)
mymodelfdcond2 <- sienaModelCreate(nsub=1, n3=50, findiff=TRUE, cond=TRUE,
								   condvar=2, seed=6)
mymodelml <- sienaModelCreate(nsub=1, n3=50, maxlike=TRUE, seed=7, mult=1)
mymodelbayes <- sienaModelCreate(seed=8, mult=1)

## create objects to drive tests of user defined interactions
effectDetails1 <- data.frame(shortName=c("egoX", "cycle3"),
							 interaction1=c("myvcov", ""),
							 interaction2=rep("", 2),
							 include=rep(TRUE, 2),
							 stringsAsFactors=FALSE)
effectDetails2 <- data.frame(shortName=c("simX", "inAct", "egoX"),
							 interaction1=c("alcohol", "", "myccov"),
							 interaction2=rep("", 3),
							 include=c(FALSE, TRUE, TRUE),
							 stringsAsFactors=FALSE)
effectDetails3 <- data.frame(shortName=c("indeg", "avAlt"),
							 interaction1=c("mynet", "mynet"),
							 interaction2=rep("", 2),
							 include=c(TRUE, TRUE),
							 stringsAsFactors=FALSE)
effectDetails4 <- data.frame(shortName=c("avAlt", "outdeg", "effFrom"),
							 interaction1=c("mynet", "mynet", "myccov"),
							 interaction2=c("", "", ""),
							 include=c(TRUE, TRUE, TRUE),
							 stringsAsFactors=FALSE)
effectDetails1ml <- data.frame(shortName=c("egoX", "cycle3"),
							 interaction1=c("myvcovml", ""),
							 interaction2=rep("", 2),
							 include=rep(TRUE, 2),
							 stringsAsFactors=FALSE)
effectDetails2ml <- data.frame(shortName=c("simX", "inAct", "egoX"),
							 interaction1=c("mybehml", "", "myccovml"),
							 interaction2=rep("", 3),
							 include=c(FALSE, TRUE, TRUE),
							 stringsAsFactors=FALSE)
effectDetails3ml <- data.frame(shortName=c("indeg", "avAlt"),
							 interaction1=c("mynetml", "mynetml"),
							 interaction2=rep("", 2),
							 include=c(TRUE, TRUE),
							 stringsAsFactors=FALSE)
effectDetails4ml <- data.frame(shortName=c("avAlt", "outdeg", "effFrom"),
							 interaction1=c("mynetml", "mynetml", "myccovml"),
							 interaction2=c("", "", ""),
							 include=c(TRUE, TRUE, TRUE),
							 stringsAsFactors=FALSE)

##create data set to test undirected models
use <- 1:20
s501s <- s501[use, use]
s502s <- s502[use, use]
s503s <- s503[use, use]
s501s <- pmin(s501s + t(s501s), 1)
s502s <- pmin(s502s + t(s502s), 1)
s503s <- pmin(s503s + t(s503s), 1)
mynets <- sienaNet(array(c(s501s, s502s, s503s),dim=c(length(use), length(use),
												3)))
mydatas <- sienaDataCreate(mynets, mybeh, myccov, myvcov, mydyccov, mydyvcov)
myeffs <- getEffects(mydatas)
myeffs$initialValue[1:2] <- c(0.5, 0.5)
mymodelAFORCE <- sienaModelCreate(modelType=2, nsub=1, n3=50, seed=8)
mymodelAAGREE <- sienaModelCreate(modelType=3, nsub=1, n3=50, seed=9,
								 cond=TRUE )
mymodelBFORCE <- sienaModelCreate(modelType=4, nsub=1, n3=50, seed=10 )
mymodelBAGREE <- sienaModelCreate(modelType=5, nsub=1, n3=50, seed=11)
mymodelBJOINT <- sienaModelCreate(modelType=6, nsub=1, n3=50, seed=12)
mymodelAFORCEfd <- sienaModelCreate(modelType=2, nsub=1, n3=50, seed=8,
									findiff=TRUE)
mymodelAAGREEfd <- sienaModelCreate(modelType=3, nsub=1, n3=50, seed=9,
									cond=TRUE, findiff=TRUE)
mymodelBFORCEfd <- sienaModelCreate(modelType=4, nsub=1, n3=50, seed=10,
									findiff=TRUE)
mymodelBAGREEfd <- sienaModelCreate(modelType=5, nsub=1, n3=50, seed=11,
									findiff=TRUE)
mymodelBJOINTfd <- sienaModelCreate(modelType=6, nsub=1, n3=50, seed=12,
									findiff=TRUE)
print01Report(mydatas, myeffs, "Siena02")


##create a data set with bipartite networks and missing data
tmp3a <- tmp3[1:20, ]
tmp4a <- tmp4[1:20, ]
tmp3b <- tmp3[, 13:32 ]
tmp4b <- tmp4[, 13:32]
mynetb1 <- sienaNet(array(c(tmp3a, tmp4a), dim=c(20, 32, 2)),
					type="bipartite", nodeSet=c("senders20", "receivers32"))
mynetb2 <-  sienaNet(array(c(tmp3b, tmp4b), dim=c(32, 20, 2)),
					 type="bipartite", nodeSet=c("senders32", "receivers20"))
senders20 <- sienaNodeSet(20, "senders20")
senders32 <- sienaNodeSet(32, "senders32")
receivers20 <- sienaNodeSet(20, "receivers20")
receivers32 <- sienaNodeSet(32, "receivers32")
var1 <- rbinom(20, 1, prob=c(0.4, 0.6))
var1[c(2, 6, 7)] <- NA
bcov20 <- coCovar(var1, nodeSet="senders20")
bcov32 <- coCovar(rbinom(32, 1, prob=c(0.4, 0.6)), nodeSet="receivers32")
bdycov2032 <- coDyadCovar(matrix(rbinom(640, 1, prob=c(.7, .3)), nrow=20),
						  nodeSet=c("senders20", "receivers32"))
bdycov3220 <- coDyadCovar(matrix(rbinom(640, 1, prob=c(.7, .3)), nrow=32),
						  nodeSet=c("senders32", "receivers20"))
mydatabm <-	sienaDataCreate(mynetb1, mynetb2, bcov20, bcov32, bdycov2032,
							bdycov3220, nodeSets=list(senders20, senders32,
										receivers20, receivers32))
myeffbm <- getEffects(mydatabm)
print01Report(mydatabm, myeffbm, "Siena03")

##create datasets for groups
tmp <- sienaDataCreateFromSession("baerveldt34.csv")
mygrp1 <- tmp$mydata
mygrpeff <- tmp$myeff
print01Report(mygrp1, mygrpeff, "sienag01")
## also use duplicated bm files
mygrp2 <- sienaGroupCreate(list(mydatabm, mydatabm, mydatabm))
mygrpeff2 <- getEffects(mygrp2)
print01Report(mygrp2, mygrpeff2, "sienag03")

## set up sampson data for valued networks
source("../sampson.r")

## create dataset for composition change
## create tmp2 to use with tmp3 and tmp4
vars <- as.matrix(read.table("VARS.DAT"))
cccov1 <- coCovar(vars[, 1])
cccov2 <- coCovar(vars[, 2])
tmp2 <- as.matrix(read.table("VRND32T2.DAT"))
tmp2[tmp2==9] <- NA
tmp2[tmp2 %in% c(4,5)] <- 0
tmp2[tmp2 > 0] <- 1
mynet3 <- sienaNet(array(c(tmp2, tmp3, tmp4), dim=c(32, 32, 3)))
comp1 <- sienaCompositionChangeFromFile("vtextexoreal.dat")
comp2 <- sienaCompositionChangeFromFile("vtextexoreal.dat", option=2)
comp3 <- sienaCompositionChangeFromFile("vtextexoreal.dat", option=3)
comp4 <- sienaCompositionChangeFromFile("vtextexoreal.dat", option=4)
cccov <- rnorm(32)
cccov <- coCovar(cccov)
cvcov <- matrix(rnorm(64), nrow=32)
cvcov <- varCovar(cvcov)
ccdyad <- matrix(rbinom(32*32, 1, 0.1), nrow=32)
ccdyad <- coDyadCovar(ccdyad)
cvdyad <- array(rbinom(32*32*2, 1, 0.2), dim=c(32,32,2))
cvdyad <- varDyadCovar(cvdyad)
mydatacc1 <- sienaDataCreate(mynet3, comp1, cccov, cvcov, cccov1, cccov2, ccdyad,
							 cvdyad)
myeffcc1 <- getEffects(mydatacc1)
print01Report(mydatacc1, myeffcc1, "sienac01")
mydatacc2 <- sienaDataCreate(mynet3, comp2, cccov, cvcov, cccov1, cccov2, ccdyad,
							 cvdyad)
myeffcc2 <- getEffects(mydatacc2)
print01Report(mydatacc2, myeffcc2, "sienac02")
mydatacc3 <- sienaDataCreate(mynet3, comp3, cccov, cvcov, cccov1, cccov2, ccdyad,
							 cvdyad)
myeffcc3 <- getEffects(mydatacc3)
print01Report(mydatacc3, myeffcc3, "sienac03")
mydatacc4 <- sienaDataCreate(mynet3, comp4, cccov, cvcov, cccov1, cccov2, ccdyad,
							 cvdyad)
myeffcc4 <- getEffects(mydatacc4)
print01Report(mydatacc4, myeffcc4, "sienac04")

tmp20 <- tmp2
tmp30 <- tmp3
tmp40 <- tmp4
tmp20[is.na(tmp2)] <- 10
tmp30[is.na(tmp3)] <- 10
tmp40[is.na(tmp4)] <- 10
diag(tmp20) <- 0
diag(tmp30) <- 0
diag(tmp40) <- 0
mynet30 <- sienaNet(array(c(tmp20, tmp30, tmp40), dim=c(32, 32, 3)))

mydatacc5 <- sienaDataCreate(mynet30, cccov, cvcov, cccov1, cccov2, ccdyad,
							 cvdyad)
myeffcc5 <- getEffects(mydatacc5)
print01Report(mydatacc5, myeffcc5, "sienac05")

## get some data for siena08
n341 <- as.matrix(read.table("N34_1.DAT"))
n343 <- as.matrix(read.table("N34_3.DAT"))
n344 <- as.matrix(read.table("N34_4.DAT"))
n346 <- as.matrix(read.table("N34_6.DAT"))
hn341 <- as.matrix(read.table("HN34_1.DAT"))
hn343 <- as.matrix(read.table("HN34_3.DAT"))
hn344 <- as.matrix(read.table("HN34_4.DAT"))
hn346 <- as.matrix(read.table("HN34_6.DAT"))
net1 <- sienaNet(array(c(n341, hn341), dim=c(45, 45, 2)))
net3 <- sienaNet(array(c(n343, hn343), dim=c(37, 37, 2)))
net4 <- sienaNet(array(c(n344, hn344), dim=c(33, 33, 2)))
net6 <- sienaNet(array(c(n346, hn346), dim=c(36, 36, 2)))

data1 <- sienaDataCreate(net1)
data3 <- sienaDataCreate(net3)
data4 <- sienaDataCreate(net4)
data6 <- sienaDataCreate(net6)
eff1 <- getEffects(data1)
eff3 <- getEffects(data3)
eff4 <- getEffects(data4)
eff6 <- getEffects(data6)
print01Report(data1, eff1, "sienaMeta1")
print01Report(data3, eff3, "sienaMeta3")
print01Report(data4, eff4, "sienaMeta4")
print01Report(data6, eff6, "sienaMeta6")

## create functions to do the tests
makefilename <- function(i, type, interaction=FALSE)
{
	if (interaction)
	{
		paste("MI", formatC(i, width=2, flag="0"),
			  type, ".out", sep="")
	}
	else
	{
		paste("M", formatC(i, width=2, flag="0"),
			  type, ".out", sep="")
	}
}


effectTestFn <- function(model, type, modelNo, effects, data, bayes, useCluster,
						 ...)
{ ##if (type =="rate" || type =='eval') return()

	myfilename1 <- makefilename(modelNo, "all")
	myfilename2 <- makefilename(modelNo, "compare")
	effectNames <- effects$shortName

	use <-
		effects$type == type & (!effects$include) & (!effectNames %in%
	                 c("unspInt", "behUnspInt"))
	count <- 0
	tot <- sum(use)
	for (i in which(use))
	{
		effects$include[i] <- TRUE
		eval <- rep(FALSE, nrow(effects))
		parms <- effects$parm[i]
		if (effects$shortName[i] %in% parmtable$shortName)
		{
			parms <- c(parms,
					   parmtable$parm[match(effects$shortName[i],
											 parmtable$shortName)])
		}
		count <- count + 1
		cat(i, ":", count, "of", tot, effects$effectNames[i],
			effects$shortName[i], "\n")
		for (parm in parms)
		{
			effects$parm[i] <- parm
			if (effects$shortName[i] %in% c("density", "linear"))
			{
				## remove the eval version
				##browser()
				eval <- effects$type == "eval" &
			    effects$name == effects$name[i] &
			    effects$shortName == effects$shortName[i]
				effects$include[eval] <- FALSE
			}
			if (effects$shortName[i] != "outRateInv" && model$modelType >= 2)
			{
				model$cconditional <- TRUE
				model$condvarno <- 1
			}
			##browser()
			## multi processes are very slow on linux, so lets only do a few
			cluster = FALSE
			if (useCluster)
			{
				if ((type == "rate" && count %in% c(1, 4)) ||
					(type =="eval" && count %in% c(2, 9, 19, 29)) ||
					(type =="endow" && count %in% c(4, 14, 39)) ||
					(type =="creation" && count %in% c(7, 19, 21, 23)))
				{
					cluster=TRUE
				}
			}
			sink(myfilename1, append=TRUE)
			if (!bayes)
			{
				ans <- try(siena07(model, verbose=TRUE, effects=effects,
								   data=data, useCluster=cluster, ...))
				if (!inherits(ans, "try-error"))
				{
					sink()
					sink(myfilename2, append=TRUE)
					print(summary(ans))
				}
			}
			else
			{
				if (type == "rate")
				{
					dfra <- diag(1, sum(effects$include))
					rate <- effects$type == "rate"
					rate <- rate[effects$include]
					dfra[rate, rate] <- 0.1
					ans <- try(bayes(data, effects, model, dfra=dfra, nwarm=10,
									 nmain=10))
				}
				else
				{
					ans <- try(bayes(data, effects, model, nwarm=10, nmain=10))
				}

				if (!inherits(ans, "try-error"))
				{
					sink()
					sink(myfilename2, append=TRUE)
					print(ans$posteriorTot)
					print(ans$posteriorMII)
				}
			}
			sink()
		}
		effects$include[i] <- FALSE
		effects$include[eval] <- TRUE
		if (model$modelType >= 2)
		{
			model$cconditional <- FALSE
			model$condvarno <- 0
		}
	}
}
## user-defined interactions
interactionTestFn <- function(model, effectDetails, effects,  name, type,
							  data, modelNo, ...)
{

	myfilename1 <- makefilename(modelNo, "all", TRUE)
	myfilename2 <- makefilename(modelNo, "compare", TRUE)
	sink(myfilename1, append=TRUE)
	effects <-
		includeInteraction(effects, effectDetails$shortName,
						   name=name,
						   type=type,
						   interaction1=effectDetails$interaction1,
						   interaction2=effectDetails$interaction2,
						   character=TRUE)
	for (i in 1:nrow(effectDetails))
	{
		effects <-
			includeEffects(effects, effectDetails$shortName[i],
						   name=name,
						   type=type,
						   interaction1=effectDetails$interaction1[i],
						   interaction2=effectDetails$interaction2[i],
						   include=effectDetails$include[i],
						   character=TRUE)
	}

	print(effects,  expand=TRUE)


	ans <- try(siena07(model, verbose=TRUE, effects=effects, data=data, ...))
	if (!inherits(ans, "try-error"))
	{
		sink()
		sink(myfilename2, append=TRUE)
		print(summary(ans))
	}
	sink()

}

testFn <- function(modelNo, model, data, effects, projname, useCluster, bayes)
{
	model$projname <- projname
	print(model)
	myfilename1 <- makefilename(modelNo, "all")
	myfilename2 <- makefilename(modelNo, "compare")
	myfilename3 <- ""
	myfilename4 <- ""
	for (type in c("rate", "eval", "endow", "creation"))
	{
		print(type)
		if (model$maxlike && type %in% c("rate"))
		{
			next
		}
		print(system.time(effectTestFn(model,type, modelNo, effects=effects,
									   useCluster=useCluster, batch=TRUE,
									   data=data, bayes=bayes)))
	}

	if ("mynet" %in% names(data$depvars))
	{
		myfilename3 <- makefilename(modelNo, "all", TRUE)
		myfilename4 <- makefilename(modelNo, "compare", TRUE)
		interactionTestFn(model, effectDetails1, effects, name="mynet",
						  type="creation", data=data, batch=TRUE, modelNo)
		interactionTestFn(model, effectDetails2, effects, name="mynet",
						  type="eval",  data=data, batch=TRUE, modelNo)
		interactionTestFn(model, effectDetails3, effects, name="alcohol",
						  type="eval", data=data, batch=TRUE, modelNo)
		interactionTestFn(model, effectDetails4, effects, name="alcohol",
						  type="endow", data=data, batch=TRUE, modelNo)
		## score test, time test
		sink(myfilename1, append=TRUE)
		effectsS <- setEffect(effects, transTrip, fix=TRUE, test=TRUE)
		ans <- try(siena07(model, verbose=TRUE, data=data, effects=effects,
						   batch=TRUE))
		if (!inherits(ans, "try-error"))
		{
			sink()
			sink(myfilename2, append=TRUE)
			print(summary(ans))
		}

		sink()
		sink(myfilename1, append=TRUE)
		effectsS <- setEffect(effectsS, indeg, name='alcohol',
							  interaction1='mynet',
							  fix=TRUE, test=TRUE)
		ans <- try(siena07(model, verbose=TRUE, data=data, batch=TRUE,
						   effects=effectsS, byWave=TRUE))
		if (!inherits(ans, "try-error"))
		{
			sink()

			sink(myfilename2, append=TRUE)
			print(summary(sienaTimeTest(ans)))
		}
		sink()
		sink(myfilename1, append=TRUE)
		effectsS <- includeTimeDummy(effects, transTrip, recip, "2")
		ans <- try(siena07(model, verbose=TRUE, data=data, batch=TRUE,
						   effects=effectsS, byWave=TRUE))
		if (!inherits(ans, "try-error"))
		{
			sink()
			sink(myfilename2, append=TRUE)
			sienaTimeTest(ans)
			sink()
			sink(myfilename1, append=TRUE)
			ans1 <- try(siena07(model, verbose=TRUE, data=data, batch=TRUE,
								effects=effectsS, prevAns=ans))
			if (!inherits(ans1, "try-error"))
			{
				sink()
				sink(myfilename2, append=TRUE)
				print(summary(ans1))
			}
		}
		sink()
	}
	if ("mynetml" %in% names(data$depvars))
	{
		myfilename3 <- makefilename(modelNo, "all", TRUE)
		myfilename4 <- makefilename(modelNo, "compare", TRUE)
		interactionTestFn(model, effectDetails1ml, effects, name="mynetml",
						  type="creation", data=data, batch=TRUE, modelNo)
		interactionTestFn(model, effectDetails2ml, effects, name="mynetml",
						  type="eval",  data=data, batch=TRUE, modelNo)
		interactionTestFn(model, effectDetails3ml, effects, name="mybehml",
						  type="eval", data=data, batch=TRUE, modelNo)
		interactionTestFn(model, effectDetails4ml, effects, name="mybehml",
						  type="endow", data=data, batch=TRUE, modelNo)
		## score test, time test
		sink(myfilename1, append=TRUE)
		effectsS <- setEffect(effects, transTrip, fix=TRUE, test=TRUE)
		ans <- try(siena07(model, verbose=TRUE, data=data, effects=effects,
						   batch=TRUE))
		if (!inherits(ans, "try-error"))
		{
			sink()
			sink(myfilename2, append=TRUE)
			print(summary(ans))
		}

		sink()
		sink(myfilename1, append=TRUE)
		effectsS <- setEffect(effectsS, indeg, name='mybehml',
							  interaction1='mynetml',
							  fix=TRUE, test=TRUE)
		ans <- try(siena07(model, verbose=TRUE, data=data, batch=TRUE,
						   effects=effectsS, byWave=TRUE))
		if (!inherits(ans, "try-error"))
		{
			sink()

			sink(myfilename2, append=TRUE)
			print(summary(sienaTimeTest(ans)))
		}
		sink()
		sink(myfilename1, append=TRUE)
		effectsS <- includeTimeDummy(effects, transTrip, recip, "2")
		ans <- try(siena07(model, verbose=TRUE, data=data, batch=TRUE,
						   effects=effectsS, byWave=TRUE))
		if (!inherits(ans, "try-error"))
		{
			sink()
			sink(myfilename2, append=TRUE)
			sienaTimeTest(ans)
			sink()
			sink(myfilename1, append=TRUE)
			ans1 <- try(siena07(model, verbose=TRUE, data=data, batch=TRUE,
							effects=effectsS, prevAns=ans))
			if (!inherits(ans1, "try-error"))
			{
				sink()
				sink(myfilename2, append=TRUE)
				print(summary(ans1))
			}
		}
		sink()
	}
}

dotest <- function(i, testdata)
  {
	  model <- get(testdata[i, "model"])
	  projname <- testdata[i, "projname"]
	  useCluster <- testdata[i, "useCluster"] == " TRUE"
	  data <- get(testdata[i, "data"])
	  effects <- get(testdata[i, "effects"])
	  bayes <- testdata[i, "bayes"] == "TRUE"
	  testFn(i, model, data, effects, projname, useCluster, bayes)
  }

siena08Fn <- function(modelNo, model, data, effects)
{

	myfilename1 <- "sienaMetaAll.out"
	myfilename2 <- "sienaMetaCompare.out"
	sink(myfilename1, append=TRUE)
	effects <- includeEffects(effects, transTrip, cycle3, inPopSqrt,
							  inActSqrt, outInv)
	effects <- setEffect(effects, transTies, test=TRUE, fix=TRUE)
	effects <- setEffect(effects, recip, type='endow', test=TRUE, fix=TRUE)
	ans <- try(siena07(model, verbose=TRUE, effects=effects,
								   data=data, batch=TRUE))
	if (!inherits(ans, "try-error"))
	{
		sink()
		sink(myfilename2, append=TRUE)
		print(summary(ans))
	}
	sink()
	ans
}


## create matrix to drive tests
testdata <-	as.matrix(data.frame(model=c(rep(c("mymodel", "mymodelcond1",
                                 "mymodelcond2", "mymodelfd", "mymodelfdcond1",
                                 "mymodelfdcond2", "mymodelml"), 2),
                                 "mymodelAFORCE", "mymodelAAGREE",
                                 "mymodelBFORCE", "mymodelBAGREE",
                                 "mymodelBJOINT", "mymodelAFORCEfd",
                                 "mymodelAAGREEfd", "mymodelBFORCEfd",
                                 "mymodelBAGREEfd", "mymodelBJOINTfd"),
                                 data=c(rep("mydata", 6), "mydataml",
                                 rep("mydatabm", 7), rep("mydatas", 10)),
                                 effects=c(rep("myeff", 6), "myeffml",
                                 rep("myeffbm", 7), rep("myeffs", 10)),
                                 projname=rep(c("siena01", "siena03",
                                 "siena02"), c(7, 7, 10)),
                                 useCluster=c(rep("FALSE", 7),
                                 rep("FALSE", 7), rep(" TRUE", 10)),
								 bayes=rep("FALSE", 24)))
## add on group rows
grptestdata <- testdata[1:14, ]
grptestdata[1:7, 2] <- "mygrp1"
grptestdata[1:7, 3] <- "mygrpeff"
grptestdata[1:7, 4] <- "sienag01"
grptestdata[1:7, 5] <- "FALSE"
grptestdata[8:14, 2] <- "mygrp2"
grptestdata[8:14, 3] <- "mygrpeff2"
grptestdata[8:14, 4] <- "sienag03"
grptestdata[8:14, 5] <- "FALSE"
grptestdata <- grptestdata[-c(3,6), ]
testdata <- rbind(testdata, grptestdata)

## add on constraints data rows
stestdata <- testdata[1:7, ]
stestdata[1:6, 2] <- "sampson234"
stestdata[1:6, 3] <- "sampson.eff234"
stestdata[1:6, 4] <- "sienasampson234"
stestdata[7, 2] <- "sampson23"
stestdata[7, 3] <- "sampson.eff23"
stestdata[7, 4] <- "sampson23"
stestdata[1:7, 5] <- "FALSE"
testdata <- rbind(testdata, stestdata)

## add on bayes rows
btestdata <- testdata[14, , drop=FALSE]
btestdata[, 6] <- "TRUE"
btestdata[, 1] <- "mymodelbayes"
btestdata[, 4] <- "sienabayes"
btestdata <- btestdata[c(1,1), ]
btestdata[2, 2] <- "mygrp1"
btestdata[2, 3] <- "mygrpeff"
testdata <- rbind(testdata, btestdata)

## add on composition change rows
ctestdata <- testdata[c(1, 2, 4, 5), ]
ctestdata[, 2] <- "mydatacc1"
ctestdata[, 3] <- "myeffcc1"
ctestdata[, 4] <- "sienac01"
ctestdata <- ctestdata[rep(1:4, 5), ]
ctestdata[5:8, 2] <- "mydatacc2"
ctestdata[5:8, 3] <- "myeffcc2"
ctestdata[5:8, 4] <- "sienac02"
ctestdata[9:12, 2] <- "mydatacc3"
ctestdata[9:12, 3] <- "myeffcc3"
ctestdata[9:12, 4] <- "sienac03"
ctestdata[13:16, 2] <- "mydatacc4"
ctestdata[13:16, 3] <- "myeffcc4"
ctestdata[13:16, 4] <- "sienac04"
ctestdata[17:20, 2] <- "mydatacc5"
ctestdata[17:20, 3] <- "myeffcc5"
ctestdata[17:20, 4] <- "sienac05"
testdata <- rbind(testdata, ctestdata)

## get some data for siena08

n341 <- as.matrix(read.table("N34_1.DAT"))
n343 <- as.matrix(read.table("N34_3.DAT"))
n344 <- as.matrix(read.table("N34_4.DAT"))
n346 <- as.matrix(read.table("N34_6.DAT"))
hn341 <- as.matrix(read.table("HN34_1.DAT"))
hn343 <- as.matrix(read.table("HN34_3.DAT"))
hn344 <- as.matrix(read.table("HN34_4.DAT"))
hn346 <- as.matrix(read.table("HN34_6.DAT"))
net1 <- sienaNet(array(c(n341, hn341), dim=c(45, 45, 2)))
net3 <- sienaNet(array(c(n343, hn343), dim=c(37, 37, 2)))
net4 <- sienaNet(array(c(n344, hn344), dim=c(33, 33, 2)))
net6 <- sienaNet(array(c(n346, hn346), dim=c(36, 36, 2)))

data1 <- sienaDataCreate(net1)
data3 <- sienaDataCreate(net3)
data4 <- sienaDataCreate(net4)
data6 <- sienaDataCreate(net6)
eff1 <- getEffects(data1)
eff3 <- getEffects(data3)
eff4 <- getEffects(data4)
eff6 <- getEffects(data6)
print01Report(data1, eff1, "sienaMeta1")
print01Report(data3, eff3, "sienaMeta3")
print01Report(data4, eff4, "sienaMeta4")
print01Report(data6, eff6, "sienaMeta6")

## create table for effect parameters
parmtable <- allEffects[allEffects$parm > 0, c("shortName", "parm")]
parmtable$oldparm <- parmtable$parm
parmtable$parm[parmtable$oldparm == 5] <- 6
parmtable[parmtable$shortName =="outTrunc", "parm"] <- 4
parmtable$parm[parmtable$oldparm == 1] <- 2
parmtable$parm[parmtable$oldparm == 2] <- 1
parmtable$parm[parmtable$oldparm == 25] <- 20
if (any(parmtable$parm == parmtable$oldparm))
{
	stop("parameter table needs fixing")
}

## run tests
if (!exists("endNo"))
{
	endNo <- nrow(testdata)
}
runNumbers <- startNo : endNo
## doRun == 0 means just source things.
if (doRun == 1)
{
	## do siena08 tests
	fits <- lapply(c(1, 3, 4, 6), function(i)
			   {
				   model <- sienaModelCreate(n3=50, nsub=1, seed=19)
				   model$projname <- "sienaMeta"
				   data <- get(paste("data", i, sep=""))
				   effects <- get(paste("eff", i, sep=""))
				   siena08Fn(i, model, data, effects)
			   }
				   )
	sink("sienaMetaAll.out", append=TRUE)
	names(fits) <- paste("fit", c(1, 3, 4, 6), sep="")
	ans <- try(do.call(siena08, fits))
	if (!inherits(ans, "try-error"))
	{
		sink()
		sink("sienaMetaCompare.out", append=TRUE)
		print(summary(ans))
	}

	sink()


	lapply(runNumbers, dotest, testdata=testdata)

	## difference results

	new <- dir(".", pattern="*ompare*")
	ref <- dir("../testrefs", pattern="*ompare*", full=TRUE)
	refs <- match(new, basename(ref))

	tmp <- lapply(1:length(new), function(x)
			  {
				  cat(x,  new[x],  ref[refs[[x]]], "\n")
				  if (!is.na(refs[x]))
				  {
					  tools::Rdiff(new[x], ref[refs[x]])
				  }
			  })
}
print(dim(testdata))
print(runNumbers)
