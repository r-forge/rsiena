##/*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: http://www.stats.ox.ac.uk/~snidjers/siena
## *
## * File: bayes.r
## *
## * Description: This file contains the code to run Bayesian simulation
## *
## ****************************************************************************/
##@bayes algorithm  fit a Bayesian model
bayes <- function(data, effects, model, nwarm=100, nmain=100, nrunMHBatches=20,
                 nrunMH=100, save=TRUE, nbrNodes=1, seed=1)
{
    createStores <- function(z)
    {
        npar <- length(z$theta)
        basicRate <- z$basicRate
        numberRows <- nmain * nrunMHBatches *
                               (z$observations - 1)
        z$posteriorTot <- rep(0, npar)
        z$posteriorMII <- matrix(0, nrow=npar, ncol=npar)
        z$lambdadist <- matrix(NA, nrow=numberRows, ncol=sum(basicRate))
        z$lambdas  <- matrix(NA, nrow=numberRows, ncol=sum(basicRate))
        z$betas <- matrix(NA, nrow=numberRows, ncol=sum(!basicRate))
        z$candidates <- matrix(NA, nrow=numberRows, ncol=sum(!basicRate))
        z$acceptances <- rep(NA, numberRows)
        z$MHacceptances <- matrix(NA, nrow=numberRows, ncol=7)
        z$MHrejections <- matrix(NA, nrow=numberRows , ncol=7)
        z$MHproportions <- matrix(NA, nrow=numberRows, ncol=7)
        z
    }
    storeData <- function()
    {
        start <- z$sub + 1
        nrun <- nrow(z$parameters)
        basicRate <- z$basicRate
        npar <- length(z$theta)
        end <- start + nrun - 1
        z$accepts <- as.logical(z$accepts)
        z$acceptances[start:end] <- z$accepts
        z$lambdas[start:end, ] <- z$parameters[, basicRate]
        z$lambdadist[start:end, ] <- z$shapes[, basicRate]
        z$candidates[start:end, ] <- z$parameters[, !basicRate]
        z$parameters[!z$accepts, !basicRate] <- NA
        z$parameters[, !basicRate] <-
            carryForward(z$parameters[, !basicRate], z$betas[1:end,])
        z$betas[start:end, ] <- z$parameters[, !basicRate]
        z$posteriorTot <- z$posteriorTot + colSums(z$parameters)
        for (i in npar)
        {
            z$posteriorMII <- z$posteriorMII +
                outer(z$parameters[i, ], z$parameters[i, ])
        }
        z$MHacceptances[start:end, ] <- z$MHaccepts
        z$MHrejections[start:end, ] <- z$MHrejects
        z$MHproportions[start:end, ] <-z$MHaccepts/ (z$MHaccepts + z$MHrejects)
        z$sub <- z$sub + nrun
        z
    }
    ## we have removed the rejected values of betas and now want to fill in
    ## the gaps by duplicating the previous value. Need to use last one from
    ## previous or original theta if we have NA in first row.
    carryForward <- function(parameters, betas)
    {
        npar <- nrow(parameters)
        nbeta <- nrow(betas)
        if (npar < nbeta)
        {
            parameters <- rbind(betas[(nbeta - npar), ], parameters)
        }
        else
        {
            parameters <- rbind(z$theta[!z$basicRate], parameters)
        }
        parameters <-
            apply(parameters, 2, function(x)
              {
                  x.pos <- which(!is.na(x))
                  if (length(x.pos) == 0 || x.pos[1] != 1)
                  {
                      x.pos <- c(1, x.pos)
                  }
                  x[rep(x.pos, c(diff(x.pos),
                                 length(x) - x.pos[length(x.pos)] + 1))]
              }
                  )
        parameters[-1, ]
    }
    improveMH <- function(z, x, tiny=1.0e-15, desired=40, maxiter=100,
                          tolerance=15)
    {
        rescaleCGD <- function(iter)
        {
            if (actual > desired)
            {
                u <-  2 - ((iter - actual) / (iter - desired))
            }
            else
            {
                u <- 1 / (2 - (actual / desired))
            }
            if (abs(actual - desired) <= tolerance)
            {
                number <<- number + 1
                if (number == 2) success <<- TRUE
            }
            else
            {
                number <<- 0
            }
            u
        }
        iter <- 0
        number <- 0
        success <- FALSE
        repeat
        {
            iter <- iter + 1
            z <- MCMCcycle(z, nrunMH=1, nrunMHBatches=100)
            actual <- z$BayesAcceptances ## acceptances
            ans <- rescaleCGD(100 * (z$observations - 1))
            z$scaleFactor <- z$scaleFactor * ans
            if (success | iter == maxiter)
            {
                break
            }
            if (z$scaleFactor < tiny)
            {
                cat('scalefactor < tiny\n')
                browser()
            }
        }
        cat('fine tuning took ', iter, ' iterations. Scalefactor:',
            z$scaleFactor, '\n')
       z
    }
    ## ################################
    ## start of function proper
    ## ################################

    ## Should we save and restore theta if we use multiple processors?
    ## Currently get separate thetas per wave (and then mix them up).
    z <- initializeBayes(data, effects, model, nbrNodes, seed)
    z <- createStores(z)

    z$sub <- 0

    z <- improveMH(z)

    if (!save)
    {
        require(lattice)
        dev.new()
        thetaplot = dev.cur()
        dev.new()
        acceptsplot = dev.cur()
    }

    for (ii in 1:nwarm)
    {
        z <- MCMCcycle(z, nrunMH=4, nrunMHBatches=20)
        numm <- z$numm
    }

    for (ii in 1:nmain)
    {
        z <- MCMCcycle(z, nrunMH=nrunMH, nrunMHBatches=nrunMHBatches)
        z <- storeData()

        numm <- z$numm

        if (ii %% 10 == 0 && !save) ## do some plots
        {
            cat('main after ii',ii,numm, '\n')
            dev.set(thetaplot)
            basicRate <- z$basicRate
            lambdas <- z$lambdas
            lambdas[lambdas == 0] <- NA
            thetadf <- data.frame(lambdas, z$betas)
            acceptsdf <- data.frame(z$MHproportions,
                                    z$acceptances)
            lambdaNames <- paste(z$effects$name[basicRate],
                                 z$effects$shortName[basicRate],
                                 z$effects$period[basicRate],
                                 z$effects$group[basicRate], sep=".")
            betaNames <- paste(z$effects$name[!basicRate],
                               z$effects$shortName[!basicRate], sep=".")
            names(thetadf) <- make.names(c(lambdaNames, betaNames),
                                         unique=TRUE)
            names(acceptsdf) <- c("InsDiag", "CancDiag", "Permute", "InsPerm",
                                  "DelPerm", "InsMissing", "DelMissing",
                                  "BayesAccepts")
            varnames <- paste(names(thetadf), sep="", collapse= " + ")
            varcall <- paste("~ ", varnames,  sep="", collapse="")
            print(histogram(as.formula(varcall), data=thetadf, scales="free",
                            outer=TRUE, breaks=NULL, type="density",
                            panel=function(x, ...)
                        {
                            panel.histogram(x, ...)
                            panel.densityplot(x, darg=list(na.rm=TRUE), ...)
                        }
                            ))
            dev.set(acceptsplot)
            varnames <- paste(names(acceptsdf), sep="", collapse= " + ")
            varcall <- paste("~ ", varnames,  sep="", collapse="")
            print(histogram(as.formula(varcall), data=acceptsdf,
                            scales=list(x="same", y="free"),
                            outer=TRUE, breaks=NULL, type="density",
                            panel=function(x, ...)
                        {
                            panel.histogram(x, ...)
                            panel.densityplot(x, darg=list(na.rm=TRUE), ...)
                        }))
        }
    }
    z
}

MCMCcycle <- function(z, nrunMH, nrunMHBatches)
{
    f <- FRANstore() ## retrieve info

    ## set up a matrix of required group/periods
    groupPeriods <- attr(f, "groupPeriods")
    callGrid <- cbind(rep(1:f$nGroup, groupPeriods - 1),
                      as.vector(unlist(sapply(groupPeriods - 1,
                                              function(x) 1:x))))

    ## z$int2 is the number of processors if iterating by period, so 1 means
    ## we are not. Can only parallelize by period at the moment.
    ## browser()
    if (nrow(callGrid) == 1)
    {
        ans <- .Call("MCMCcycle", PACKAGE=pkgname, f$pData, f$pModel,
                     f$myeffects, as.integer(1), as.integer(1),
                     z$scaleFactor, nrunMH, nrunMHBatches)
    }
    else
    {
        if (z$int2 == 1)
        {
            anss <- apply(callGrid, 1, doMCMCcycle, z$scaleFactor,
                          nrunMH, nrunMHBatches)
        }
        else
        {
            use <- 1:(min(nrow(callGrid), z$int2))
            anss <- parRapply(z$cl[use], callGrid, doMCMCcycle, z$scaleFactor,
                          nrunMH, nrunMHBatches)
        }
         ## reorganize the anss so it looks like the normal one
        ans <- NULL
        ans[[1]] <- c(sapply(anss, "[[", 1))## acceptances
        ans[[2]] <- do.call(rbind, lapply(anss, "[[", 2))
        ans[[3]] <- do.call(rbind, lapply(anss, "[[", 3))
        ans[[4]] <- do.call(rbind, lapply(anss, "[[", 4))
        ans[[5]] <- do.call(rbind, lapply(anss, "[[", 5))
         ans[[6]] <- do.call(c, lapply(anss, "[[", 6))
   }
    ## process the return values
    z$BayesAcceptances <- sum(ans[[1]])
    z$accepts <- ans[[1]]
    z$parameters <- ans[[2]]
    z$shapes <- ans[[3]]
    z$MHaccepts <- ans[[4]]
    z$MHrejects <- ans[[5]]
    z$chains <- ans[[6]]
    z
}

doMCMCcycle <- function(x, scaleFactor, nrunMH, nrunMHBatches)
{
    group <- x[1]
    period <- x[2]
    f <- FRANstore()
    .Call("MCMCcycle", PACKAGE=pkgname, f$pData, f$pModel,
                 f$myeffects, as.integer(period),
                 as.integer(group),
                 scaleFactor, nrunMH, nrunMHBatches)
}

initializeBayes <- function(data, effects, model, nbrNodes, seed)
{
    ## initialise
    set.seed(seed)
    Report(openfiles=TRUE, type="n") #initialise with no file
    z  <-  NULL
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
    }
    model$FRAN <- getFromNamespace(model$FRANname, pos=grep("RSiena",
                                                   search())[1])
    z <- model$FRAN(z, model, INIT=TRUE, data=data, effects=effects)

    is.batch(TRUE)

    WriteOutTheta(z)

    if (nbrNodes > 1)
    {
        require(snow)
        require(rlecuyer)
        clusterString <- rep("localhost", nbrNodes)
        z$cl <- makeCluster(clusterString, type = "SOCK",
                            outfile = "cluster.out")
        clusterCall(z$cl, library, pkgname, character.only = TRUE)
        clusterCall(z$cl, storeinFRANstore,  FRANstore())
        clusterCall(z$cl, FRANstore)
        clusterCall(z$cl, initializeFRAN, z, model,
                            initC = TRUE, profileData=FALSE, returnDeps=FALSE)
        clusterSetupRNG(z$cl,
                        seed = as.integer(runif(6, max=.Machine$integer.max)))
    }

    z$numm <- 20
    z$scaleFactor <- 1


    z
}
