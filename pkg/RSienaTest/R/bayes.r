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
                  nrunMH=100, save=TRUE, nbrNodes=1, seed=1, dfra=NULL,
                  priorSigma=NULL)
{
    createStores <- function(z)
    {
        npar <- length(z$theta)
        basicRate <- z$basicRate
        numberRows <- nmain * nrunMHBatches
        z$posteriorTot <- matrix(0, nrow=z$nGroup, ncol=npar)
        z$posteriorMII <- array(0, dim=c(z$nGroup, npar, npar))
        #z$lambdadist <- matrix(NA, nrow=numberRows, ncol=sum(basicRate))
        #z$lambdas  <- matrix(NA, nrow=numberRows, ncol=sum(basicRate))
        #z$betas <- matrix(NA, nrow=numberRows, ncol=sum(!basicRate))
        z$candidates <- array(NA, dim=c(numberRows, z$nGroup, npar))
        z$acceptances <- matrix(NA, nrow=z$nGroup, ncol=numberRows)
        z$MHacceptances <- array(NA, dim=c(z$nGroup, numberRows, 7))
        z$MHrejections <- array(NA, dim=c(z$nGroup, numberRows, 7))
        z$MHproportions <- array(NA, dim=c(z$nGroup, numberRows, 7))
        z
    }
    storeData <- function()
    {
        start <- z$sub + 1
        nrun <- nrow(z$parameters)
        basicRate <- z$basicRate
        npar <- length(z$theta)
        end <- start + nrun - 1
        z$acceptances[start:end] <- z$accepts
       # z$lambdas[start:end, ] <- z$parameters[, basicRate]
       # z$lambdadist[start:end, ] <- z$shapes[, basicRate]
        z$candidates[start:end, ] <- z$parameters
      #  z$parameters[start:end, ] <- NA
      #  z$parameters[, !basicRate] <-
      #      carryForward(z$parameters[, !basicRate], z$betas[1:end,])
      #  z$betas[start:end, ] <- z$parameters[, !basicRate]
        browser()
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
            z$scaleFactors <- z$scaleFactors * ans
            if (success | iter == maxiter)
            {
                break
            }
            if (z$scaleFactors < tiny)
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
    z <- initializeBayes(data, effects, model, nbrNodes, seed, priorSigma)
    z <- createStores(z)

    z$sub <- 0

    if (is.null(dfra))
    {
        z <- getDFRA(z, 10)
    }
    else
    {
        ## maybe some validation here one day
        z$dfra <- dfra
    }
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
    z$FRAN <- NULL
    z
}

MCMCcycle <- function(z, nrunMH, nrunMHBatches)
{
    z$accepts <- matrix(NA, nrow=z$nGroup, nrunMHBatches)
    z$parameters <- array(NA, dim=c(nrunMHBatches, z$nGroup, z$pp))
    for (i in 1:nrunMHBatches)
    {
        for (j in 1:nrunMH)
        {
            ans <- z$FRAN(z, returnChains=TRUE, byGroup=TRUE)
        }
        z$chain <- ans$chain
        z <- sampleParameters(z)
        z$accepts[, i] <- z$accept
        z$parameters[i, , ] <- z$thetaMat
    }
    z$BayesAcceptances <- sum(z$accepts)
    z
}
sampleParameters <- function(z)
{
    #browser()
    f <- FRANstore()
    ## get a multivariate normal with covariance matrix dfra multiplied by a
    ## scale factor which varies between groups
    require(MASS)
    require(mvtnorm)
    thetaChanges <- t(sapply(z$scaleFactors, function(x)
                           mvrnorm(1, mu=rep(0, z$pp),
                           Sigma= x*z$dfra)))
    priorOld <- apply(z$thetaMat, 1, dmvnorm, mean=rep(0, z$pp),
                      sigma=z$priorSigma,
                      log=TRUE)
    priorNew<- apply(z$thetaMat + thetaChanges, 1, dmvnorm, mean=rep(0, z$pp),
                     sigma=z$priorSigma, log=TRUE)
    logpOld <- sapply(z$chain, function(x)
                  {
                      sum(sapply(x, function(x1) sum(x1$LogOptionSetProb +
                                     x1$LogChoiceProb)))
                  }
                      )
    thetas <- z$thetaMat + thetaChanges
    ans <- lapply(1:length(z$chain),  function(i)
              {
                  chaini <- z$chain[[i]]
                  lapply(1:length(chaini), function(i1)
                     {
                         ch <- chaini[[i1]]
                         .Call("getChainProbabilities", ch,
                               PACKAGE = pkgname,
                               f$pData, f$pModel,
                               as.integer(i),
                               as.integer(i1),
                               f$myeffects, thetas[i, ])
                     }
                         )
              }
                  )
    logpNew <- sapply(ans, function(x)
                  {
                      sum(sapply(x, function(x1) sum(x1$LogOptionSetProb +
                                                     x1$LogChoiceProb)))
                  }
                      )
    proposalProbability <- priorNew - priorOld + logpNew - logpOld
    if (log(runif(1)) < proposalProbability)
    {
        z$accept <- TRUE
        z$thetaMat <- thetas
    }
    else
    {
        z$accept <- FALSE
    }
    cat(thetas, priorNew, priorOld, logpNew, logpOld, exp(proposalProbability),
        z$accept,'\n')
    z
}


initializeBayes <- function(data, effects, model, nbrNodes, seed, priorSigma)
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
    z$FRAN <- getFromNamespace(model$FRANname, pos=grep("RSiena",
                                                   search())[1])
    z <- z$FRAN(z, model, INIT=TRUE, data=data, effects=effects)

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
    z$scaleFactors <- rep(1, z$nGroup)
    z$returnDataFrame <- TRUE # chains come back as data frames not lists
    if (is.null(priorSigma))
    {
        z$priorSigma <- diag(z$pp) * 10000
    }
    else
    {
        z$priorSigma <- priorSigma
    }
    z
}

getDFRA <- function(z, n)
{
    # do n MLmodelsteps with the initial thetas and get
    # derivs
    z$sdf <- array(0, dim=c(n, z$pp, z$pp))
    z$Phase <- 1
    z$Deriv <- TRUE
    for (i in 1:n)
    {
       ans <- z$FRAN(z)
       z$sdf[i, , ] <- ans$dff
   }
    z$dfra <-  t(apply(z$sdf, c(2, 3), mean))
    z
}

