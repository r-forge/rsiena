#******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snijders/siena/
# *
# * File: algorithms.r
# *
# * Description: Use to call simstats outside the Robbins-Monro framework.
# * Only works for one period at the moment, although the chains are returned
# * inside a group/period list structure for future expansion.
# *****************************************************************************/
#/**
##@algorithms algorithms use to call simstats outside R-M
##options useC - store change contributions in C
## useSlowC - don't store change contributions
## if neither, return change contributions to R.
pkgname <- "RSiena"
algorithms <- function(data, effects, x, ...)
{
    ## initialize
    z <- algorithmsInitialize(data, effects, x, ...)

    z$iter  <-  0

    finalLoop <- FALSE
    library(lattice)

    z$thetasub <- z$thetasub + 1
    z$thetaHistory[z$thetasub, ] <- z$theta

    repeat ## done twice once for main body and once for final loop
    {
        repeat ## done numiter times plus final loop
        {
            ########################################################
            ## initialize for a new set of samples
            ########################################################
            z$iter <- z$iter + 1
            ##print(z$theta)

            if (z$useOptim)
            {
                if (!finalLoop)
                {
                    z$nIter = z$optimSchedule[z$iter]
                    ## cat(z$iter, z$nIter, '\n')
                }
                z$importanceSamplingWeights <- rep(1, z$nIter)
            }
            cat(z$iter, z$nIter, '\n')
            #########################################################
            ## get the next set of samples
            #########################################################
            if (finalLoop && z$variedThetas)
            {
                thetas <- data.matrix(z$thetaHistory)[1:z$thetasub, ]
                if (z$thetasub > 20)
                {
                    thetas <- thetas[-c(1:10), ]
                }
                else
                {
                    thetas <- thetas[-c(1:2), ]
                }
                z <- getSamples(z, x, z$nIter, thetas=thetas)
            }
            else
            {
                z <- getSamples(z, x, z$nIter)
            }
            predictions <- colSums(colMeans(z$fra)) - z$targets

            iiter <- 0
            #######################################################
            ## do update steps with these samples
            ######################################################
            z$maxVar <- var(c(rep(0, z$nIter - 1), 1)) / 10
            repeat
            {
                iiter <- iiter + 1
                cat('iiter', iiter,'\n')
                if (z$useOptim)
                {
                    ## store old theta
                    z$oldTheta <- z$theta

                    z <- doOptimStep(z, x)
                    if (z$thetasub < nrow(z$thetaHistory))
                    {
                        z$thetasub <- z$thetasub + 1
                        z$thetaHistory[z$thetasub, ] <- z$theta
                    }

                    tmp <- getPredictions(z$theta, z, predictions=FALSE)
                    z$importanceSamplingWeights <- tmp$ps
                    varPs <-  tmp$varPs

                    if (iiter > z$maxiiter || varPs > z$maxVar ||
                        z$optimFn == 2)
                    {
                        break
                    }
                    if (all(abs(z$oldTheta - z$theta) < 0.001))
                    {
                        break
                    }
                }
                else if (z$responseSurface)
                {
                    z <- responseSurfaceChangeStep(z)
                    print(z$theta)
                    if (z$thetasub < nrow(z$thetaHistory))
                    {
                        z$thetasub <- z$thetasub + 1
                        z$thetaHistory[z$thetasub, ] <- z$theta
                    }
                    break
                }
                else
                {
                    z <- doAlgorithmChangeStep(z, x, predictions)
                    print(z$theta)
                    if (z$thetasub < nrow(z$thetaHistory))
                    {
                        z$thetasub <- z$thetasub + 1
                        z$thetaHistory[z$thetasub, ] <- z$theta
                    }
                        cat(z$oldTheta - z$theta, '\n')

                    if (finalLoop && all(abs(z$oldTheta - z$theta) < 0.01))
                    {
                        ##browser()
                        cat(z$oldTheta - z$theta, '\n')
                        break
                    }
                    else
                    {
                        tmp <- getPredictions(z$theta, z)
                        predictions <- tmp$predictions
                        varPs <- tmp$varPs
                        if (iiter > z$maxiiter || is.na(predictions))
                        {
                            break

                        }
                    }
                }
            }
            #######################################################
            ## finished with these samples
            #######################################################
            if (z$thetasub > 1)
            {
                myformula <- as.formula(paste(z$formula, z$thetasub))
                print(xyplot(myformula, data=z$thetaHistory[1:z$thetasub, ],
                             outer=TRUE,  scales="free",
                             type="l"))
            }
            ## clear the storage
            if (z$useC)
            {
                f <- RSiena:::FRANstore()
                .Call("clearStoredChains", PACKAGE=pkgname, f$pModel, x$maxlike)
            }
            if (z$useHistory) #&& z$iter > 4)
            {
                z$storedZZZ[[z$iter]] <- z$zzz
                ##     z$storedZZZ[[z$iter - 4]] <- z$zzz
                ##   if (z$iter > 14)
                if (z$iter > 3)
                {
                    z$storedZZZ[z$iter - 3] <- list(NULL)
                }
            }
            ##z$zzz <- NULL
            ##gc()
            #print(object.size(z))
            ##cat(varPs, z$maxVar,'\n')
            if (z$iter == z$numiter || finalLoop)
                ## (iiter > z$maxiiter && varPs < z$maxVar / 2) || finalLoop )
            {
                iter <- 0
                break
            }
        }  # end of inner repeat
        ## only final loop to do now
        if (finalLoop)
        {
            break
        }
        z$diag <- FALSE
        z$gain <- z$gain / 2
        z$nIter <- z$finalIter
        z$numiter <- 1
        finalLoop <- TRUE
        if (z$useOptim)
        {
            z$maxiiter <- 0
        }
        else
        {
            z$maxiiter <- 20
        }
        z$maxit <- 100
        z$optimFn <- z$optimFinal
        z$useOptim <- z$useOptim && z$useOptimFinal
    }  ## end of outer repeat loop
    if (z$thetasub > 1)
    {
        myformula <- as.formula(paste(z$formula, z$thetasub))
        print(xyplot(myformula, data=z$thetaHistory[1:z$thetasub, ],
                     outer=TRUE,  scales="free",
                     type="l"))
    }

    ans <- z$FRAN(z, x, TERM=TRUE)
    if (z$nbrNodes > 1)
    {
        stopCluster(z$cl)
    }
    if (z$useHistory)
    {
        z$storedZZZ <- NULL ## reduce size
    }
 #   z$zzz <- NULL
    z$FRAN <- NULL
    ## if FRAN is there you will load RSiena on startup
    ## if the return object is in your
    ## workspace. Then it is difficult to recreate RSiena.
    ##print(z$theta)
    z
}

getProbsAlgs <- function(x, theta, getScores, ratePar, nactors)
{
    f <- RSiena:::FRANstore()
    group <- attr(x, "group")
    period <-attr(x, "period")
    k <- ratePar[[group]][[period]]
    resp <- .Call("getChainProbabilitiesList",
                  PACKAGE = pkgname, x,
                  f$pData, f$pModel, as.integer(group),
                  as.integer(period), f$myeffects,
                  theta, getScores)
    lik <- getLikelihoodAlgs(resp[[1]])#, nactors[[group]],
                         ##theta[k])
    if (getScores)
    {
        sc <- resp[[2]]
    }
    else
    {
        sc <- NULL
    }
    list(lik=lik, sc=sc)
}

##@getProbabilities algorithms Recalculates change contributions
getProbabilities <- function(chain, theta, nactors, rateParameterPosition,
                             evalParameterPosition, endowParameterPosition,
                             getScores=FALSE, cl=NULL)
{
    f <-  RSiena:::FRANstore()

    groupPeriods <- attr(f, "groupPeriods")
    callGrid <- cbind(rep(1:f$nGroup, groupPeriods - 1),
                      as.vector(unlist(sapply(groupPeriods - 1,
                                              function(x) 1:x))))
    nIter <- length(chain) / nrow(callGrid)
    if (!is.null(cl) )
    {
        use <- 1:min(length(chain), length(cl))
        tmp <- parLapply(cl[use], chain, getProbsAlgs, theta=theta,
                         getScores=getScores,
                         ratePar=rateParameterPosition, nactors=nactors
                         )
    }
    else
    {
        tmp <- lapply(chain, getProbsAlgs, theta=theta, getScores=getScores,
                     ratePar=rateParameterPosition, nactors=nactors)
    }
    lik <- sapply(tmp, function(x)x[[1]])
    lik <- matrix(lik, ncol=nIter)
    if (getScores)
    {
        sc <- t(sapply(tmp, function(x)x[[2]]))
        dim(sc) <- c(nrow(lik), nIter, length(theta))
    }
    else
    {
        sc <- NULL
    }
    list(lik=lik, sc=sc)
}

##getProbabilitiesFromCStore algorithms Uses stored changed contributions in C
## probably does not work. will not get scores.
getProbabilitiesFromCStore <- function(chain, theta, nactors,
                                       rateParameterPosition,
                                       evalParameterPosition,
                                       endowParameterPosition)
{
    f <-  RSiena:::FRANstore()

    for (i in 1:length(chain)) ## group
    {
        nactors <- nactors[[i]]
        ratePar <- rateParameterPosition[[i]]
        for (j in 1:length(chain[[i]])) ## period
        {
            thisRatePar <- ratePar[[j]]
            rateth <- theta[thisRatePar]
            sumrates <- sum(rateth * nactors)
            rates <- log(rateth / sumrates)
            names(rates) <- names(thisRatePar)
            names(rateth) <- names(thisRatePar)
            ans <- .Call("getStoredChainProbabilities",
                         PACKAGE = pkgname,
                         f$pData, f$pModel,
                         as.integer(i),
                         as.integer(j), f$myeffects, theta)
            logOptionSetProbs <- rep(0, length(chain[[i]][[j]]))
            nc <- matrix(0, nrow=length(chain[[i]][[j]]),
                         ncol=length(rates))
            for (k in 1:length(chain[[i]][[j]]))
            {
                vars <- sapply(chain[[i]][[j]][[k]], function(x) x[[3]])
                nc[k, ] <- table(vars)
                logOptionSetProbs[k] <- sum(rates[vars])
            }
            colnames(nc) <- names(table(vars))
            tmp <- getLikelihood3(ans, logOptionSetProbs, nactors, rateth, nc)
        }
    }
    list(lik=tmp)
}

##@getProbabilitiesR algorithms Uses stored changes in R
## probably does not work. will not calculate scores.
getProbabilitiesR <- function(chain, theta, nactors, rateParameterPosition,
                     evalParameterPosition, endowParameterPosition)
{
    for (i in 1:length(chain)) #group
    {
        nactors <- nactors[[i]]
        for (j in 1:length(chain[[i]])) #period
        {
            thisRatePar <- rateParameterPosition[[i]][[j]]
            rateth <- theta[thisRatePar]
            sumrates <- sum(rateth * nactors)
            rates <- log(rateth/sumrates)
            names(rates) <- names(thisRatePar)
            mysum <- rep(0, length(chain[[i]][[j]]))
            for (kk in 1:length(chain[[i]][[j]])) # chain
            {
                thisChain <- chain[[i]][[j]][[kk]]
                chainNumerics <-
                    sapply(thisChain,function(x)c(x[[5]], x[[6]]))
                chainStrings <-
                    sapply(thisChain,function(x)c(x[[1]], x[[3]]))
                choice <- rep(0, length(thisChain))
                subs <- ifelse(chainStrings[1, ] == "Network",
                               chainNumerics[1, ] + 1,
                               chainNumerics[2, ] + 2)
                ratesubs <- match(chainStrings[2, ], names(rates))
                opts <- rates[ratesubs]
                orderedRates <- ratesubs[order(ratesubs)]
                key <- ratesubs[1]
                evalth <- theta[evalParameterPosition[[key]]]
                endowth <- theta[endowParameterPosition[[key]]]
                for (k in order(ratesubs))
                {
                    if (ratesubs[k] != key)
                    {
                        evalth <- theta[evalParameterPosition[[ratesubs[k]]]]
                        key <- ratesubs[k]
                    }
                    evalchange <- thisChain[[k]][[10]]
                    endowchange <- thisChain[[k]][[11]]
                    ps <- exp((evalchange %*% evalth)[,1] +
                              (endowchange %*% endowth)[, 1])
                    choice[k] <- log(ps[subs[k]] / sum(ps, na.rm=TRUE))
                }
                mysum[kk] <- getLikelihood2(choice, opts, nactors, rateth,
                                            table(chainStrings[2,]))
            }

        }
    }
    list(lik=mysum)
}

getLikelihoodAlgs <- function(chain)#, nactors, lambda)
{
    loglik <- 0
   # ncvals <- sapply(chain, function(x)x[[3]])
   # nc <- nactors
   # nc[] <- 0
   # ncvals <- table(ncvals)
   # nc[names(ncvals)] <- ncvals
    logChoiceProb <- sapply(chain, function(x)x[[9]])
    logOptionSetProb <- sapply(chain, function(x)x[[8]])
    loglik <- sum(logChoiceProb)  + sum(logOptionSetProb)
    #print(sum(logOptionSetProb))
    #loglik <- loglik - sum(nactors * lambda) + sum(nc * log(lambda))-
    #    sum(lfactorial(nc))
    mu <- attr(chain, "mu")
    sigma <- sqrt(attr(chain, "sigma2"))
    finalReciprocalRate <- attr(chain, "finalReciprocalRate")
    loglik <- loglik + dnorm(1, mu, sigma, log=TRUE) + log(finalReciprocalRate)
    loglik
}

getLikelihood2 <- function(logChoiceProb,logOptionSetProb, nactors, lambda,
                           nc)
{
    loglik <- 0
    nc <- nc[match(names(nactors), names(nc))]
    loglik <- sum(logChoiceProb)  #+ sum(logOptionSetProb)
  # cat (sum(logChoiceProb),'\n')
    loglik <- loglik - sum(nactors * lambda) + sum(nc * log(lambda)) -
        sum(lfactorial(nc))
    loglik
}

getLikelihood3 <- function(sumLogChoiceProb, sumlogOptionSetProb, nactors,
                           lambda, nc)
{
    loglik <- 0
    loglik <- sumLogChoiceProb  #+ sumlogOptionSetProb
    nc <- nc[, match(names(lambda), colnames(nc)), drop=FALSE]
    loglik <- loglik - sum(nactors * lambda) +
        rowSums(nc * rep(log(lambda), each=length(sumLogChoiceProb))) -
        rowSums(lfactorial(nc))
    loglik
}

doOptimStep <- function(z, x)
{
    theta <- z$theta
    theta[z$posj] <- log(theta[z$posj])
 ## browser()
 ##   ss <- diag(z$dfra)
 ##  ss[z$posj] <- ss[z$posj] * theta[z$posj]
   ## browser()
   ## tmp <- optim(theta, optimFn, gr=derivFn, z=z,  method="BFGS", control=list(maxit=1000, trace=1000,
    ##                            parscale=1/ss))
    if (z$optimFn == 1)
    {
        tmp <- optim(theta, optimFn, gr=derivFn, z=z,
                     method='BFGS',
                     control=list(maxit=z$maxit, trace=100))
    }
    else
    {
        tmp <- optim(theta, optimFn2, gr=derivFn2,
                     z=z, method='BFGS',
                     control=list(maxit=z$maxit, trace=100))
    }
  #print(tmp)
  #print(tmp)
    #browser()
    move <- sqrt(sum((theta - tmp$par)^2))
   # cat(move, diag(z$dfra), '\n')
    maxmove <- z$gain
   # browser()
    if (move > maxmove)
    {
        move <- z$gain/move
    }
    else
    {
        move <- 1
    }
    z$par <- theta -  move * (theta - tmp$par)
    z$par[z$posj] <- exp(z$par[z$posj])
    z$theta <- z$par
    #cat(z$theta, move,'\n')
    z
}

derivFn <- function(theta, z)
{
    theta[z$posj] <- exp(theta[z$posj])
    sc <- getLikelihoods(theta, z$zzz, z, getScores=TRUE)$sc
    sc <- colMeans(sc)
    sc[z$posj] <- sc[z$posj] * theta[z$posj]
    if (z$useHistory && z$iter > 1) #was 5
    {
        prevScores <- sapply(z$storedZZZ, function(x)
                         {
                             if (is.null(x))
                             {
                                 rep(0, length(z$theta))
                             }
                             else
                             {
                                 prevSc <- getLikelihoods(theta, x, z,
                                                getScores=TRUE)$sc
                                 prevSc <- colMeans(prevSc)
                                 prevSc [z$posj] <- prevSc[z$posj] *
                                     theta[z$posj]
                                 prevSc
                             }
                         }
                             )
        use <- max(1, z$iter - 3) : (z$iter - 1)
        use1 <- 1 : min(3, z$iter - 1)
        useWeights <- z$optimWeights[use1]
        sumWeights <- sum(useWeights)
        weights <- useWeights * z$optimWeight / sumWeights
        weight <- 1 - z$optimWeight
        weights <- rev(weights)
       # browser()
        sc <- weight * sc + rowSums(rep(weights, each=length(z$theta)) *
            prevScores[, use, drop=FALSE])
    }
   # print(sc)

    -sc
}
derivFn2 <- function(theta, z)
{
    opt <- optimFn2(theta,z)
    theta[z$posj] <- exp(theta[z$posj])
    tmp <- getLikelihoods(theta, z$zzz, z, getScores=TRUE)
    sc <- tmp$sc
    sc[, z$posj] <- sc[, z$posj] * rep(theta[z$posj], each=nrow(sc))
    loglik <- tmp$lik
    ps <- exp(loglik - z$lik0)
    if (z$verbose)
    {
        stem(ps/sum(ps))
    }
    - 1/exp(opt) /z$nIter* colSums(sc*ps)
}
optimFn <- function(theta, z)
{
    theta[z$posj] <- exp(theta[z$posj])
    loglik <- getLikelihoods(theta, z$zzz, z, getScores=FALSE)$lik
    loglik <- weighted.mean(loglik, z$importanceSamplingWeights)
    if (z$useHistory && z$iter > 1)
    {
        prevLoglik <- sapply(z$storedZZZ, function(x)
                         {
                             if (is.null(x))
                             {
                                 0
                             }
                             else
                             {
                                 mean(getLikelihoods(theta, x, z,
                                                     getScores=FALSE)$lik)
                             }
                         }
                             )
    #    use <- 1 : min(10, (z$iter - 5))
        use1 <- 1 : min(3, (z$iter - 1))
        use <- max(1, z$iter - 3) : (z$iter - 1)
        useWeights <- z$optimWeights[use1]
        sumWeights <- sum(useWeights)
        weights <- useWeights * z$optimWeight / sumWeights
        weight <- 1 - z$optimWeight
        weights <- rev(weights)
        #browser()
       # print(sumWeights)
       # print(sum(weights))
       # print(weights)
       # print(weight)
        loglik <- weight * loglik + sum(weights * prevLoglik[use])
    }
    - loglik
}
optimFn2 <- function(theta, z)
{
    theta[z$posj] <- exp(theta[z$posj])
    ##print(theta)
    loglik <- getLikelihoods(theta, z$zzz, z, getScores=FALSE)$lik
    ps <- exp(loglik - z$lik0)
    ##print(-log(mean(ps)))
    -  log(mean(ps))
}
getPredictions <- function(theta, z, predictions=TRUE)
{
    ps <- getLikelihoods(theta, z$zzz, z)$lik
   ##browser()
   # ps <- colSums(ps)
   #  cat(ps,"\n")
   # browser()
   # ps <- exp(unlist(ps) - z$lik0)
    ps <- exp(ps - z$lik0)
    rawps <- ps
    ps <- ps / sum(ps)
    if (z$verbose && any(!is.na(ps)))
   {
        stem(ps)
    }
    if (predictions)
    {
        if (is.na(var(ps)) || var(ps) > z$maxVar)
        {
            predictions <- NA
        }
        else
        {
            fras <- apply(z$fra, 3, function(x) colSums(x*ps))
            ## dim fras is nwaves-1 by number of parameters (vector if nwaves=2)
            ## so recreate it
            dim(fras) <- dim(z$fra)[2:3]
            predictions <- colSums(fras) - z$targets
        }
        list(predictions=predictions, varPs=var(ps), ps)
    }
    else
    {
        list(varPs=var(ps), ps=ps)
    }
}

algorithmsInitialize <-
    function(data, effects, x, useC=TRUE, useSlowC=TRUE, scale=0.2, nIter=10,
             verbose=TRUE, numiter=10, maxiiter=10, useOptim=FALSE, diag=TRUE,
             finalIter=100, optimMaxit=1, optimWeight=0.7, useHistory=FALSE,
             optimSchedule=c(rep(c(20, 50, 100), c(10, 10, 10))),
             responseSurface=FALSE, variedThetas=FALSE,
             optimFinal=2, useOptimFinal=TRUE, nbrNodes=1,
             returnDataFrame=FALSE)
{
    ##Report(openfiles=TRUE, type="n", silent=TRUE) #initialise with no file
    z  <-  NULL
    ## get the function: simstats or maxlikec
    z$FRAN <- getFromNamespace(x$FRANname, pos=grep("RSiena",
                                                   search())[1])
    z$maxlike <- x$maxlike
    x$cconditional <-  FALSE
    z$gain <- scale
    z$diag <- diag
    z$nIter <- nIter
    z$finalIter <- finalIter
    z$verbose <- verbose
    z$useOptim <- useOptim
    z$useHistory <- useHistory
    z$maxiiter <- maxiiter
    z$numiter <-  numiter
    z$print <- FALSE
    z$useC <- useC
    z$useSlowC <- useSlowC
    z$responseSurface <- responseSurface
    z$maxit <- optimMaxit
    z$optimWeight <- optimWeight
    z$optimFinal <- optimFinal
    z$optimFn <- 1
    z$useOptimFinal <- useOptimFinal
    z$thetasub <- 0
    z$variedThetas <- variedThetas
    z$nbrNodes <- nbrNodes
    z$returnDataFrame <- returnDataFrame
    if (!is.null(x$randomSeed))
    {
        set.seed(x$randomSeed)
    }
    else
    {
        if (exists(".Random.seed"))
        {
            rm(.Random.seed, pos=1)
        }
    }

    z <- RSiena:::initializeFRAN(z, x, data, effects, prevAns=NULL,
                                 initC=FALSE, profileData=FALSE,
                                 returnDeps=FALSE)


    if (z$nbrNodes > 1)
    {
        require(snow)
        require(rlecuyer)
        clusterString <- rep("localhost", nbrNodes)
        z$cl <- makeCluster(clusterString, type = "SOCK",
                            outfile = "cluster.out")
   ##     clusterCall(z$cl, function(x)print(ls(envir=.GlobalEnv)))
    ##    clusterCall(z$cl, function(x)print(search()))
     ##   clusterCall(z$cl, function(x)print(gc()))
        clusterCall(z$cl, library, pkgname, character.only = TRUE)
        clusterCall(z$cl, RSiena:::storeinFRANstore,  RSiena:::FRANstore())
        ans <- clusterCall(z$cl, RSiena:::FRANstore)
        ans <-  clusterCall(z$cl, RSiena:::initializeFRAN, z, x,
                            initC = TRUE,
                            profileData=FALSE, returnDeps=FALSE)
        ans <- clusterCall(z$cl, source,
                           "~ruth/nuffieldsiena/inst/examples/algorithms.r")
        clusterSetupRNG(z$cl,
                        seed = as.integer(runif(6, max=.Machine$integer.max)))
        if (z$maxlike)
        {
            z$int2 <- nbrNodes
        }
    }


    ## more 'parameters' expected by simstats!
    z$Deriv <- TRUE
    if (z$nbrNodes == 1)
    {
        z$cl <- NULL
    }
    z$Phase <- 1

    if (useC)
    {
        z$addChainToStore <- TRUE
        z$needChangeContributions <- TRUE
    }
    else if (useSlowC)
    {
        z$addChainToStore <- FALSE
        z$needChangeContributions <- FALSE
    }
    else
    {
        z$addChainToStore <- FALSE
        z$needChangeContributions <- TRUE
    }
    z$thetaHistory <- matrix(NA, nrow=z$numiter*z$maxiiter + 100,
                             ncol=length(z$theta))
    colnames(z$thetaHistory) <- z$effects$shortName
    z$thetaHistory <- data.frame(z$thetaHistory)
    if (z$useOptim)
    {
        z$optimSchedule <- optimSchedule
        z$numiter <- length(z$optimSchedule)
        z$optimWeights <- z$optimWeight ^ (1:z$numiter)
        if (z$useHistory)
        {
            z$storedZZZ <- vector("list", z$numiter)
        }
    }
    z$formula <- paste(names(z$thetaHistory), collapse="+")
    z$formula <- paste(z$formula, "~ 1:")
    z
}
doSimstats <- function(i, z)
{
    RSiena:::simstats0c(z, NULL, returnDeps=FALSE, returnChains=TRUE)
}

getSamples <- function(z, x, nIter, thetas=NULL)
{
    ## get a set of z$nIter samples
    sdf <- array(0, dim=c(nIter, length(z$theta), length(z$theta)))
    if (is.null(thetas))
    {
        doDeriv <- TRUE
        thetas <- matrix(z$theta, nrow=1)
    }
    else
    {
        doDeriv <- FALSE
    }
    if (z$nbrNodes > 1 && (!z$maxlike || nrow(thetas) > 1 || z$observations > 2))
    {
        zsmall <- RSiena:::makeZsmall(z)
        if (nrow(thetas) == 1) # straight repeats with same theta
        {
            if (z$maxlike)
            {
                ## not a parLapply as we only parallelize by wave and
                ## this is done in maxlikec automatically controlled by z$int2
                tmp <- lapply(1:nIter, function(i, z)
                          {
                              RSiena:::maxlikec(z, NULL, returnDeps=FALSE,
                                                returnChains=TRUE)

                          },
                              z=zsmall
                              )
                sdf <- t(sapply(tmp, function(x) x$dff))
                dim(sdf) <- c(nIter, length(z$theta), length(z$theta))
                sc <- array(0, dim=c(nIter, z$observations - 1, length(z$theta)))
            }
            else
            {
                tmp <- parLapply(z$cl, 1:nIter, doSimstats, z=zsmall)
                sc <- t(sapply(tmp, function(x)x$sc))
                nobs <- dim(tmp[[1]]$fra)[1]
                dim(sc) <- c(nIter, nobs, z$pp)
                sdf <- array(0, dim=c(nIter, length(z$theta), length(z$theta)))
            }
            fra <- t(sapply(tmp, function(x)x$fra))
            nobs <- dim(tmp[[1]]$fra)[1]
            dim(fra) <- c(nIter, nobs, z$pp)
            z$fra <- fra
            zzz <- lapply(tmp, function(x)x$chain)
            lik <- sapply(zzz, getLikelihoodLocal, theta=z$theta,
                                     z=z)
            ## dim(lik) will be number of waves - 1 by number of samples.
            dim(lik) <- c(nobs, nIter)
            z$lik0 <- colSums(lik)
        }
        else ## multiple different thetas. Keep same thetas on one process
        {
            sc <- array(0, dim=c(nIter, z$observations - 1, length(z$theta)))
            niter <- ceiling(nIter / nrow(thetas))
            oldtheta <- z$theta
                   ##alternative: turn off z$int2 and use parRapply here.
                ## need another function defined outside this one to reduce
                ##traffic
           tmp <- apply(thetas, 1, function(th, z, n)
                         {
                             z$theta <- th
                             tt <- lapply(1:n, function(i, z1)
                                {
                                    RSiena:::maxlikec(z1, NULL,
                                                      returnChains=TRUE)
                                }, z1=z)
                             tt
                         }, z=zsmall, n=niter)

            ## likelihood
            z$lik0 <-
                c(sapply(1:length(tmp), function(i, tt, th)
                     {
                         tt <- tt[[i]]
                         th <- th[i, ]
                         zzz <- lapply(tt, function(x) x$chain)
                         lik <- sapply(zzz, function(xx)
                                        getLikelihoodLocal(xx, theta=th, z=z))
                         dim(lik) <- c(z$observations - 1, niter)
                         colSums(lik)
                     }, tt=tmp, th=thetas))
            ##fra
            fra <- lapply(tmp, function(i) lapply(i, function(x)x$fra))
            fra <- do.call(c, fra)
            fra <- array(unlist(fra),
                         dim=c(z$observations - 1, z$pp, length(fra)))
            fra <- aperm(fra, c(3, 1, 2))
            ##deriv
            sdf <- array(0, dim=c(niter*nrow(thetas),
                            length(z$theta), length(z$theta)))
            for (i in 1:nrow(thetas))
            {
                tmpa <- tmp[[i]]
                tmp1 <- t(sapply(tmpa, function(x)x$dff))
                dim(tmp1) <- c(niter, z$pp, z$pp)
                sdf[((i-1) * niter + 1):(i * niter), , ] <- tmp1
            }

            ## chains
            zzz <- lapply(tmp, function(y)
                          lapply(y, function(x)x$chain))
           zzz <- do.call(c, zzz)
            z$theta <- oldtheta
        }
    }
    else
    {
        zzz <- vector("list", nIter)
        fra <- array(0, dim=c(nIter, z$observations - 1, length(z$theta)))
        sc <- array(0, dim=c(nIter, z$observations - 1, length(z$theta)))
        sdf <- array(0, dim=c(nIter, length(z$theta), length(z$theta)))
        z$startTheta <- z$theta
        thetasub <- 1
        lik0 <- matrix(NA, ncol=nIter, nrow=z$observations - 1)
        oldtheta <- z$theta
        for (i in (1:nIter))
        {
            z$theta <- thetas[thetasub, ]
            if (thetasub < nrow(thetas))
            {
                thetasub <- thetasub + 1
            }
            else
            {
                thetasub <- 1
            }
            ans <- myfran(z, x, returnChains=TRUE)
            fra[i, , ] <- ans$fra ## statistics
            if (x$maxlike)
            {
                sdf[i, , ] <- ans$dff
            }
            else
            {
                sc[i, , ] <- ans$sc
            }
            zzz[[i]] <- ans$chain
            if (nrow(thetas) > 1)
            {
                lik0[, i] <- getLikelihoodLocal(ans$chain,z, z$theta)
            }
        }
        if (nrow(thetas) == 1)
        {
           lik0[] <- sapply(zzz, getLikelihoodLocal, theta=z$theta, z=z)
      }
        z$lik0 <- colSums(lik0)
        z$theta <- oldtheta
    }
    z$zzz <- zzz
    z$fra <- fra
    z$sc <- sc
    z$sd <- apply(z$fra[, 1, ], 2, sd)
    z$sdf <- sdf
    if (doDeriv)
    {
        if (!z$maxlike)
        {
            z$dfra <- RSiena:::derivativeFromScoresAndDeviations(z$sc, z$fra)
        }
        else
        {
            z$dfra <- t(apply(z$sdf, c(2, 3), mean))
        }
        if (inherits(try(z$dinv <- solve(z$dfra)), 'try-error'))
        {
            z$dinv <- NULL
        }
    }
    else
    {
  #      z$dfra <- NULL
  #      z$dinv <- NULL
    }
    ## flatten list of chains and add group and period as attributes
    z$zzz <- flattenChains(z$zzz)
    z
}

flattenChains <- function(zz)
{
    for (i in 1:length(zz)) ## iter
    {
        for (j in 1:length(zz[[i]])) ##group
        {
            for (k in 1:length(zz[[i]][[j]])) ## period
            {
                attr(zz[[i]][[j]][[k]], "group") <- j
                attr(zz[[i]][[j]][[k]], "period") <- k
            }
        }
    }
    zz <- do.call(c, zz)
    zz <- do.call(c, zz)
    zz
}

getLikelihoodLocal <- function(chains, z, theta)
{
    sapply(1:z$nGroup, function(x, chains, periods, nactors, ratePar)
       {
           chains2 <- chains[[x]]
           nactors <- nactors[[x]]
           ratePar <- ratePar[[x]]
           sapply(1:periods[x], function(y)
              {
                  getLikelihoodAlgs(chains2[[y]])#, nactors,
                                ##z$theta[ratePar[[y]]])
              })
       }, chains=chains, periods=z$groupPeriods,
           nactors=z$nactors,
           ratePar=z$rateParameterPosition)
}
getLikelihoods <- function(newTheta, zzz, z, getScores=FALSE)
{
    if (z$useC)
    {
        ps <- getProbabilitiesFromCStore(zzz, newTheta, z$nactors,
                                z$rateParameterPosition,
                                z$evalParameterPosition,
                                z$endowParameterPosition)

    }
    else if (z$useSlowC)
    {
        ps <-  getProbabilities(zzz, newTheta, z$nactors,
                                z$rateParameterPosition,
                                z$evalParameterPosition,
                                z$endowParameterPosition,
                                getScores=getScores, z$cl)
    }
    else
    {
        ps <- getProbabilitiesR(zzz, newTheta, z$nactors,
                       z$rateParameterPosition,
                       z$evalParameterPosition,
                       z$endowParameterPosition)
    }
    ps$lik <- colSums(ps$lik)
    if (!is.null(ps$sc))
    {
        ps$sc <- colSums(ps$sc)
    }
    ps
}
##@myfran models wrapper so profiler will include this separately
myfran <- function(z, x, ...)
{
    z$FRAN(z, x, ...)
}

doAlgorithmChangeStep <- function(z, x, fra)
{
    z$oldTheta <- z$theta
    x$diag <- z$diag
    z <- RSiena:::doChangeStep(z, x, fra)
    z
}


responseSurfaceChangeStep <- function(z, eps=0.1)
{
    npar <- length(z$theta)

    ## make a grid
    mygrid <- matrix(1, ncol=npar, nrow=npar * 4)
    mygrid[cbind(seq(1, (4 * npar), 4), 1:npar)] <- 0.95
    mygrid[cbind(seq(2, (4 * npar), 4), 1:npar)] <- 1.05
    mygrid[cbind(seq(3, (4 * npar), 4), 1:npar)] <- 0.90
    mygrid[cbind(seq(4, (4 * npar), 4), 1:npar)] <- 1.10

    thetaForGrid <- z$theta

    thiseps <- ifelse(abs(z$theta) < 1e-10, eps, 0)
    myvals <- mygrid * rep(thetaForGrid + thiseps, each=nrow(mygrid))
    colnames(myvals) <- paste("theta", 1:length(z$theta), sep=".")
    oldmaxVar <- z$maxVar
    z$maxVar <- 1
    predictedStats <- t(apply(myvals, 1, function(i, z)
                            getPredictions(i, z)$predictions, z=z))
    names(predictedStats) <- paste("prediction", 1:length(z$theta), sep=".")
    mydf<- data.frame(myvals, pred=predictedStats, weights=rep(1, nrow(mygrid)))
    if (z$iter > 1)
    {
        mydf <- rbind(z$df2, mydf)
    }
    respCols <- (npar+1) : (2*npar)
    mylm <- lm(as.matrix(mydf[, respCols]) ~ . -weights ,
               data=mydf[, -respCols], weights=weights)
    coefs <- coef(mylm)[-1, ]
    newvals <- solve(t(coefs), - mylm$coef[1, ])
    move <- sqrt(sum((z$theta - newvals)^2))
    cat(move, diag(z$dfra), '\n')
    maxmove <- z$gain
    if (move > maxmove)
    {
        move <- z$gain/move
    }
    else
    {
        move <- 1
    }
    mydf$weights <- mydf$weights * z$optimWeight
    z$df2 <- mydf
    z$theta <- z$theta -  move * (z$theta - newvals)
    z
}
responseSurfaceQuadChangeStep <- function(z, eps=0.1)
{
    npar <- length(z$theta)

    ## make a grid
    mygrid <- matrix(0, ncol=npar, nrow=npar * npar * 2)
    diag(mygrid[1:npar, ]) <- 1
    diag(mygrid[npar+1:npar, ]) <- -1
    mysubs <- matrix(0, nrow=(nrow(mygrid) - 2 * npar) * 2, ncol=2)
    mysubs[, 1] <- rep((2*npar + 1):nrow(mygrid), each=2)
    myrep <- rep(rep(1:(npar - 1), (npar - 1) : 1), 2)
    mysubs[, 2] <- t(cbind(myrep, myrep + rep(sequence((npar - 1):1), 4)))
    myrep <- c(rep(c(1.2, 0.8), each=nrow(mygrid)/ 2 - npar),
               rep(c(0.8, 1.2), (nrow(mygrid)/2 - npar)/2),
               rep(c(1.2, 0.8), (nrow(mygrid)/2 - npar)/2))
    myrep <- c(rep(c(1, -1), each=nrow(mygrid)/ 2 - npar),
               rep(c(-1, 1), (nrow(mygrid)/2 - npar)/2),
               rep(c(1, -1), (nrow(mygrid)/2 - npar)/2))
    mygrid[mysubs] <- myrep
    for (i in 1:npar)
    {
        ztheta <- z$theta
        mygrid[,i] <- mygrid[,i] * max(.1, (ztheta[i]+ eps)/10)
    }
    thetaForGrid <- z$theta
    #myvals <- mygrid * rep(thetaForGrid + eps, each=nrow(mygrid))
         myvals <-  mygrid + rep(thetaForGrid, each=nrow(mygrid))
   colnames(myvals) <- paste("theta", 1:length(z$theta), sep=".")
    oldmaxVar <- z$maxVar
    z$maxVar <- 1
    predictedStats <- t(apply(myvals, 1, function(i, z)
                            getPredictions(i, z)$predictions, z=z))
    names(predictedStats) <- paste("prediction", 1:length(z$theta), sep=".")
    ssdev <- apply(predictedStats, 1, function(x) sum(x^2))
    mydf<- data.frame(myvals, pred=ssdev/1000, weights=rep(1, nrow(mygrid)))
    if (z$iter > 1)
    {
        mydf <- rbind(z$df2, mydf)
    }
   myxx <- jitter(data.matrix(mydf[, 1:npar]))
    myxxx <- poly(myxx, degree=2, raw=TRUE)
    myxmat <- as.matrix(myxxx)
    mydf <- cbind(mydf, xx=myxmat)
    names(mydf)[(npar+3):ncol(mydf)] <-
        paste('x', 1:(ncol(mydf)-npar -2),sep='')
    mylm <- lm(pred ~ . -weights, data=mydf[,-(1:npar)], weights=mydf$weights)
    coefs <- coef(mylm)[-1]
    degrees <- attributes(myxxx)$degree
    ## split the coefficients into A and b
    b <- coefs[degrees==1]
    AA <- coefs[degrees==2]
    A <- matrix(0, npar, npar)
    A[upper.tri(A, diag=TRUE)] <- AA
    A <- A + t(A)
    diag(A) <- diag(A)/2
    eigend <- eigen(A)
    if (any(eigend$values < 0))
    {
       mat1 <- eigend$vectors
       mat2 <- diag(eigend$values)
       mat2[mat2 < 0] <- 0
       aaa <- mat1 %*% mat2 %*% t(mat1)
       require(MASS)
       newvals <- -ginv(aaa) %*% b/2
    }
    else
    {
        Ainv <- solve(A)
        newvals <- -Ainv %*% b/2
    }
    newvals <- newvals[, 1]
    move <- sqrt(sum((z$theta - newvals)^2))
   maxmove <- z$gain
    if (move > maxmove)
    {
        move <- z$gain/move
    }
    else
    {
        move <- 1
    }
    mydf$weights <- mydf$weights * z$optimWeight
    z$df2 <- mydf[, 1:(2*npar)]
    z$theta <- z$theta -  move * (z$theta - newvals)
   cat(move, newvals, z$theta,'\n')
    z
}

profileLikelihoods <- function(resp, x, data, effects,
                               i, j=NULL, gridl=c(0.8, 1.2), seqlen=5,
                               maxit=2, method="BFGS", trace=0, init=TRUE,
                               nIter=100, ...)
{
    ## initialize
    if (init)
    {
        z <- algorithmsInitialize(data, effects, x, useC=FALSE, ...)
        z$theta <- resp$theta
        if (is.null(resp$zzz))
        {
            z <- getSamples(z, x, nIter)
        }
        else
        {
            z$zzz <- resp$zzz
        }
        z$importanceSamplingWeights <- rep(1, nIter)
    }
    theta <- z$theta
    theta[z$posj] <- log(theta[z$posj])
    thetaFix <- theta
    fix <- rep(FALSE, length(theta))
    grid1 <- sort(seq(thetaFix[i]*gridl[1], thetaFix[i] * gridl[2], len=seqlen))
    fix[i] <- TRUE
    if (is.null(j)) # one dimensional density
    {
        zz <- sapply(grid1, function(x)
                 {
                     thetaFix[fix] <- x
                     tmp <- optim(theta, profOptimFn, gr=profDerivFn,
                                  z=z, fix=fix, thetaFix=thetaFix,
                                  method=method,
                                  control=list(maxit=maxit, trace=trace))
                     tmp$value
                 }
                     )
        if (z$posj[i])
        {
            grid1 <- exp(grid1)
        }
        plot(zz ~ grid1, type='b')
    }
    else # contour for 2
    {
        grid2 <- sort(seq(thetaFix[j]*gridl[1], thetaFix[j] * gridl[2],
                          len=seqlen))
        fix[j] <- TRUE
        xy <- expand.grid(grid1, grid2)
        zz <- apply(xy, 1, function(x)
                {
                    thetaFix[fix] <- x
                    tmp <- optim(theta, profOptimFn, gr=profDerivFn,
                                 z=z, fix=fix, thetaFix=thetaFix,
                                 method=method,
                                 control=list(maxit=maxit, trace=trace))
                    tmp$value
                }
                    )
        zmat <- matrix(zz, nrow=length(grid1))
        if (z$posj[i])
        {
            grid1 <- exp(grid1)
        }
        if (z$posj[j])
        {
            grid2 <- exp(grid2)
        }
        contour(grid1, grid2, zmat)
        points(grid1, grid2)
    }
    if (init)
    {
        ans <- z$FRAN(z, x, TERM=TRUE)
        if (z$nbrNodes > 1)
        {
            stopCluster(z$cl)
        }
    }
    zz
}

profOptimFn <- function(theta, z, fix, thetaFix)
{
    theta[fix] <- thetaFix[fix]
    optimFn(theta, z)
}

profDerivFn <- function(theta, z, fix, thetaFix)
{
    theta[fix] <- thetaFix[fix]
    derivFn(theta, z)
}
