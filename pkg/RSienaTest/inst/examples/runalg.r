library(RSiena)

source('d:/sienasvn/siena/inst/examples/algorithms.r')

mynet1 <- sienaNet(array(c(s501,s502,s503), dim=c(50, 50, 3)))
mynet2 <- sienaNet(s50a[,1:2], type='behavior')
mydata <- sienaDataCreate(mynet1, mynet2)
mynet1 <- sienaNet(array(c(s501,s502), dim=c(50, 50,2)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, recip, include=TRUE)
myeff <- includeEffects(myeff, transTrip, inPopSqrt)
myeff <- includeEffects(myeff, altX, interaction1='mynet2')
myeff <- includeEffects(myeff, outRateInv, type="rate", include=FALSE)

mynetm1 <- sienaNet(array(c(tmp3, tmp4), dim=c(32,32,2)))
mydatam <- sienaDataCreate(mynetm1)
myeffm <- getEffects(mydatam)
MOMmodel <- sienaModelCreate(cond=FALSE, maxlike=FALSE, n3=100, seed=1)
MLmodel <- sienaModelCreate(cond=FALSE,  maxlike=TRUE, n3=100, seed=1)

ans <- siena07(MOMmodel, data=mydata, effects=myeff, useCluster=TRUE,
                verbose=TRUE)
ans50ml <- siena07(MLmodel, data=mydata, effects=myeff, useCluster=FALSE,
                verbose=TRUE)
ansm <- siena07(MOMmodel, data=mydatam, effects=myeffm, useCluster=TRUE,
               verbose=TRUE)
ansml <- siena07(MLmodel, data=mydatam, effects=myeffm, useCluster=FALSE,
                verbose=TRUE)

## MOM statistics / Gelman few/ few many/all
print("resp1")
resp1 <- algorithms(mydata, myeff, MOMmodel, nIter=20, numiter=20, useC=FALSE,
                    nbrNodes=1)
myeff <- includeEffects(myeff, outRateInv, type="rate", include=FALSE)

## ML scores / Gelman few/ few many/all
print("resp2")
resp2 <- algorithms(mydata, myeff, MLmodel, nIter=20, numiter=20, useC=FALSE,
                    nbrNodes=1)

## SEM + Gelman
print("resp3")
resp3 <- algorithms(mydata, myeff, MLmodel, useC=FALSE, finalIter=50,
                    useOptim=TRUE, optimFinal=1, optimSchedule=rep(50, 10),
                    nbrNodes=1)
## EM varied sampling rate
print("resp4")
resp4 <- algorithms(mydata, myeff, MLmodel, useC=FALSE, useOptim=TRUE,
                    optimFinal=1,
                    scale=1, nbrNodes=2)
## EM reuse history
print("resp5")
resp5 <- algorithms(mydata, myeff, MLmodel, useC=FALSE, useOptim=TRUE,
                    optimFinal=1, useHistory=TRUE, optimWeight=0.25,
                    scale=1, nbrNodes=2)

## Geyer in final loop
print("resp6")
resp6 <- algorithms(mydata, myeff, MLmodel, useC=FALSE, useOptim=TRUE,
                    optimFinal=2, scale=0.5, nbrNodes=2)

## Geyer with varied theta
print("resp7")
resp7 <- algorithms(mydata, myeff, MLmodel, useC=FALSE, useOptim=TRUE,
                    optimFinal=2, variedThetas=TRUE, scale=0.5,
                    finalIter=100, optimSchedule=c(10, 20, 20, 20, 50,50),
                    nbrNodes=2)


## response surface
print("resp8")
resp8 <- algorithms(mydata, myeff, MOMmodel, nIter=100, numiter=10, useC=FALSE,
                    finalIter=100, responseSurface=TRUE, scale=1,
                    optimWeight=0.75, nbrNodes=2)

## profileLikelihoods
profileLikelihoods(list(theta=c(5.2, 5.2,-2.5, 2, .7, -0.1, -0.03, 1, 1, .4, -0.05)),
                   MLmodel, mydata, myeff,
                   1,3, gridl=c(0.9,1.1),nIter=100, nbrNodes=2)

##bayes
print("resp9")
resp9 <- RSiena:::bayes(mydata, myeff, MLmodel, nrunMH=200, nbrNodes=1,save=FALSE)


resp1 <- algorithms(mydatam, myeffm, MOMmodel, nIter=20, numiter=20, useC=FALSE,
                    nbrNodes=1)

## ML scores / Gelman few/ few many/all
print("resp2")
resp2 <- algorithms(mydatam, myeffm, MLmodel, nIter=20, numiter=20, useC=FALSE,
                    nbrNodes=1)

## SEM + Gelman
print("resp3")
resp3 <- algorithms(mydatam, myeffm, MLmodel, useC=FALSE, finalIter=50,
                    useOptim=TRUE, optimFinal=1, optimSchedule=rep(50, 10),
                    nbrNodes=2)
## EM varied sampling rate
print("resp4")
resp4 <- algorithms(mydatam, myeffm, MLmodel, useC=FALSE, useOptim=TRUE,
                    optimFinal=1, optimSchedule=c(10, 20, 20, 20, 50),
                    scale=1, nbrNodes=2)
## EM reuse history
print("resp5")
resp5 <- algorithms(mydatam, myeffm, MLmodel, useC=FALSE, useOptim=TRUE,
                    optimFinal=1, useHistory=TRUE, scale=0.75, nbrNodes=2)
## Geyer in final loop
print("resp6")
resp6 <- algorithms(mydatam, myeffm, MLmodel, useC=FALSE, useOptim=TRUE,
                    optimFinal=2, scale=0.5, nbrNodes=2)

## Geyer with varied theta
print("resp7")
resp7 <- algorithms(mydatam, myeffm, MLmodel, useC=FALSE, useOptim=TRUE,
                    optimFinal=2, variedThetas=TRUE, scale=0.5,
                    finalIter=100, optimSchedule=c(10, 20, 20, 20, 50,50),
                    nbrNodes=2)


## response surface
print("resp8")
resp8 <- algorithms(mydatam, myeffm, MOMmodel, nIter=100, numiter=10, useC=FALSE,
                    finalIter=100, responseSurface=TRUE, scale=1,
                    optimWeight=0.75, nbrNodes=2)

## profileLikelihoods
profileLikelihoods(list(theta=c(5.2, 5.2,-2.5, 2, .7, -0.1, -0.03, 1, 1, .4, -0.05)),
                   MLmodel, mydata, myeff,
                   1,3, gridl=c(0.9,1.1),nIter=100, nbrNodes=2)

##bayes
print("resp9")
resp9 <- RSienaTest:::bayes(mydatam, myeffm, MLmodel, nrunMH=200, nbrNodes=2,save=FALSE)

resp9 <- RSienaTest:::bayes(mydata, myeff, mymodel, nrunMH=200, nbrNodes=1,save=FALSE)


acc <- resp9$MHacceptances[1:10,]
rej <- resp9$MHrejections[1:10,]


rdf <- data.frame(outcome=rep("Reject", nrow(rej)), rej)
adf <- data.frame(outcome=rep("Accept", nrow(rej)), acc)
ardf <- rbind(adf,rdf)
names(ardf) <- c("Outcome", "InsDiag", "CancDiag", "Permute", "InsPerm",
                                  "DelPerm", "InsMissing", "DelMissing")
varnames <- paste(names(ardf), sep="", collapse= " + ")
varcall <- paste("~ ", varnames,  sep="", collapse="")
artab <- xtabs(as.formula(varcall), data=ardf)

