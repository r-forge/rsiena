library(RSiena)

source('algorithms.r')

mynet1 <- sienaNet(array(c(s501,s502), dim=c(50, 50, 2)))
mynet2 <- sienaNet(s50a[, 1:2], type='behavior')
mydata <- sienaDataCreate(mynet1, mynet2)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, recip, include=TRUE)
myeff <- includeEffects(myeff, transTrip, inPopSqrt)
myeff <- includeEffects(myeff, altX, interaction1='mynet2')

MOMmodel <- sienaModelCreate(cond=FALSE, maxlike=FALSE, n3=100)
MLmodel <- sienaModelCreate(cond=FALSE,  maxlike=TRUE, n3=100)

ans <- siena07(MOMmodel, data=mydata, effects=myeff, useCluster=TRUE,
               initC=TRUE)

## MOM statistics / Gelman few/ few many/all
print("resp1")
resp1 <- algorithms(mydata, myeff, MOMmodel, nIter=20, numiter=20, useC=FALSE,                    nbrNodes=2)

## ML scores / Gelman few/ few many/all
print("resp2")
resp2 <- algorithms(mydata, myeff, MLmodel, nIter=20, numiter=20, useC=FALSE,
                    nbrNodes=2)

## SEM + Gelman
print("resp3")
resp3 <- algorithms(mydata, myeff, MLmodel, useC=FALSE, finalIter=50,
                    useOptim=TRUE, optimFinal=1, optimSchedule=rep(50, 10),
                    nbrNodes=2)
## EM varied sampling rate
print("resp4")
resp4 <- algorithms(mydata, myeff, MLmodel, useC=FALSE, useOptim=TRUE,
                    optimFinal=1, optimSchedule=c(10, 20, 20, 20, 50),
                    scale=1, nbrNodes=2)
## EM reuse history
print("resp5")
resp5 <- algorithms(mydata, myeff, MLmodel, useC=FALSE, useOptim=TRUE,
                    optimFinal=1, useHistory=TRUE, scale=0.75, nbrNodes=2)
## Geyer in final loop
print("resp6")
resp6 <- algorithms(mydata, myeff, MLmodel, useC=FALSE, useOptim=TRUE,
                    optimFinal=2, scale=0.5, nbrNodes=2)

## Geyer with varied theta
print("resp7")
resp7 <- algorithms(mydata, myeff, MLmodel, useC=FALSE, useOptim=TRUE,
                    optimFinal=2, variedThetas=TRUE, scale=0.5,
                    finalIter=500, nbrNodes=2)


## response surface
print("resp8")
resp8 <- algorithms(mydata, myeff, MOMmodel, nIter=100, numiter=10, useC=FALSE,
                    finalIter=100, responseSurface=TRUE, scale=1,
                    optimWeight=0.75, nbrNodes=2)

## profileLikelihoods
profileLikelihoods(list(theta=c(4.2, -2, 2, .7, -0.1, -0.03, 1, .4, -0.05)),
                   MLmodel, mydata, myeff,
                   1, nIter=10, nbrNodes=2)

##bayes - only one variable at the moment
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, recip, include=TRUE)
myeff <- includeEffects(myeff, transTrip, inPopSqrt)
print("resp9")
resp9 <- RSiena:::bayes(mydata, myeff, MLmodel, nrunMH=200)


