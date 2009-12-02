library(RSienaTest)
print(packageDescription("RSienaTest",fields="Repository/R-Forge/Revision"))

##test1
print('test1')
mynet1 <- sienaNet(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
mymodel<- model.create(findiff=TRUE, fn=simstats0c, projname='test1',
                       cond=FALSE)
ans <- siena07(mymodel, data=mydata, effects=myeff, batch=TRUE,parallelTesting=TRUE, verbose=TRUE)#,dll='../siena/src/RSiena.dll')
##test2
print('test2')
mymodel2 <- mymodel
mymodel2$cconditional <- TRUE
mymodel2$condvarno <- 1
mymodel2$projname <- 'test2'
ans <- siena07(mymodel2, data=mydata, effects=myeff, batch=TRUE,parallelTesting=TRUE, verbose=TRUE)#,dll='../siena/src/RSiena.dll')
##test3
mynet1 <- sienaNet(array(c(tmp3,tmp4),dim=c(32,32,2)))
mydata <- sienaDataCreate(mynet1)
myeff<- getEffects(mydata)
mymodel<- model.create(findiff=TRUE, fn = simstats0c, projname='test3',
                       cond=FALSE)
print('test3')
ans<- siena07(mymodel, data=mydata, effects=myeff,  batch=TRUE, parallelTesting=TRUE, verbose=TRUE)#,dll='../siena/src/RSiena.dll')
##test4
mymodel$projname <- 'test4'
mymodel$cconditional <- TRUE
mymodel$condvarno<- 1
print('test4')
ans<- siena07(mymodel, data=mydata, effects=myeff,  batch=TRUE, parallelTesting=TRUE, verbose=TRUE)#,dll='../siena/src/RSiena.dll')
##test5
mynet1 <- sienaNet(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
mymodel<- model.create(fn=simstats0c, projname='test5',
                       cond=FALSE)
print('test5')
ans <- siena07(mymodel, data=mydata, effects=myeff, batch=TRUE,parallelTesting=TRUE, verbose=TRUE)#,dll='../siena/src/RSiena.dll')
##test6
mymodel2 <- mymodel
mymodel2$cconditional <- TRUE
mymodel2$condvarno <- 1
mymodel2$projname <- 'test6'
print('test6')
ans <- siena07(mymodel2, data=mydata, effects=myeff, batch=TRUE,parallelTesting=TRUE, verbose=TRUE)#,dll='../siena/src/RSiena.dll')
##test7
mynet1 <- sienaNet(array(c(tmp3,tmp4),dim=c(32,32,2)))
mydata <- sienaDataCreate(mynet1)
myeff<- getEffects(mydata)
mymodel<- model.create(fn = simstats0c, projname='test7',
                       cond=FALSE)
print('test7')
ans<- siena07(mymodel, data=mydata, effects=myeff,  batch=TRUE, parallelTesting=TRUE, verbose=TRUE)#,dll='../siena/src/RSiena.dll')
##test8
mymodel$projname <- 'test8'
mymodel$cconditional <- TRUE
mymodel$condvarno<- 1
print('test8')
ans<- siena07(mymodel, data=mydata, effects=myeff,  batch=TRUE, parallelTesting=TRUE, verbose=TRUE)#,dll='../siena/src/RSiena.dll')
##test9

print('test9')
mynet1 <- sienaNet(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mynet2 <- sienaNet(s50a,type='behavior')
mydata <- sienaDataCreate(mynet1, mynet2)
myeff <- getEffects(mydata)
myeff$initialValue[96] <- 0.34699930338 ## siena3 starting values differ
mymodel<- model.create(findiff=FALSE, fn=simstats0c, projname='test9',
                       cond=FALSE)
ans<- siena07(mymodel, data=mydata, effects=myeff,  batch=TRUE, parallelTesting=TRUE, verbose=TRUE)
##test10
print('test10')
mymodel$projname <- 'test10'
mymodel$cconditional <- TRUE
mymodel$condvarno<- 1
ans<- siena07(mymodel, data=mydata, effects=myeff,  batch=TRUE, parallelTesting=TRUE, verbose=TRUE)
##test11
print('test11')
data501 <- sienaDataCreateFromSession("s50.csv", modelName="s50")
data501e <- sienaDataCreateFromSession("s50e.csv", modelName="s50e")
data501paj <- sienaDataCreateFromSession("s50paj.csv", modelName="s50paj")

model501 <- model.create( projname="s50",  cond=FALSE)
model501e <- model.create( projname="s50e", cond=FALSE )
model501paj <- model.create(projname="s50paj", cond=FALSE )
ans501 <- siena07(model501, data=data501$mydata, effects=data501$myeff,
                  parallelTesting=TRUE, batch=TRUE, verbose=TRUE)
ans501e <- siena07(model501e, data=data501e$mydata, effects=data501e$myeff,
                   parallelTesting=TRUE, batch=TRUE, verbose=TRUE)
ans501paj <- siena07(model501paj, data=data501paj$mydata,
                     effects=data501paj$myeff,
                  parallelTesting=TRUE, batch=TRUE, verbose=TRUE)
## compare with outputs in parallelchecked/
