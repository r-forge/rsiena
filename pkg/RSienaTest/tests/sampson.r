sampson_t2 <- as.matrix(read.table("Sampson_t2.txt"))
sampson_t3 <- as.matrix(read.table("Sampson_t3.txt"))
sampson_t4 <- as.matrix(read.table("Sampson_t4.txt"))

n <- 25
m2 <- matrix(0, n, n)
# put edge values in desired places
m2[sampson_t2[, 1:2]] <- sampson_t2[, 3]
m3 <- matrix(0, n, n)
m3[sampson_t3[, 1:2]] <- sampson_t3[, 3]
m4 <- matrix(0, n, n)
m4[sampson_t4[, 1:2]] <- sampson_t4[, 3]

mplus2 <- m2
mplus2[mplus2 < 0] <- 0
mplus2[mplus2 > 0] <- 1
mplus3 <- m3
mplus3[mplus3 < 0] <- 0
mplus3[mplus3 > 0] <- 1
mplus4 <- m4
mplus4[mplus4 < 0] <- 0
mplus4[mplus4 > 0] <- 1

mmin2 <- m2
mmin2[mmin2 > 0] <- 0
mmin2[mmin2 < 0] <- 1
mmin3 <- m3
mmin3[mmin3 > 0] <- 0
mmin3[mmin3 < 0] <- 1
mmin4 <- m4
mmin4[mmin4 > 0] <- 0
mmin4[mmin4 < 0] <- 1

degs1 <- rowSums(mmin2+mmin3+mmin4+mplus2+mplus3+mplus4)
degs2 <- colSums(mmin2+mmin3+mmin4+mplus2+mplus3+mplus4)
# actors 1-5, 9, 14 have all degrees 0
absent <- (degs1 + degs2 == 0)


mplus2 <- mplus2[!absent,!absent]
mplus3 <- mplus3[!absent,!absent]
mplus4 <- mplus4[!absent,!absent]
mmin2 <- mmin2[!absent,!absent]
mmin3 <- mmin3[!absent,!absent]
mmin4 <- mmin4[!absent,!absent]

plus234 <- sienaNet(array(c(mplus2, mplus3, mplus4),
                               dim=c(18, 18, 3)))
min234 <- sienaNet(array(c(mmin2, mmin3, mmin4),
                               dim=c(18, 18, 3)))
sampson234 <- sienaDataCreate(plus234, min234)
sampson.eff234 <- getEffects(sampson234)
s.model234 <- sienaModelCreate(projname = "sampson234")
print01Report(sampson234, sampson.eff234, "sampson234")
plus23 <- sienaNet(array(c(mplus2, mplus3),
                               dim=c(18, 18, 2)))
min23 <- sienaNet(array(c(mmin2, mmin3),
                               dim=c(18, 18, 2)))
sampson23 <- sienaDataCreate(plus23, min23)
sampson.eff23 <- getEffects(sampson23)
s.model23 <- sienaModelCreate(projname = "sampson23")
