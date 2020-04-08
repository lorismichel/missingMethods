######################### predict x1 from stationary distributions ###########################
##############################################################################################

######################### simple example: 1 torus

d <- genData(n = 800, d=2, dataset = "parabola", pNA = 0.05)
et <- extraTorus(X = rbind(X.NA.grid,d$X.NA), nb.nodes = 10)
ct <- closedTree(X = rbind(X.NA.grid,d$X.NA), depth = 4)


st_et <- stationaryDistr(et,X= d$X.NA, prob = prob, method = "eigen")
st_ct_all <- stationaryDistr(ct, d$X.NA, prob = prob, method = "eigen")

########################## ensemble example: many torus

d <- genData(n = 800, d=2, dataset = "parabola", pNA = 0.05)
B <- 10
for (b in 1:B){
  ct[[b]]<- closedTree(X = d$X.NA, depth = 5)
  et[[b]] <- extraTorus(X = d$X.NA, nb.nodes = 7)
}

st_et <- getSampleWeightsEnsemble(et, X = d$X.NA,
                                 prob = c(0.05,0.05,1-2*0.05))$pi.mat
#st_ct_all <- getSampleWeightsEnsemble(ct, X = d$X.NA,
#                                     prob = c(0.05,0.05,1-2*0.05), subset = FALSE)$pi.mat
#st_ct_leaves <- getSampleWeightsEnsemble(ct, X = d$X.NA,
#                                        prob = c(0.05,0.05,1-2*0.05), subset = TRUE)$pi.mat

X1 <- d$X.NA[,1]

ind_NA_inX1 <- which(is.na(X1))
true_X1 <- d$X[ind_NA_inX1]


dat_noNA <- data.frame(X1 = X1[-ind_NA_inX1],feat = st_et[-ind_NA_inX1,])

# comparing ranges of quantiles (unclear how to do practically, CI?)
rf <- ranger::ranger(formula = X1~., data = dat_noNA, quantreg=T, num.trees = 2000, mtry=20)
imputedX1 <- predict(rf, data= data.frame(feat=st_et[ind_NA_inX1,]), type="quantiles",
                     quantiles = seq(0.01,0.99, length.out = 100))$predictions


##############################################################################################
#################################### plotting ################################################
##############################################################################################
X2 <- d$X.NA[,2]
X2_points <- X2[ind_NA_inX1]
X2_points # these are the points where we try to impute the corresponding x1
X2_points

################ plot some stationary distributions ########################################
par(mfrow=c(3,2))

ind_1 <- which(d$X.NA[,2]==X2_points[4])
plot(st_et[ind_1,], main = paste0("StatDistr for X2 = ", round(X2_points[4],5)))
ind_2 <- which(d$X.NA[,2]==X2_points[5])
plot(st_et[ind_2,], main = paste0("StatDistr for X2 = ", round(X2_points[5],5)))

ind_3 <- which(d$X.NA[,2]==X2_points[12])
plot(st_et[ind_3,], main = paste0("StatDistr for X2 = ", round(X2_points[12],5)))
ind_4 <- which(d$X.NA[,2]==X2_points[19])
plot(st_et[ind_4,], main = paste0("StatDistr for X2 = ", round(X2_points[19],5)))

ind_5 <- which(d$X.NA[,2]==X2_points[13])
plot(st_et[ind_5,], main = paste0("StatDistr for X2 = ", round(X2_points[13],5)))
ind_6 <- which(d$X.NA[,2]==X2_points[14])
plot(st_et[ind_6,], main = paste0("StatDistr for X2 = ", round(X2_points[14],5)))



for (i in 1:nrow(imputedX1)){
#plot(d$X)
#plot(density(imputedX1[i,], bw= 0.1))
hist(imputedX1[i,], breaks=50, xlim=c(-1,1))
}




############################################################################################
############################# predict x2 from stationary distributions #####################
############################################################################################

X2 <- d$X.NA[,2]
ind_NA_inX2 <- which(is.na(X2))

dat_noNA <- data.frame(X2 = X2[-ind_NA_inX2],feat = st_et[-ind_NA_inX2,])


# version 1: comparing ranges of quantiles (unclear how to do practically, CI?)
rf <- ranger::ranger(formula = X2~., data = dat_noNA, quantreg=T)
imputedX2 <- predict(rf, data= data.frame(feat=st_et[ind_NA_inX2,]), type="quantiles",
                     quantiles = seq(0.1,0.9, length.out = 100))$predictions

# version 2: conditional mean imputation
rf <- ranger::ranger(formula = X2~., data = dat_noNA)
imputedX1 <- predict(rf, data= data.frame(feat=st_et[ind_NA_inX2,]))$predictions

plot(d$X.NA)
imputed_points <- cbind(d$X.NA[ind_NA_inX2,1],imputedX2)
points(na.omit(imputed_points), col="red")

