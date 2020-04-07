#### predict x1 from stationary distributions

st_et <- stationaryDistr(et, d$X.NA, prob = prob, method = "eigen")
X1 <- d$X.NA[,1]

ind_NA_inX1 <- which(is.na(X1))

dat_noNA <- data.frame(X1 = X1[-ind_NA_inX1],feat = st_et[-ind_NA_inX1,])

# two versions of prediction
# version 1: comparing ranges of quantiles (unclear how to do practically, CI?)
rf <- ranger::ranger(formula = X1~., data = dat_noNA, quantreg=T)
imputedX1 <- predict(rf, data= data.frame(feat=st_et[ind_NA_inX1,]), type="quantiles",
                     quantiles = seq(0.1,0.9, by=0.05))$predictions

# version 2: conditional mean imputation
rf <- ranger::ranger(formula = X1~., data = dat_noNA)
imputedX1 <- predict(rf, data= data.frame(feat=st_et[ind_NA_inX1,]))$predictions

plot(d$X.NA)
imputed_points <- cbind(imputedX1, d$X.NA[ind_NA_inX1,2])
points(na.omit(imputed_points), col="red")

#### predict x2 from stationary distributions

X2 <- d$X.NA[,2]
ind_NA_inX2 <- which(is.na(X2))

dat_noNA <- data.frame(X2 = X2[-ind_NA_inX2],feat = st_et[-ind_NA_inX2,])


# version 1: comparing ranges of quantiles (unclear how to do practically, CI?)
rf <- ranger::ranger(formula = X2~., data = dat_noNA, quantreg=T)
imputedX2 <- predict(rf, data= data.frame(feat=st_et[ind_NA_inX2,]), type="quantiles",
                     quantiles = seq(0.1,0.9, by=0.05))$predictions

# version 2: conditional mean imputation
rf <- ranger::ranger(formula = X2~., data = dat_noNA)
imputedX1 <- predict(rf, data= data.frame(feat=st_et[ind_NA_inX2,]))$predictions

plot(d$X.NA)
imputed_points <- cbind(d$X.NA[ind_NA_inX2,1],imputedX2)
points(na.omit(imputed_points), col="red")
