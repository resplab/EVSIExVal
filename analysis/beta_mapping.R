

#####################Test bette mapper things#######################
source("include.R")
val_data <- val_data[sample(1:(dim(data_us)[1]),500,F),]

tt <- t.test(val_data$Y)
prev <- c(point=tt$estimate, upper_ci=tt$conf.int[2])
require(pROC)
tt <- roc(val_data$Y~val_data$pi, ci=T)
cs <- c(point=tt$auc, upper_ci=tt$ci[3])

val_data$lpi <- log(val_data$pi/(1-val_data$pi))
reg <- glm(Y~lpi, data=val_data, family=binomial(link="logit"))
tmp <- summary(reg)
A <- c(point=tmp$coefficients[1,1], sd=tmp$coefficients[1,2])
B <- c(point=tmp$coefficients[2,1], sd=tmp$coefficients[2,2])


z <- 0.02
samples <- gen_triplets(1000, z, prev, cs, A, B)
future_sample_sizes <- c(500,1000,2000,4000,8000,16000)
EVPI <- mean(apply(cbind(0,samples$NBm,samples$NBa),1,max))-max(0,mean(samples$NBm),mean(samples$NBa))
res_gf <- EVSI_gf(as.matrix(samples[,1:3]), z, future_sample_sizes, 100)
#res_g <- EVSI_generic(samples[,1:3],z, c(500,1000,2000,4000,8000,16000),1000)
plot(c(0,future_sample_sizes), c(0,res_gf$EVSI), ylim=c(0,res_gf$EVPI), type='l')
lines(c(0,max(future_sample_sizes)),rep(res_gf$EVPI,2),col='gray')
