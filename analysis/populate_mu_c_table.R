logit <- function(x) {log(x/(1-x))}
expit <- function(x) {1/(1+exp(-x))}



calc_c_from_beta <- function(a,b)
{
  res <- integrate(function(x){dbeta(x,a,b+1)*(1-pbeta(x,a+1,b))}, lower=0, upper=1)
  res$value
}





calc_c_from_beta_sim <- function(a,b, n_sim=10^6)
{
  pi <- rbeta(n_sim, a, b)
  Y <- rbinom(n_sim, 1, pi)
  require(pROC)
  roc(Y~pi)$auc
}




solve_beta_given_prev_cs <- function(prev,cs)
{
  f <- function(a)
  {
    b <- a*(1-prev)/prev
    (calc_c_from_beta(a,b)-cs)^2
  }
  
  res <- optimize(f,c(0.01,100))
  a <- res$minimum
  b <- a*(1-prev)/prev
  
  c(a=a,b=b)
}






#Maps an input [mu, upper bound] to 
solve_beta_given_mean_upper_ci <- function(mu, upper_ci)
{
  f <- function(x) {(qbeta(0.975, x, x*(1-mu)/mu)-upper_ci)^2}
  
  res <- optimize(f,c(0.00001,100000))
  
  a <- res$minimum
  b <- a*(1-mu)/mu
  
  list(a=a,b=b)
}



solve_beta_given_mean_var <- function(mu, v)
{
  a <- (mu*(1-mu)/v-1)*mu
  b <- (mu*(1-mu)/v-1)*(1-mu)
  
  list(a=a,b=b)
}















gen_triplets <- function(n, z, prev=c(point=0.5,upper_ci=0.6), cs=c(point=0.75, upper_ci=0.80), A=c(point=0, sd=0.2), B=c(point=1, sd=0.25))
{
  tmp <- solve_beta_given_mean_upper_ci(prev[1], prev[2])
  dist_prev <-list(a=tmp$a, b=tmp$b)
  
  tmp <- solve_beta_given_mean_upper_ci(cs[1], cs[2])
  dist_cs <-list(a=tmp$a, b=tmp$b)
  
  prevs <- rbeta(n, dist_prev$a, dist_prev$b)
  css <- rbeta(n, dist_cs$a, dist_cs$b)
  
  dists_risk <- t(apply(cbind(prevs,css), 1, function(x) {solve_beta_given_prev_cs(x[1],x[2])}))
  
  As <- rnorm(n,A[1],A[2])
  Bs <- rnorm(n,B[1],B[2])
  
  ses <- 1-pbeta(expit(As+Bs*logit(z)),dists_risk[,1]+1,dists_risk[,2])
  sps <- pbeta(expit(As+Bs*logit(z)),dists_risk[,1],dists_risk[,2]+1)
  
  data.frame(prev=prevs,
             se=ses,sp=sps, 
             risk_a=dists_risk[,1], 
             risk_b=dists_risk[,2], 
             A=As, 
             B=Bs,
             NBm=prevs*ses-(1-prevs)*(1-sps)*z/(1-z),
             NBa=prevs-(1-prevs)*z/(1-z))
}



#####################Test things#######################
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
